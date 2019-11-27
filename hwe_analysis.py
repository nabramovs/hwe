import os
import sys
import progressbar
import pymongo
import csv
import time
import re
import numpy as np
from collections import OrderedDict
from common import MongoDB, get_xpos, write_table_to_csv, get_confidence_intervals
from hardy_weinberg import Hardy_Weinberg_Equilibrium_exact_test_user_Kantale as hwe #(obs_hets, obs_hom1, obs_hom2)
from gnomad_api import get_variant_details, get_variant_details_gnomad_3
from scipy.stats import fisher_exact, skew
from figures import calculate_single_pop_significant_af_thresholds_data

HWE_P_VALUE_THRESHOLD = 0.05
RARE_HET_EXCESS_MAX_AF = 0.1
ALLELE_BALANCE_MID_RATIO_THRESHOLD = 0.5

POP_SIZES = OrderedDict()
POP_SIZES['NFE'] = 64603
POP_SIZES['AMR'] = 17720
POP_SIZES['SAS'] = 15308
POP_SIZES['FIN'] = 12562
POP_SIZES['AFR'] = 12487
POP_SIZES['EAS'] = 9977
POP_SIZES['ASJ'] = 5185

POP_MAX_AN = OrderedDict()
POP_MAX_AN['NFE'] = 64603 * 2
POP_MAX_AN['AMR'] = 17720 * 2
POP_MAX_AN['SAS'] = 15308 * 2
POP_MAX_AN['FIN'] = 12562 * 2
POP_MAX_AN['AFR'] = 12487 * 2
POP_MAX_AN['EAS'] = 9977 * 2
POP_MAX_AN['ASJ'] = 5185 * 2

VALID_CSQS = set(["transcript_ablation",
				  "splice_acceptor_variant",
				  "splice_donor_variant",
				  "stop_gained",
				  "frameshift_variant",
				  "stop_lost",
				  "start_lost",
				  "transcript_amplification",
				  "inframe_insertion",
				  "inframe_deletion",
				  "missense_variant",
				  "protein_altering_variant",
				  "splice_region_variant",
				  "incomplete_terminal_codon_variant",
				  "start_retained_variant",
				  "stop_retained_variant",
				  "synonymous_variant",
			])

# FILTER:
# Minimum population AF to be analysed for deviations from HWE
MIN_POP_AF = 0.001

OUTPUT_FOLDER = './tables/'

class VariantPop():
	def __init__(self):
		self.variant_id = ''
		self.rsid = ''
		self.chrom = ''
		self.xpos = 0
		self.pos = 0
		self.csq = ''
		self.flags = []
		self.genomes_and_exomes = False
		self.gdit_gene_name = ''
		self.gene_id = ''
		self.canonical_transcript = ''

		self.pop = ''
		self.ac = 0
		self.an = 0
		self.an_ratio = 0.0
		self.af = 0.0
		self.het = 0
		self.hom = 0
		self.ref_hom = 0
		self.exp_hom = 0
		self.exp_hom_af = 0.0
		self.exp_het = 0.0
		self.inbreeding_coef = 0.0
		self.hwe_p_value = 0.0
		self.balance = ''

		self.max_pop = ''
		self.max_pop_af = 0.0
		self.all_pop_an_ratio_pass = False
		self.all_pop_negative_balance = False

	def get_dictionary(self):
		dictionary = OrderedDict()

		dictionary['variant_id'] = self.variant_id
		dictionary['rsid'] = self.rsid
		dictionary['chrom'] = self.chrom
		dictionary['xpos'] = self.xpos
		dictionary['pos'] = self.pos
		dictionary['csq'] = self.csq
		dictionary['flags'] = self.flags
		dictionary['genomes_and_exomes'] = self.genomes_and_exomes
		dictionary['gdit_gene_name'] = self.gdit_gene_name
		dictionary['gene_id'] = self.gene_id
		dictionary['canonical_transcript'] = self.canonical_transcript

		dictionary['pop'] = self.pop
		dictionary['ac'] = self.ac
		dictionary['an'] = self.an
		dictionary['an_ratio'] = self.an_ratio
		dictionary['af'] = self.af
		dictionary['het'] = self.het
		dictionary['hom'] = self.hom
		dictionary['ref_hom'] = self.ref_hom
		dictionary['exp_hom'] = self.exp_hom
		dictionary['exp_hom_af'] = self.exp_hom_af
		dictionary['exp_het'] = self.exp_het
		dictionary['inbreeding_coef'] = self.inbreeding_coef
		dictionary['hwe_p_value'] = self.hwe_p_value
		dictionary['balance'] = self.balance

		dictionary['max_pop'] = self.max_pop
		dictionary['max_pop_af'] = self.max_pop_af
		dictionary['all_pop_an_ratio_pass'] = self.all_pop_an_ratio_pass
		dictionary['all_pop_negative_balance'] = self.all_pop_negative_balance

		return dictionary


def analyse_gnomad_variants_deviations_from_hwe(db):
	db.hw.variants_hwe_pop.drop()

	variants = db.hw.gnomad_variants.find({})

	variants_buffer = []

	total_variants = variants.count()
	variant_num = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for variant in variants:
		variant_num += 1
		bar.update((variant_num + 0.0) / total_variants)

		# FILTER: exclude variants on allosomes
		if variant['chrom'] == 'X' or variant['chrom'] == 'Y':
			continue

		# FILTER: exclude non-canonical variants
		if not variant['isCanon']:
			continue

		# FILTER: exclude non-coding variants
		if variant['consequence'] not in VALID_CSQS:
			continue

		# FILTER: exclude non-pass quality variants
		filters = []
		if variant['genome']:
			filters += variant['genome']['filters']
		if variant['exome']:
			filters += variant['exome']['filters']
		if len(filters) > 0:
			continue

		pop_acs = {}
		pop_ans = {}
		pop_homs = {}
		genomes_and_exomes = False
		if variant['genome']:
			for pop_stats in variant['genome']['populations']:
				pop = pop_stats['id']
				pop_acs[pop] = pop_stats['ac']
				pop_ans[pop] = pop_stats['an']
				pop_homs[pop] = pop_stats['ac_hom']

		if variant['exome'] and variant['genome']:
			genomes_and_exomes = True
			for pop_stats in variant['exome']['populations']:
				pop = pop_stats['id']
				pop_acs[pop] += pop_stats['ac']
				pop_ans[pop] += pop_stats['an']
				pop_homs[pop] += pop_stats['ac_hom']

		elif variant['exome']:
			for pop_stats in variant['exome']['populations']:
				pop = pop_stats['id']
				pop_acs[pop] = pop_stats['ac']
				pop_ans[pop] = pop_stats['an']
				pop_homs[pop] = pop_stats['ac_hom']

		max_pop = ''
		max_pop_af = 0.0
		all_pop_an_ratio_pass = True
		all_pop_negative_balance = True

		variant_pops = []
		for pop, max_an in POP_MAX_AN.iteritems():
			ac = pop_acs[pop]
			an = pop_ans[pop]

			# FILTER: exclude variants with no alleles
			if an == 0:
				all_pop_an_ratio_pass = False
				continue

			hom = pop_homs[pop]
			het = ac - 2 * hom
			af = float(ac) / an

			if af > max_pop_af:
				max_pop = pop
				max_pop_af = af

			# FILTER: exclude extremely rare variants (ON POP LVL)
			if af < MIN_POP_AF:
				continue

			variant_pop = VariantPop()

			variant_pop.variant_id = variant['variant_id']
			variant_pop.rsid = variant['rsid']
			variant_pop.chrom = variant['chrom']
			variant_pop.xpos = variant['xpos']
			variant_pop.pos = variant['pos']
			variant_pop.csq = variant['consequence']
			variant_pop.flags = variant['flags']
			variant_pop.genomes_and_exomes = genomes_and_exomes
			variant_pop.gdit_gene_name = variant['gdit_gene_name']
			variant_pop.gene_id = variant['gene_id']
			variant_pop.canonical_transcript = variant['canonical_transcript']

			variant_pop.pop = pop
			variant_pop.ac = ac
			variant_pop.an = an
			variant_pop.an_ratio = float(an) / max_an

			if variant_pop.an_ratio < 0.8:
				all_pop_an_ratio_pass = False

			variant_pop.af = af
			variant_pop.het = het
			variant_pop.hom = hom
			variant_pop.ref_hom = an / 2 - het - hom
			variant_pop.exp_hom = af * af * an / 2 # q^2 * indiviudlas number (an / 2)
			variant_pop.exp_hom_af = af * af

			variant_pop.exp_het = 2 * (1 - af) * af * an / 2 # 2pq * indiviudlas number (an / 2); (p = 1 - q)

			if af == 1:
				variant_pop.inbreeding_coef = 1
			else:
				variant_pop.inbreeding_coef = 1 - float(variant_pop.het) / variant_pop.exp_het

			variant_pop.hwe_p_value = hwe(variant_pop.het, variant_pop.hom, variant_pop.ref_hom)

			if variant_pop.hom < variant_pop.exp_hom:
				variant_pop.balance = '-'
			else:
				variant_pop.balance = '+'
				all_pop_negative_balance = False

			variant_pops.append(variant_pop)

		for variant_pop in variant_pops:
			variant_pop.max_pop = max_pop
			variant_pop.max_pop_af = max_pop_af
			variant_pop.all_pop_an_ratio_pass = all_pop_an_ratio_pass
			variant_pop.all_pop_negative_balance = all_pop_negative_balance

			variants_buffer.append(variant_pop.get_dictionary())

		if len(variants_buffer) > 10000:
			db.hw.variants_hwe_pop.insert_many(variants_buffer) # variants_buffer
			variants_buffer = []

	bar.finish()

	if len(variants_buffer) > 0:
		db.hw.variants_hwe_pop.insert_many(variants_buffer)

	db.hw.variants_hwe_pop.create_index([('variant_id', pymongo.ASCENDING)], name='variant_id_1')
	db.hw.variants_hwe_pop.create_index([('pop', pymongo.ASCENDING)], name='pop_1')
	db.hw.variants_hwe_pop.create_index([('xpos', pymongo.ASCENDING)], name='xpos_1')


def get_variant_alt_af(db, variant_id, xpos, pop):
	gnomad_variants = db.hw.gnomad_variants.find({'xpos': xpos})
	alt_af = 0.0
	alt_var_ids = []
	for gnomad_variant in gnomad_variants:
		if gnomad_variant['variant_id'] == variant_id:
			continue

		pop_ac = 0
		pop_an = 0

		if gnomad_variant['exome']:
			populations = gnomad_variant['exome']['populations']
			for population in populations:
				if population['id'] == pop:
					pop_ac += population['ac']
					pop_an += population['an']
					break

		if gnomad_variant['genome']:
			populations = gnomad_variant['genome']['populations']
			for population in populations:
				if population['id'] == pop:
					pop_ac += population['ac']
					pop_an += population['an']
					break

		var_af = 0.0
		if pop_an > 0 and pop_ac > 0:
			var_af = float(pop_ac) / pop_an
			alt_var_ids.append(gnomad_variant['variant_id'])

		alt_af += var_af

	return alt_af, alt_var_ids


class AltVar():
	def __init__(self):
		self.alt_af = False
		self.alt_var_ids = False

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['alt_af'] = self.alt_af
		dictionary['alt_var_ids'] = self.alt_var_ids
		return dictionary


def add_alt_af_data(db):
	variants = db.hw.variants_hwe_pop.find({})

	alt_vars = {}
	total_lines = variants.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for variant in variants:
		var_id = variant['variant_id']
		xpos = variant['xpos']
		pop = variant['pop']
		alt_var = AltVar()
		alt_var.alt_af, alt_var.alt_var_ids = get_variant_alt_af(db, var_id, xpos, pop)
		alt_vars[var_id] = alt_var

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	bulk = db.hw.variants_hwe_pop.initialize_unordered_bulk_op()
	counter = 0

	total_lines = len(alt_vars)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for var_id, alt_data in alt_vars.iteritems():

		# process in bulk
		bulk.find({ 'variant_id': var_id }).update({ '$set': { 'alt_af': alt_data.alt_af, 'alt_vars': alt_data.alt_var_ids } })
		counter += 1

		if (counter % 500 == 0):
			bulk.execute()
			bulk = db.hw.variants_hwe_pop.initialize_ordered_bulk_op()

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	if (counter % 500 != 0):
		bulk.execute()


###################################################
### HWE POP VARIANTS SEG DUP and TANDEM REPEATS ###
###################################################

def create_hwe_variants_regions(db):
	db.hw.variants_hwe_regions.drop()

	xposes = set([])
	variants = db.hw.variants_hwe_pop.find({})
	for variant in variants:
		xposes.add(variant['xpos'])

	total_lines = len(xposes)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for xpos in xposes:
		xpos_stats = OrderedDict()
		xpos_stats['_id'] = xpos
		seg_dup = db.hw.seg_dups.find_one({ "xstart": { "$lte": xpos }, "xend": { "$gte": xpos } })
		if seg_dup:
			xpos_stats['seg_dup'] = True
		else:
			xpos_stats['seg_dup'] = False

		tandem_repeat = db.hw.tandem_repeats.find_one({"_id": xpos})

		if tandem_repeat:
			xpos_stats['tandem_repeat'] = True
		else:
			xpos_stats['tandem_repeat'] = False

		db.hw.variants_hwe_regions.insert(xpos_stats)

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()


def update_variants_hwe_pop_with_region_data(db):
	xpos_seg_dups = {}
	xpos_tandem_repeats = {}
	xposes = []
	regions = db.hw.variants_hwe_regions.find({})
	for region in regions:
		xpos = region['_id']
		xposes.append(xpos) 
		xpos_seg_dups[xpos] = region['seg_dup']
		xpos_tandem_repeats[xpos] = region['tandem_repeat']

	
	bulk = db.hw.variants_hwe_pop.initialize_unordered_bulk_op()
	counter = 0

	total_lines = len(xposes)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for xpos in xposes:
		# process in bulk
		bulk.find({ 'xpos': xpos }).update({ '$set': { 'seg_dup': xpos_seg_dups[xpos], 'tandem_repeat': xpos_tandem_repeats[xpos] } })
		counter += 1

		if (counter % 500 == 0):
			bulk.execute()
			bulk = db.hw.variants_hwe_pop.initialize_ordered_bulk_op()

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	if (counter % 500 != 0):
		bulk.execute()

#########################################################
### OBTAINING VARIANT ALLELE BALANCE DATA FROM GNOMAD ###
#########################################################

class AlleleBalance():
	def __init__(self):
		self.variant_id = ''
		self.af = 0.0
		self.flags = []
		self.filters = []
		self.exome_and_genome = False
		self.is_canonical = False
		self.ab_mid_ratio = -1
		self.ab_values = []
		self.e_site_qm = {}
		self.g_site_qm = {}

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['_id'] = self.variant_id
		dictionary['variant_id'] = self.variant_id
		dictionary['af'] = self.af
		dictionary['flags'] = self.flags
		dictionary['filters'] = self.filters

		if len(self.filters) == 0:
			dictionary['pass'] = True
		else:
			dictionary['pass'] = False

		dictionary['is_canonical'] = self.is_canonical
		dictionary['exome_and_genome'] = self.exome_and_genome
		
		dictionary['ab_mid_ratio'] = self.ab_mid_ratio
		dictionary['ab_values'] = self.ab_values
		dictionary['e_site_qm'] = self.e_site_qm
		dictionary['g_site_qm'] = self.g_site_qm
		return dictionary


# Allele Balance Calculation
def parse_gnomad_variant_ab_data(variant, genome=False):
	if genome:
		gnomad_variant = variant['genome']
	else:
		gnomad_variant = variant['exome']

	if not gnomad_variant:
		return {}, {}

	site_qm = gnomad_variant['qualityMetrics']['siteQualityMetrics']

	ab_bin_edges = gnomad_variant['qualityMetrics']['alleleBalance']['alt']['bin_edges']
	ab_bin_freq = gnomad_variant['qualityMetrics']['alleleBalance']['alt']['bin_freq']

	allele_balance = OrderedDict()
	for x in range(0, len(ab_bin_freq)):
		allele_balance[ab_bin_edges[x + 1]] = ab_bin_freq[x]
	return allele_balance, site_qm


def get_variant_ab_balance(variant_id):
	try:
		variant = get_variant_details(variant_id, timeout=60)
	except:
		return -1, {}, {}, []

	if not variant:
		print variant_id

	e_allele_balance, e_site_qm = parse_gnomad_variant_ab_data(variant, genome=False)
	g_allele_balance, g_site_qm = parse_gnomad_variant_ab_data(variant, genome=True)

	allele_balance = OrderedDict()

	if len(e_allele_balance) > 0 and len(g_allele_balance) > 0:
		for value_bin in e_allele_balance.keys():
			allele_balance[value_bin] = e_allele_balance[value_bin] + g_allele_balance[value_bin]
	elif len(e_allele_balance) > 0:
		allele_balance = e_allele_balance
	elif len(g_allele_balance) > 0:
		allele_balance = g_allele_balance
	else:
		return -1, e_site_qm, g_site_qm, []

	if sum(allele_balance.values()) > 0:
		ab_mid_ratio = float(allele_balance[0.45] + allele_balance[0.5] + allele_balance[0.55]) / sum(allele_balance.values())
	else:
		return 0, e_site_qm, g_site_qm, []

	return ab_mid_ratio, e_site_qm, g_site_qm, allele_balance.values()

####################################
### RARE VARIANTS ALLELE BALANCE ###
####################################

class AlleleBalanceBase():
	def __init__(self, variant_id):
		self.variant_id = variant_id
		self.ab_mid_ratio = -1
		self.ab_values = []
		self.e_site_qm = {}
		self.g_site_qm = {}

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['_id'] = self.variant_id
		dictionary['variant_id'] = self.variant_id
		dictionary['ab_mid_ratio'] = self.ab_mid_ratio
		dictionary['ab_values'] = self.ab_values
		dictionary['e_site_qm'] = self.e_site_qm
		dictionary['g_site_qm'] = self.g_site_qm
		return dictionary


def create_rare_variants_ab(db):
	db.hw.rare_variants_ab.drop()
	pop_het_excess_sign_thresholds = {}
	for pop, individuals_num in POP_SIZES.iteritems():
		pop_het_excess_sign_thresholds[pop] = calculate_single_pop_significant_af_thresholds_data(individuals_num, 0)['af']

	unique_variants = set([])

	variants = db.hw.variants_hwe_pop.find({ "all_pop_an_ratio_pass": True, "alt_af": { "$lt": 0.001 }})
	total_lines = variants.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for variant in variants:
		pop = variant['pop']
		af = variant['af']

		if af >= pop_het_excess_sign_thresholds[pop] and af <= RARE_HET_EXCESS_MAX_AF:
			unique_variants.add(variant['variant_id'])

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	for variant_id in unique_variants:
		variant_ab_base = AlleleBalanceBase(variant_id)
		variant_ab_base = variant_ab_base.get_dictionary()
		db.hw.rare_variants_ab.insert(variant_ab_base)


	# Update gnomAD variants Allele Balance (AB) using gnomAD api

	variants = db.hw.rare_variants_ab.find({ "ab_mid_ratio": -1.0 })

	vars_ab = {}
	total_lines = variants.count()
	counter = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for variant in variants:
		var_id = variant['variant_id']

		ab = AlleleBalance()
		ab.ab_mid_ratio, ab.e_site_qm, ab.g_site_qm, ab.ab_values = get_variant_ab_balance(var_id)
		vars_ab[var_id] = ab
		counter += 1
		time.sleep(1)
		if (counter % 50 == 0):
			update_rare_variant_ab_in_db(db, vars_ab)
			vars_ab = {}

		bar.update((counter + 0.0) / total_lines)
	bar.finish()

	if (counter % 50 != 0):
		update_rare_variant_ab_in_db(db, vars_ab)


def update_rare_variant_ab_in_db(db, vars_ab):
	bulk = db.hw.rare_variants_ab.initialize_unordered_bulk_op()
	for var_id, ab_data in vars_ab.iteritems():
		# process in bulk
		bulk.find({ '_id': var_id }).update({ '$set': { 'ab_mid_ratio': ab_data.ab_mid_ratio,
															   'ab_values': ab_data.ab_values,
															   'e_site_qm': ab_data.e_site_qm,
															   'g_site_qm': ab_data.g_site_qm} })
	bulk.execute()

################################
### RARE HET EXCESS VARIANTS ###
################################

def add_clin_var_data_to_variant(db, variant):
	gene_name = variant['gdit_gene_name']
	clinvar_gene = db.hw.clinvar_variants.find_one({'_id': gene_name})
	clinvar_variants = clinvar_gene['variants']
	variant['clinvar'] = 'N'
	variant['clinvar_significance'] = ''
	variant['clinvar_stars'] = ''
	for clinvar_variant in clinvar_variants:
		if clinvar_variant['variantId'] == variant['variant_id']:
			variant['clinvar'] = 'Y'
			variant['clinvar_significance'] = clinvar_variant['clinicalSignificance']
			variant['clinvar_stars'] = clinvar_variant['goldStars']


def add_omim_data_to_variant(db, variant):
	gdit_gene = db.hw.gdit_genes.find_one({'gene': variant['gdit_gene_name']})
	if gdit_gene:
		variant['omim'] = gdit_gene['omim']
		variant['Inheritance_pattern'] = gdit_gene['Inheritance_pattern']
		variant['phenotype'] = gdit_gene['phenotype']
	else:
		variant['omim'] = 'N'
		variant['Inheritance_pattern'] = ''
		variant['phenotype'] = ''


def add_quality_metrics_to_variant(db, variant):
	variant_ab = db.hw.rare_variants_ab.find_one({'variant_id': variant['variant_id']})
	if variant_ab:
		variant['ab_mid_ratio'] = variant_ab['ab_mid_ratio']
		variant['ab_values'] = variant_ab['ab_values']
		variant['e_site_qm'] = variant_ab['e_site_qm']
		variant['g_site_qm'] = variant_ab['g_site_qm']
	else:
		ab_mid_ratio, e_site_qm, g_site_qm, allele_balance_values = get_variant_ab_balance(variant['variant_id'])	
		variant['ab_mid_ratio'] = ab_mid_ratio
		variant['ab_values'] = allele_balance_values
		variant['e_site_qm'] = e_site_qm
		variant['g_site_qm'] = g_site_qm
		time.sleep(1)


def add_seg_dup_to_variant(db, variant):
	xpos = variant['xpos']
	seg_dup = db.hw.seg_dups.find_one({ "xstart": { "$lte": xpos }, "xend": { "$gte": xpos } })
	if seg_dup:
		variant['seg_dup'] = True
	else:
		variant['seg_dup'] = False


def add_tandem_repeat_to_variant(db, variant):
	xpos = variant['xpos']
	tandem_repeat = db.hw.tandem_repeats.find_one({"_id": xpos})

	if tandem_repeat:
		variant['tandem_repeat'] = True
	else:
		variant['tandem_repeat'] = False


def create_rare_het_excess_variants(db):
	db.hw.rare_het_excess_variants.drop()

	# FILTERS: Initial filters for rare heterozygote excess variants
	variants = db.hw.variants_hwe_pop.find({"balance": "-",
											"hwe_p_value": { "$lte": HWE_P_VALUE_THRESHOLD },
											"all_pop_an_ratio_pass": True,
											"all_pop_negative_balance": True,
											"alt_af": { "$lt": 0.001 },
											"max_pop_af": { "$lte": RARE_HET_EXCESS_MAX_AF }, # "af" can be used as alternative 
											"tandem_repeat": False, 
											"seg_dup": False,
											})

	total_lines = variants.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for variant in variants:
		line_number += 1
		
		add_clin_var_data_to_variant(db, variant)
		add_omim_data_to_variant(db, variant)
		add_quality_metrics_to_variant(db, variant)
		# FILTER: Exclude variants with low allele balance
		if variant['ab_mid_ratio'] < ALLELE_BALANCE_MID_RATIO_THRESHOLD:
			continue

		db.hw.rare_het_excess_variants.insert(variant)

		bar.update((line_number + 0.0) / total_lines)
	bar.finish()


#################
### GNOMAD v3 ###
#################

def export_rare_het_exc_variants_for_lift_over(db):
	variants = db.hw.rare_het_excess_variants.find({})
	table = []
	for variant in variants:
		variant_id = variant['variant_id']
		chrom, pos, ref, alt = variant_id.split('-')
		row = ['chr' + chrom + ':' + str(pos) + '-' + str(pos), variant_id]
		table.append(row)

	output_csv = OUTPUT_FOLDER + 'het_exc_variants_bed.csv'
	write_table_to_csv(table, output_csv)


def import_rare_het_exc_variants_lift_over_results(db):
	input_file = open(OUTPUT_FOLDER + 'het_exc_variants_bed.csv', 'rt')
	reader = csv.reader(input_file)
	variant_ids = []
	for row in reader:
		variant_ids.append(row[1])

	input_file = open(OUTPUT_FOLDER + 'het_exc_variants_build_38.csv', 'rt')
	reader = csv.reader(input_file)

	variant_poses = []
	for row in reader:
		chrom, start_stop = row[0].split(':')
		start, stop = start_stop.split('-')
		variant_poses.append(start)

	variant_id_37_to_38 = OrderedDict()
	for x in range(0, len(variant_ids)):
		chrom, pos, ref, alt = variant_ids[x].split('-')
		variant_id_38 = '-'.join([chrom, variant_poses[x], ref, alt])
		variant_id_37_to_38[variant_ids[x]] = variant_id_38

	data = {'_id': 'rare_het_exc_variants_38', 'variant_id_37_to_38': variant_id_37_to_38}
	db.hw.stats.delete_one({'_id': 'rare_het_exc_variants_38'})
	db.hw.stats.insert(data)


def update_rare_het_exc_variants_with_gnomad_3_data(db):
	pops_max_an_g3 = OrderedDict()
	pops_max_an_g3['AFR'] = 21042 * 2
	pops_max_an_g3['AMR'] = 6835 * 2
	pops_max_an_g3['ASJ'] = 1662 * 2
	pops_max_an_g3['EAS'] = 1567 * 2
	pops_max_an_g3['FIN'] = 5244 * 2
	pops_max_an_g3['NFE'] = 32299 * 2
	pops_max_an_g3['SAS'] = 1526 * 2

	variant_ids = db.hw.stats.find_one({"_id": "rare_het_exc_variants_38"})
	variant_ids = variant_ids["variant_id_37_to_38"]
	variant_pops = {}

	for variant_id in variant_ids:
		het_exc_variants = db.hw.rare_het_excess_variants.find({'variant_id': variant_id})
		variant_pops[variant_id] = []
		for het_exc_variant in het_exc_variants:
			variant_pops[variant_id].append(het_exc_variant['pop'])

	total_lines = len(variant_ids)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for variant_id, variant_id_38 in variant_ids.iteritems():
		pops = variant_pops[variant_id]
		variant_38 = get_variant_details_gnomad_3(variant_id_38, timeout=120)
		variant_38 = variant_38['genome']
		filters = variant_38['filters']
		if len(filters) == 0:
			filters = 'PASS'
		else:
			filters = ', '.join(filters)
		pops_data = variant_38['populations']
		for pop_data in pops_data:
			pop = pop_data['id']
			if pop not in pops:
				continue
			ac = pop_data['ac']
			an = pop_data['an']
			hom = pop_data['ac_hom']
			an_ratio = float(an) / pops_max_an_g3[pop]

			het = ac - 2 * hom
			af = float(ac) / an
			ref_hom = an / 2 - het - hom
			exp_hom = af * af * an / 2 # q^2 * indiviudlas number (an / 2)
			exp_het = 2 * (1 - af) * af * an / 2 # 2pq * indiviudlas number (an / 2); (p = 1 - q)
			hwe_p_value = hwe(het, hom, ref_hom)

			if hom < exp_hom:
				balance = '-'
			else:
				balance = '+'

			db.hw.rare_het_excess_variants.update({"variant_id": variant_id, "pop": pop}, 
												  { '$set': { 'ac_g3': ac,
															  'an_g3': an,
															  'an_ratio_g3': an_ratio,
															  'af_g3': af,
															  'het_g3': het,
															  'exp_het_g3': exp_het,
															  'hom_g3': hom,
															  'exp_hom_g3': exp_hom,
															  'hwe_p_value_g3': hwe_p_value,
															  'balance_g3': balance,
															  'filters_g3': filters }})
		time.sleep(1)
		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()


def main():
	db = MongoDB()
	# For the first run, please carefully read the comments and uncomment all of the following functions:

	# These functions create "variants_hwe_pop" collection which stores 
	# data related to variant deviation from Hardy-Weinberg Equilibrium in each population.
	# "variants_hwe_pop" collection is then updated with aggregated frequencies of alternative alleles (alt_af),
	# which are used for further variant filtering (high alt_af might compromise HWE analysis results).
	#analyse_gnomad_variants_deviations_from_hwe(db)
	#add_alt_af_data(db)

	# These functions create temporary "variants_hwe_regions" collection
	# and add flags to "variants_hwe_pop" which indicate whether variants
	# are located in tandem repeat and segmental duplication regions or not.
	#create_hwe_variants_regions(db)
	#update_variants_hwe_pop_with_region_data(db)

	# This function obtains Allele Balance (AB) data from gnomAD for rare variants (i.e. 0.001<=AF<=RARE_HET_EXCESS_MAX_AF (0.1))
	# It creates "rare_variants_ab" collection, which is used for further variant filtering (low AB might be a sign of sequencing errors).
	# !!! IMPORTANT !!! 
	# It works with gnomAD API (i.e. requires internet connection) and can take up to a couple of days to run!!!
	#create_rare_variants_ab(db)

	# This function creates a dataset ("rare_het_excess_variants" collection) of variants with heterozygous excess (HetExc)
	#create_rare_het_excess_variants(db)

	# These functions export HetExc variant chromosomal coordinates for LiftOver conversion (this has to be done manually),
	# import converted variants (./tables/het_exc_variants_build_38.csv) back to the "rare_het_excess_variants" collection
	# and use these new coordinates to get allele data from gnomAD v3 via API.	
	#export_rare_het_exc_variants_for_lift_over(db)
	#import_rare_het_exc_variants_lift_over_results(db)
	#update_rare_het_exc_variants_with_gnomad_3_data(db)

if __name__ == "__main__":
	sys.exit(main())