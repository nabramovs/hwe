import os
import sys
import progressbar
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from scipy import stats
from scipy.stats import fisher_exact
from scipy.stats import shapiro
from collections import OrderedDict
from common import MongoDB, write_table_to_csv, get_confidence_intervals
from hardy_weinberg import Hardy_Weinberg_Equilibrium_exact_test_user_Kantale as hwe

HWE_P_VALUE_THRESHOLD = 0.05
RARE_HET_EXCESS_MAX_AF = 0.05
ALLELE_BALANCE_MID_RATIO_THRESHOLD = 0.5

# NOTE: if you recreate "rare_het_excess_variants" dataset with different filters
# To check clinvar consequences in rare het excess variants set run:
# report_clinvar_statuses_in_rare_het_excess_genes(db)


CLINVAR_BENIGN_STATUSES = set(['Benign', 'Likely_benign', 'Benign/Likely_benign'])
CLINVAR_CONFLICTING_STATUSES = set(['Conflicting_interpretations_of_pathogenicity'])
CLINVAR_PATHOGENIC_STATUSES = set(['Pathogenic'])
CLINVAR_UNKNOWN_STATUSES = set(['', 'not_provided'])

OUTPUT_FOLDER = './tables/'
FIGURES_FOLDER = './figures/'

#SUBPLOT_LETTERS = {1: 'a', 2: 'b', 3: 'c', 4: 'd', 5: 'e', 6: 'f', 7: 'g', 8: 'h', 9: 'i'}
SUBPLOT_LETTERS = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 5: 'E', 6: 'F', 7: 'G', 8: 'H', 9: 'I'}

# NOTE: population names in POP_SIZES and POP_COLOURS should be in the same order:
POP_ORDER = ['NFE', 'AMR', 'SAS', 'FIN', 'AFR', 'EAS', 'ASJ']

POP_SIZES = OrderedDict()
POP_SIZES['NFE'] = 64603
POP_SIZES['AMR'] = 17720
POP_SIZES['SAS'] = 15308
POP_SIZES['FIN'] = 12562
POP_SIZES['AFR'] = 12487
POP_SIZES['EAS'] = 9977
POP_SIZES['ASJ'] = 5185

GNOMAD_INDIVIDUALS_NUM = 141456

POP_COLOURS = OrderedDict()
POP_COLOURS['NFE'] = '#6aa5cd'
POP_COLOURS['AMR'] = '#ed2324'
POP_COLOURS['SAS'] = '#ff9912'
POP_COLOURS['FIN'] = '#1b77b7'
POP_COLOURS['AFR'] = '#942194'
POP_COLOURS['EAS'] = '#108c43'
POP_COLOURS['ASJ'] = '#ff7f50'

POP_FULL_NAMES = OrderedDict()
POP_FULL_NAMES['NFE'] = 'Non-Finnish European (NFE)'
POP_FULL_NAMES['AMR'] = 'Latino/Admixed American (AMR)'
POP_FULL_NAMES['SAS'] = 'South Asian (SAS)'
POP_FULL_NAMES['FIN'] = 'Finnish (FIN)'
POP_FULL_NAMES['AFR'] = 'African/African American (AFR)'
POP_FULL_NAMES['EAS'] = 'East Asian (EAS)'
POP_FULL_NAMES['ASJ'] = 'Ashkenazi Jewish (ASJ)'

C_LIGHT_GRAY = '#d9d9d9'
C_RED = '#ff3300'
C_GREEN = '#00cc00'
C_ORANGE = '#ff7f0e' # '#ff8c1a'
C_BLUE = '#1f77b4' # '#0066ff'
C_BLACK = '#000000'


def export_fisher_test_results(table_name, title, group_names, variable_names, group_1_values, group_2_values, fold_enrichemnt, p_value, subtitle=''):
	table = [[title, '', ''],
			 [subtitle] + variable_names,
			 [group_names[0]] + group_1_values,
			 [group_names[1]] + group_2_values,
			 ['fold_enrichemnt:', fold_enrichemnt, ''],
			 ['p-value:', p_value, ''],
			] 
	output_csv = OUTPUT_FOLDER + table_name
	write_table_to_csv(table, output_csv)


def bar_chart_get_bar_width_and_x_paddings(n_bars, base_width=0.8):
	bar_width = base_width / n_bars

	if n_bars % 2 == 0:
		x_paddings = []
		left_bars_n = n_bars / 2
		for x in range(0, n_bars / 2):
			x_paddings.append(-1 * bar_width * left_bars_n + 0.5 * bar_width)
			left_bars_n -= 1

		right_bars_n = 1
		for x in range(0, n_bars / 2):
			x_paddings.append(bar_width * right_bars_n - 0.5 * bar_width)
			right_bars_n += 1
	else:
		x_paddings = []
		left_bars_n = (n_bars - 1) / 2
		for x in range(0, (n_bars - 1) / 2):
			x_paddings.append(-1 * bar_width * left_bars_n)
			left_bars_n -= 1

		# middle bar
		x_paddings.append(0)
		right_bars_n = 1
		for x in range(0, n_bars / 2):
			x_paddings.append(bar_width * right_bars_n)
			right_bars_n += 1
	return bar_width, x_paddings


# Significance line drawing code was taken from here:
# https://stackoverflow.com/questions/11517986/indicating-the-statistically-significant-difference-in-bar-graph
def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None, linewidth=1, maxasterix=None):
	""" 
	Annotate barplot with p-values.

	:param num1: number of left bar to put bracket over
	:param num2: number of right bar to put bracket over
	:param data: string to write or number for generating asterixes
	:param center: centers of all bars (like plt.bar() input)
	:param height: heights of all bars (like plt.bar() input)
	:param yerr: yerrs of all bars (like plt.bar() input)
	:param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
	:param barh: bar height in axes coordinates (0 to 1)
	:param fs: font size
	:param linewidth: width of the significance line
	:param maxasterix: maximum number of asterixes to write (for very small p-values)
	"""

	if type(data) is str:
		text = data
	else:
		# * is p <= 0.05
		# ** is p < 0.01
		# *** is p < 0.001
		# etc.
		text = ''
		p = .01

		if data <= 0.05:
			text += '*'

		while data < p:
			text += '*'
			p /= 10.

			if maxasterix and len(text) == maxasterix:
				break

		if len(text) == 0:
			text = 'n.s.'

	lx, ly = center[num1], height[num1]
	rx, ry = center[num2], height[num2]

	if yerr:
		ly += yerr[num1]
		ry += yerr[num2]

	ax_y0, ax_y1 = plt.gca().get_ylim()
	dh *= (ax_y1 - ax_y0)
	barh *= (ax_y1 - ax_y0)

	y = max(ly, ry) + dh

	barx = [lx, lx, rx, rx]
	bary = [y, y+barh, y+barh, y]
	mid = ((lx+rx)/2.0, y+barh)

	plt.plot(barx, bary, c='black', linewidth=linewidth)

	kwargs = dict(ha='center', va='bottom')
	if fs is not None:
		kwargs['fontsize'] = fs

	plt.text(mid[0], mid[1], text, **kwargs)


##############################################
### POP HET EXCESS SIGNIFICANCE THRESHOLDS ###
##############################################

def calculate_single_pop_significant_af_thresholds_data(individuals_num, rare_hom, p_value_threshold=HWE_P_VALUE_THRESHOLD):
	het = 0
	af = 0.0
	common_hom = individuals_num - rare_hom - het

	reached_threshold = False

	while not reached_threshold:
		het += 1
		common_hom -= 1
		p_value = hwe(het, rare_hom, common_hom, mid_p=True)
		af = het / float(individuals_num * 2)

		if p_value < p_value_threshold and not reached_threshold:
			reached_threshold = True
			pop_data = OrderedDict()
			pop_data['individuals'] = individuals_num
			pop_data['het'] = het
			pop_data['af'] = af
			pop_data['p_value'] = p_value

	return pop_data


def calculate_minimum_population_number_required_to_achive_significance_for_af(af, p_value_threshold=HWE_P_VALUE_THRESHOLD):
	individuals_num = 100
	pop_increase_step = 100
	rare_hom = 0
	reached_threshold = False
	while not reached_threshold:
		het = int(individuals_num * af)
		common_hom = individuals_num - het
		p_value = hwe(het, rare_hom, common_hom, mid_p=True)
		if p_value < p_value_threshold and not reached_threshold:
			reached_threshold = True
			pop_data = OrderedDict()
			pop_data['individuals'] = individuals_num
			pop_data['het'] = het
			pop_data['af'] = af
			pop_data['p_value'] = p_value

		individuals_num += pop_increase_step

	return pop_data			


################
### FIGURE 1 ###
################

F1_SUBPLOT_LABEL_FONTSIZE = 12
F1_TITLE_FONTSIZE = 10
F1_NORMAL_FONTSIZE = 8
F1_AXIS_FONTSIZE = 8
F1_LEGEND_FONTSIZE = 10

def caculate_f1_data(db):
	print '--- Calculating data for Figure 1... ---'
	pop_variants = OrderedDict()
	pop_het_excess_variants = OrderedDict()
	pop_het_deficiency_variants = OrderedDict()
	pop_inbreeding_coef = OrderedDict()
	pop_het_excess_sign_thresholds = OrderedDict()

	for pop, individuals_num in POP_SIZES.iteritems():
		pop_variants[pop] = 0
		pop_het_excess_variants[pop] = 0
		pop_het_deficiency_variants[pop] = 0
		pop_inbreeding_coef[pop] = []
		pop_het_excess_sign_thresholds[pop] = calculate_single_pop_significant_af_thresholds_data(individuals_num, 0)

	# FILTER: exclude variants which site is not covered in at least 80% of each population and have common alternative alleles
	variants = db.hw.variants_hwe_pop.find({ "all_pop_an_ratio_pass": True, "alt_af": { "$lt": 0.001 }})

	total_lines = variants.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for variant in variants:
		pop = variant['pop']
		pop_variants[pop] += 1
		# FILTER: for inbreeding coefficient include only relatively common variants
		pop_inbreeding_coef[pop].append(variant['inbreeding_coef'])

		if variant['hwe_p_value'] < HWE_P_VALUE_THRESHOLD:
			if variant['balance'] == '+':
				pop_het_deficiency_variants[pop] += 1
			elif variant['balance'] == '-':
				pop_het_excess_variants[pop] += 1
		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	f1_data = OrderedDict()
	f1_data['_id'] = 'f1_data'
	f1_data['pop_sizes'] = POP_SIZES
	f1_data['pop_variants'] = pop_variants
	f1_data['pop_het_excess_variants'] = pop_het_excess_variants
	f1_data['pop_het_deficiency_variants'] = pop_het_deficiency_variants
	f1_data['pop_inbreeding_coef'] = pop_inbreeding_coef
	f1_data['pop_het_excess_sign_thresholds'] = pop_het_excess_sign_thresholds

	db.hw.stats.delete_one({'_id': f1_data['_id']})
	db.hw.stats.insert(f1_data)


def draw_f1_pop_bar_subplot(subplot_num, xs, y_label, yerr=[], title='', annotation_type='', xticks=[], x_offset=0, subplot_letter_x_offset=-0.25):
	ax = plt.subplot(3,2,subplot_num)
	ax.text(subplot_letter_x_offset, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=F1_SUBPLOT_LABEL_FONTSIZE, fontweight='bold', va='top', ha='right')
	
	ys = POP_ORDER
	colours = []
	for pop in POP_ORDER:
		colours.append(POP_COLOURS[pop])

	width = 0.8 # the width of the bars 
	ind = np.arange(len(ys))  # the y locations for the groups

	if yerr:
		ax.bar(ind, xs, width, yerr=yerr, capsize=7, color=colours)
	else:
		ax.bar(ind, xs, width, color=colours)
	
	for i, v in enumerate(xs):
		if annotation_type == 'int':
			if v >= 100000:
				y_text_offset = max(xs) * 0.11
			elif v >= 10000:
				y_text_offset = max(xs) * 0.09
			else:
				y_text_offset = max(xs) * 0.07
			ax.text(i-0.6, v+y_text_offset, "  {:,}".format(v), va='center', fontsize=F1_NORMAL_FONTSIZE, rotation=45)
		elif annotation_type == 'float':
			y_text_offset = max(xs) * 0.09
			ax.text(i-0.6, v+y_text_offset, "  {0:.4f}".format(v), va='center', fontsize=F1_NORMAL_FONTSIZE, rotation=45)
		elif annotation_type == 'percent':
			y_text_offset = max(xs) * 0.03
			ax.text(i, v+y_text_offset, " {0:.1f}%".format(v), va='center', ha='center', fontsize=F1_NORMAL_FONTSIZE) # , rotation=45
		else:
			ax.text(i-0.3, v+y_text_offset, "  {}".format(v), va='center', ha='center', fontsize=F1_NORMAL_FONTSIZE, rotation=45)

	if title:
		ax.text(4.5, max(xs), title, va='center', fontsize=F1_TITLE_FONTSIZE, fontweight='bold')

	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)

	ax.set_xticks(ind)
	ax.tick_params(axis='both', which='major', labelsize=F1_AXIS_FONTSIZE)
	ax.set_xticklabels(ys, minor=False, fontsize=F1_AXIS_FONTSIZE)
	ax.set_ylabel(y_label, fontsize=F1_AXIS_FONTSIZE)
	ax.set_ylim(ymax=max(xs) + x_offset)

	ax.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=True) # labels along the bottom edge are off

	if xticks:
		ax.set_xticks(xticks)

	if annotation_type == 'int':
		ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))


def draw_f1(db):
	print '--- Drawing Figure 1... ---'
	fig = plt.figure(1)

	f1_data = db.hw.stats.find_one({'_id': 'f1_data'})

	pop_sizes = []
	pop_variants = []
	pop_af_limits = []
	pop_het_deficiency_percents = []
	pop_het_excess_percents = []

	for pop in POP_ORDER:
		pop_sizes.append(f1_data['pop_sizes'][pop])
		pop_variants.append(f1_data['pop_variants'][pop])
		pop_af_limits.append(f1_data['pop_het_excess_sign_thresholds'][pop]['af'])
		pop_het_deficiency_percents.append(f1_data['pop_het_deficiency_variants'][pop] * 100 / float(f1_data['pop_variants'][pop]))
		pop_het_excess_percents.append(f1_data['pop_het_excess_variants'][pop] * 100 / float(f1_data['pop_variants'][pop]))

	# Figure A Population Sizes
	subplot_num = 1
	draw_f1_pop_bar_subplot(subplot_num, pop_sizes, 'Number of individuals', annotation_type='int', x_offset=10000, subplot_letter_x_offset=-0.24)

	# Figure B Population Variants
	subplot_num += 1	
	draw_f1_pop_bar_subplot(subplot_num, pop_variants, 'Number of variants', annotation_type='int', x_offset=30000, subplot_letter_x_offset=-0.28)

	# Figure C Population Het Excess Significance Thresholds
	subplot_num += 1
	draw_f1_pop_bar_subplot(subplot_num, pop_af_limits, 'Allele Frequency (AF)', annotation_type='float', x_offset=0.005, subplot_letter_x_offset=-0.24)

	# Figure D Population Heterozygote Deficiency
	subplot_num += 1
	draw_f1_pop_bar_subplot(subplot_num, pop_het_deficiency_percents, 'Variants (%)', annotation_type='percent', x_offset=3, subplot_letter_x_offset=-0.27, title='HetDef')	

	# Figure E Population Heterozygote Excess
	subplot_num += 1
	draw_f1_pop_bar_subplot(subplot_num, pop_het_excess_percents, 'Variants (%)', annotation_type='percent', x_offset=0.1, subplot_letter_x_offset=-0.24, title='HetExc')		

	# Figure Legend (Full population names)
	pop_legend = OrderedDict()
	for pop in POP_ORDER:
		pop_legend[POP_FULL_NAMES[pop]] = mpatches.Patch(color=POP_COLOURS[pop])
	
	fig.legend(pop_legend.values(), pop_legend.keys(), loc=[0.53,0.07], frameon=False, ncol=1, fontsize=F1_LEGEND_FONTSIZE)

	fig = plt.figure(1)
	fig.set_size_inches(7, 8)
	plt.tight_layout(rect=[0.02, 0.02, 0.99, 0.99], h_pad=2)
	plt.savefig(FIGURES_FOLDER + 'F1.png', format='png', dpi=300)
	plt.close(fig)


################
### FIGURE 2 ###
################

def get_list_of_values_from_dict_by_keys(keys, dictionary):
	values = []
	for key in keys:
		values.append(dictionary[key])
	return values


def calculate_f2_data(db):
	print '--- Calculating data for Figure 2... ---'
	variants = db.hw.rare_variants_ab.find({ "min_alt_af": { "$lt": 0.001 }, "all_pop_an_ratio_pass": True, "max_pop_af": { "$lte": RARE_HET_EXCESS_MAX_AF }})

	variants_ab = {}
	all_variants = set([])
	seg_dup_variants = set([])
	tandem_repeat_variants = set([])
	het_excess_variants = set([])

	ab_proportions = OrderedDict()
	for x in range(0, 20):
		ab_proportions[x] = []

	total_lines = variants.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for variant in variants:
		variant_id = variant['variant_id']
		all_variants.add(variant_id)
		if variant['seg_dup']:
			seg_dup_variants.add(variant_id)
		if variant['tandem_repeat']:
			tandem_repeat_variants.add(variant_id)
		if variant['min_het_exc_p_value'] < 0.05 and variant['all_pop_negative_balance']:
			het_excess_variants.add(variant_id)

		variants_ab[variant_id] = variant['ab_mid_ratio'] * 100

		ab_values = variant['ab_values']
		ab_samples = float(sum(ab_values))
		# Estimate mid AB ratios using only variants in non seg dup and tandem repeat regions
		if not variant['seg_dup'] and not variant['tandem_repeat']:
			for x in range(0, len(ab_values)):
				ab_proportions[x].append(ab_values[x] * 100 / ab_samples)

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	not_in_repeat_variants = all_variants - seg_dup_variants - tandem_repeat_variants
	not_in_repeat_het_excess_variants = not_in_repeat_variants & het_excess_variants
	seg_dup_het_excess_variants = seg_dup_variants & het_excess_variants
	tandem_repeat_het_excess_variants = tandem_repeat_variants & het_excess_variants

	not_in_repeat_variants_abs = get_list_of_values_from_dict_by_keys(not_in_repeat_variants, variants_ab)
	not_in_repeat_het_excess_variants_abs = get_list_of_values_from_dict_by_keys(not_in_repeat_het_excess_variants, variants_ab)
	seg_dup_variants_abs = get_list_of_values_from_dict_by_keys(seg_dup_variants, variants_ab)
	tandem_repeat_variants_abs = get_list_of_values_from_dict_by_keys(tandem_repeat_variants, variants_ab)

	not_in_repeat_variants_ab_low = 0
	for ab in not_in_repeat_variants_abs:
		if ab < ALLELE_BALANCE_MID_RATIO_THRESHOLD * 100:
			not_in_repeat_variants_ab_low += 1
	not_in_repeat_variants_ab_low_percent = float(not_in_repeat_variants_ab_low) / len(not_in_repeat_variants_abs)

	not_in_repeat_het_excess_variants_ab_low = 0
	for ab in not_in_repeat_het_excess_variants_abs:
		if ab < ALLELE_BALANCE_MID_RATIO_THRESHOLD * 100:
			not_in_repeat_het_excess_variants_ab_low += 1
	not_in_repeat_het_excess_variants_ab_low_percent = float(not_in_repeat_het_excess_variants_ab_low) / len(not_in_repeat_het_excess_variants_abs)

	mean_abs = OrderedDict()
	err_abs = OrderedDict()

	low_lim = 0
	step = 0.05
	for x in range(0, 20):
		bin_name = str(low_lim) + '-' + str(low_lim+step)
		bin_name = bin_name.replace('.','_')
		mean_abs[bin_name] = np.mean(ab_proportions[x])
		err_abs[bin_name] = np.std(ab_proportions[x])
		low_lim += step

	f2_data = OrderedDict()
	f2_data['_id'] = 'f2_data'
	f2_data['not_in_repeat_variants'] = len(not_in_repeat_variants)
	f2_data['not_in_repeat_het_excess_variants'] = len(not_in_repeat_het_excess_variants)
	f2_data['not_in_repeat_het_excess_variants_percent'] = len(not_in_repeat_het_excess_variants) / float(len(not_in_repeat_variants))
	
	f2_data['seg_dup_variants'] = len(seg_dup_variants)
	f2_data['seg_dup_het_excess_variants'] = len(seg_dup_het_excess_variants)
	f2_data['seg_dup_het_excess_variants_percent'] = len(seg_dup_het_excess_variants) / float(len(seg_dup_variants))

	f2_data['tandem_repeat_variants'] = len(tandem_repeat_variants)
	f2_data['tandem_repeat_het_excess_variants'] = len(tandem_repeat_het_excess_variants)
	f2_data['tandem_repeat_het_excess_variants_percent'] = len(tandem_repeat_het_excess_variants) / float(len(tandem_repeat_variants))

	f2_data['all_mean_bin_abs'] = mean_abs
	f2_data['all_err_bin_abs'] = err_abs

	f2_data['not_in_repeat_variants_abs'] = not_in_repeat_variants_abs
	f2_data['not_in_repeat_het_excess_variants_abs'] = not_in_repeat_het_excess_variants_abs
	f2_data['seg_dup_variants_abs'] = seg_dup_variants_abs
	f2_data['tandem_repeat_variants_abs'] = tandem_repeat_variants_abs

	f2_data['not_in_repeat_variants_ab_low'] = not_in_repeat_variants_ab_low
	f2_data['not_in_repeat_variants_ab_low_percent'] = not_in_repeat_variants_ab_low_percent
	f2_data['not_in_repeat_het_excess_variants_ab_low'] = not_in_repeat_het_excess_variants_ab_low	
	f2_data['not_in_repeat_het_excess_variants_ab_low_percent'] = not_in_repeat_het_excess_variants_ab_low_percent

	db.hw.stats.delete_one({'_id': f2_data['_id']})
	db.hw.stats.insert(f2_data)

	#get_confidence_intervals(stats, stats_name='', alpha=0.95)

F2_SUBPLOT_LABEL_FONTSIZE = 12
F2_NORMAL_FONTSIZE = 8
F2_AXIS_FONTSIZE = 8


def draw_f2(db):
	print '--- Drawing Figure 2... ---'
	fig = plt.figure(1)
	f2_data = db.hw.stats.find_one({'_id': 'f2_data'})
	subplot_num = 1
	draw_f2_repeats(subplot_num, f2_data)
	subplot_num += 1
	draw_f2_ab_distribution(subplot_num, f2_data)	
	subplot_num += 1
	draw_f2_repeat_variants_ab(subplot_num, f2_data)
	subplot_num += 1
	draw_f2_low_ab(subplot_num, f2_data)

	fig = plt.figure(1)
	fig.set_size_inches(7, 6)
	plt.tight_layout(rect=[0.02, 0.02, 0.99, 0.99])
	plt.savefig(FIGURES_FOLDER + 'F2.png', format='png', dpi=300)
	plt.close(fig)


def draw_f2_repeat_variants_ab(subplot_num, f2_data):
	ax = plt.subplot(2,2,subplot_num)
	ax.text(-0.17, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=F2_SUBPLOT_LABEL_FONTSIZE, fontweight='bold', va='top', ha='right')

	ref_data = f2_data['not_in_repeat_variants_abs']
	rep_data = f2_data['tandem_repeat_variants_abs'] + f2_data['seg_dup_variants_abs']
	tr_data = f2_data['tandem_repeat_variants_abs']
	sd_data = f2_data['seg_dup_variants_abs']

	get_confidence_intervals(ref_data, stats_name='Ref', alpha=0.80)
	get_confidence_intervals(rep_data, stats_name='Segdup+tandem repeat', alpha=0.80)
	get_confidence_intervals(tr_data, stats_name='Tandem repeat', alpha=0.80)
	get_confidence_intervals(sd_data, stats_name='Segdup', alpha=0.80)
	get_confidence_intervals(ref_data, stats_name='Ref', alpha=0.95)

	ref_weights = np.ones_like(ref_data)/(float(len(ref_data)) / 100)

	sd_weights = np.ones_like(sd_data)/(float(len(sd_data)) / 100)
	tr_weights = np.ones_like(tr_data)/(float(len(tr_data)) / 100)

	ref_name = 'Ref'
	sd_name = 'Segmental\nduplication'
	tr_name = 'Tandem repeat'

	seg_dup_weights = np.ones_like(rep_data)/(float(len(rep_data)) / 100)

	bins = np.linspace(0, 100, 20)
	plt.hist([ref_data, sd_data, tr_data], bins=bins, label=[ref_name, sd_name, tr_name], weights=[ref_weights, sd_weights, tr_weights])

	plt.legend(loc='upper left', frameon=False, ncol=1, fontsize=F2_NORMAL_FONTSIZE)
	ax.set_xlabel('Carriers with normal allele balance (%)', fontsize=F2_AXIS_FONTSIZE)
	ax.set_ylabel('Variants (%)', fontsize=F2_AXIS_FONTSIZE)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.tick_params(axis='both', which='major', labelsize=F2_AXIS_FONTSIZE)
	get_confidence_intervals(f2_data['not_in_repeat_variants_abs'], stats_name='not_in_repeat_variants', alpha=0.95)


def draw_f2_ab_distribution(subplot_num, f2_data):
	ax = plt.subplot(2,2,subplot_num)
	ax.text(-0.15, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=F2_SUBPLOT_LABEL_FONTSIZE, fontweight='bold', va='top', ha='right')

	ab_means = f2_data['all_mean_bin_abs']
	ab_errs = f2_data['all_err_bin_abs'].values()
	mid_abs_ys = [ab_means['0_4-0_45'], ab_means['0_45-0_5'], ab_means['0_5-0_55']]
	bin_names = ab_means.keys()
	mid_ab_xs = [bin_names.index('0_4-0_45'), bin_names.index('0_45-0_5'), bin_names.index('0_5-0_55')]

	ys = ab_means.values()
	xs = [i for i, _ in enumerate(ys)]

	x_labels = []
	for bin_name in ab_means.keys():
		bin_name = bin_name.replace('_', '.')
		x_labels.append(bin_name)

	for i in range(0, len(mid_ab_xs)):
		x = mid_ab_xs[i]
		y = 0.1
		v = mid_abs_ys[i]
		ax.text(x, 12, "{0:.1f}%".format(v), ha='center', rotation=90, fontsize=F2_NORMAL_FONTSIZE)

	plt.bar(xs, ys, yerr=ab_errs, capsize=3, width=0.9)
	plt.bar(mid_ab_xs, mid_abs_ys, width=0.9)

	plt.xticks(xs, x_labels, rotation=90)
	ax.set_xlabel('Allele balance', fontsize=F2_AXIS_FONTSIZE)
	ax.set_ylabel('Variant carriers (%)', fontsize=F2_AXIS_FONTSIZE)

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.tick_params(axis='both', which='major', labelsize=F2_AXIS_FONTSIZE)

	legend = OrderedDict()
	legend['Normal'] = mpatches.Patch(color=C_ORANGE)
	legend['Other'] = mpatches.Patch(color=C_BLUE)
	plt.legend(legend.values(), legend.keys(), loc='upper left', frameon=False, ncol=1, fontsize=F2_NORMAL_FONTSIZE)


def draw_f2_low_ab(subplot_num, f2_data):
	ax = plt.subplot(2,2,subplot_num)
	ax.text(-0.15, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=F2_SUBPLOT_LABEL_FONTSIZE, fontweight='bold', va='top', ha='right')

	nir_ab_low = f2_data['not_in_repeat_variants_ab_low']
	nir_ab_low_percent = f2_data['not_in_repeat_variants_ab_low_percent']
	he_nir_ab_low = f2_data['not_in_repeat_het_excess_variants_ab_low']
	he_nir_ab_low_percent = f2_data['not_in_repeat_het_excess_variants_ab_low_percent']

	xs = [0, 1]
	ys = [nir_ab_low_percent * 100, he_nir_ab_low_percent * 100]
	vals = [nir_ab_low, he_nir_ab_low]
	x_labels = ['Ref All', 'Ref HetExc']
	plt.bar(xs, ys)
	for i in range(0, len(xs)):
		ax.text(xs[i], ys[i] / 2.0, "{:,}".format(vals[i]), ha='center', va='center', color='white', fontsize=F2_NORMAL_FONTSIZE)
		ax.text(xs[i], ys[i] + 1, "{0:.1f}%".format(ys[i]), ha='center', va='center', fontsize=F2_NORMAL_FONTSIZE)

	nir_low_ab_stats = [nir_ab_low, f2_data['not_in_repeat_variants']]
	he_nir_low_ab_stats = [he_nir_ab_low, f2_data['not_in_repeat_het_excess_variants']]
	low_ab_fe, low_ab_p_value = fisher_exact([he_nir_low_ab_stats, nir_low_ab_stats])

	title = 'Comparison of variants with VCNAB < 50%% in the whole Ref group and a subset of variants with statistically significant excess of heterozygotes (HetExc) in Ref group.'
	export_fisher_test_results('f2_d.csv', title, ['Ref', 'Ref HetExc'], ['VCNAB < 50%', 'All'], nir_low_ab_stats, he_nir_low_ab_stats, low_ab_fe, low_ab_p_value)

	print 'LOW AB Ref/HetExc', low_ab_fe, low_ab_p_value

	# To report raw number in scientific notation use this: '%.2E' % low_ab_p_value
	barplot_annotate_brackets(0, 1, low_ab_p_value, xs, ys, dh=.15, barh=.03, fs=F2_NORMAL_FONTSIZE, maxasterix=4)

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)

	ax.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=True) # labels along the bottom edge are off

	ax.set_ylabel('Variants (%)', fontsize=F2_AXIS_FONTSIZE)
	ax.set_ylim(ymax=max(ys) + 5)
	plt.xticks(xs, x_labels)
	ax.tick_params(axis='both', which='major', labelsize=F2_AXIS_FONTSIZE)


def draw_f2_repeats(subplot_num, f2_data):
	ax = plt.subplot(2,2,subplot_num)
	ax.text(-0.17, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=F2_SUBPLOT_LABEL_FONTSIZE, fontweight='bold', va='top', ha='right')

	ys = [f2_data['not_in_repeat_het_excess_variants_percent'] * 100,
		  f2_data['tandem_repeat_het_excess_variants_percent'] * 100,
		  f2_data['seg_dup_het_excess_variants_percent'] * 100,
		 ]
	vals = [f2_data['not_in_repeat_het_excess_variants'],
			f2_data['tandem_repeat_het_excess_variants'],
			f2_data['seg_dup_het_excess_variants'],
		   ]
	xs = range(0, 3)
	x_labels = ['Ref', 'Tandem\nrepeat', 'Segmental\nduplication',]
	plt.bar(xs, ys)
	plt.xticks(xs, x_labels, fontsize=F2_AXIS_FONTSIZE)

	for i in range(0, len(xs)):
		ax.text(xs[i], ys[i] / 2.0, "{:,}".format(vals[i]), ha='center', va='center', color='white', fontsize=F2_NORMAL_FONTSIZE)
		ax.text(xs[i], ys[i] + 0.2, "{0:.1f}%".format(ys[i]), ha='center', va='center', fontsize=F2_NORMAL_FONTSIZE)

	ref_stats = [f2_data['not_in_repeat_het_excess_variants'], f2_data['not_in_repeat_variants']]
	seg_dup_stats = [f2_data['seg_dup_het_excess_variants'], f2_data['seg_dup_variants']]
	tandem_repeat_stats = [f2_data['tandem_repeat_het_excess_variants'], f2_data['tandem_repeat_variants']]
	seg_dup_fe, seg_dup_p_value = fisher_exact([seg_dup_stats, ref_stats])
	title = 'Comparison of variants deviating from HWE due to HetExc which are located in segmental duplication regions and reference (Ref) group (i.e. all other regions except segmental duplications and tandem repeats)'
	export_fisher_test_results('f2_a_1.csv', title, ['Segmental Duplication', 'Ref'], ['HetExc', 'All'], seg_dup_stats, ref_stats, seg_dup_fe, seg_dup_p_value)
	tandem_repeat_fe, tandem_repeat_p_value = fisher_exact([tandem_repeat_stats, ref_stats])
	title = 'Comparison of variants deviating from HWE due to HetExc which are located in tandem repeat regions and reference (Ref) group (i.e. all other regions except segmental duplications and tandem repeats)'
	export_fisher_test_results('f2_a_2.csv', title, ['Tandem Repeat', 'Ref'], ['HetExc', 'All'], tandem_repeat_stats, ref_stats, tandem_repeat_fe, tandem_repeat_p_value)

	print 'Tandem Repeat FE:', tandem_repeat_fe, 'p-value:', tandem_repeat_p_value
	print 'Segmental Duplication FE:', seg_dup_fe, 'p-value:', seg_dup_p_value

	# To report raw number in scientific notation use this: '%.2E' % tandem_repeat_p_value
	barplot_annotate_brackets(0, 1, tandem_repeat_p_value, xs, ys, dh=.15, barh=.03, fs=F2_NORMAL_FONTSIZE, maxasterix=4)
	barplot_annotate_brackets(0, 2, seg_dup_p_value, xs, ys, dh=.15, barh=.03, fs=F2_NORMAL_FONTSIZE, maxasterix=4)

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)

	ax.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=True)  # labels along the bottom edge are off

	ax.set_ylabel('Variants (%)', fontsize=F2_AXIS_FONTSIZE)
	ax.set_ylim(ymax=max(ys) + 1)
	ax.tick_params(axis='both', which='major', labelsize=F2_AXIS_FONTSIZE)


################
### FIGURE 3 ###
################


def report_clinvar_statuses_in_rare_het_excess_genes(db):
	clinvar_statuses = {}
	variants = db.hw.rare_het_excess_variants.find({ "ab_adjusted_hwe_het_exc": True })
	for variant in variants:
		clinvar_status = variant['clinvar_significance']
		if not clinvar_status:
			clinvar_status = 'Empty'
		if clinvar_status not in clinvar_statuses:
			clinvar_statuses[clinvar_status] = 1
		else:
			clinvar_statuses[clinvar_status] += 1

	for clinvar_status, variants_num in clinvar_statuses.iteritems():
		print clinvar_status, variants_num


def caculate_f3_data(db):
	print '--- Calculating data for Figure 3... ---'

	#################
	### HetExcess ###
	#################

	pop_het_excess_variants = OrderedDict()
	for pop in POP_ORDER:
		pop_het_excess_variants[pop] = 0
	pop_het_excess_genes = set([])
	het_excess_unique_variants = set([])
	het_excess_csqs = {}

	pop_het_excess_clinvar = OrderedDict()
	pop_het_excess_clinvar['Pathogenic'] = set([])
	pop_het_excess_clinvar['Conflicting Interpretations of Pathogenicity'] = set([])
	pop_het_excess_clinvar['Benign or Likely_benign'] = set([])
	pop_het_excess_clinvar['Unknown'] = set([])
	
	rare_het_excess_variants = db.hw.rare_het_excess_variants.find({ "ab_adjusted_hwe_het_exc": True })

	for rare_het_excess_variant in rare_het_excess_variants:

		# Delete lists/dictionary fields with extra stats
		del rare_het_excess_variant['ab_values']
		del rare_het_excess_variant['e_site_qm']
		del rare_het_excess_variant['g_site_qm']

		# Update statistics
		pop = rare_het_excess_variant['pop']
		pop_het_excess_variants[pop] += 1

		gene_name = rare_het_excess_variant['gdit_gene_name']
		pop_het_excess_genes.add(gene_name)

		variant_id = rare_het_excess_variant['variant_id']
		if variant_id not in het_excess_unique_variants:
			csq = rare_het_excess_variant['csq']
			if csq in het_excess_csqs:
				het_excess_csqs[csq] += 1
			else:
				het_excess_csqs[csq] = 1


		het_excess_unique_variants.add(variant_id)

		clinvar_significance = rare_het_excess_variant['clinvar_significance']
		if clinvar_significance in CLINVAR_UNKNOWN_STATUSES:
			pop_het_excess_clinvar['Unknown'].add(variant_id)
		elif clinvar_significance in CLINVAR_BENIGN_STATUSES:
			pop_het_excess_clinvar['Benign or Likely_benign'].add(variant_id)
		elif clinvar_significance in CLINVAR_CONFLICTING_STATUSES:
			pop_het_excess_clinvar['Conflicting Interpretations of Pathogenicity'].add(variant_id)
		elif clinvar_significance in CLINVAR_PATHOGENIC_STATUSES:
			pop_het_excess_clinvar['Pathogenic'].add(variant_id)
		else:
			print '!!!WARNING!!!'
			print 'Variant ' + variant_id + ' ClinVar status (' + clinvar_significance + ') is not present in analysed groups'
			print 'You have to redefine ClinVar groups!!!'

	for clivar_group, variant_list in pop_het_excess_clinvar.iteritems():
		pop_het_excess_clinvar[clivar_group] = len(variant_list)

	#####################
	### NOT HetExcess ###
	#####################

	pop_het_excess_sign_thresholds = {}
	for pop, individuals_num in POP_SIZES.iteritems():
		pop_het_excess_sign_thresholds[pop] = calculate_single_pop_significant_af_thresholds_data(individuals_num, 0)['af']

	pop_not_het_excess_variants = OrderedDict()
	for pop in POP_ORDER:
		pop_not_het_excess_variants[pop] = 0
	pop_not_het_excess_genes = set([])

	variant_ab_mid_ratios = {}

	variant_abs = db.hw.rare_variants_ab.find({})
	for variant_ab in variant_abs:
		variant_ab_mid_ratios[variant_ab['variant_id']] = variant_ab['ab_mid_ratio']

	not_het_excess_unique_variants = set([])
	not_het_excess_csqs = {}
	# FILTER: calculate number and genes of not het excess variants
	variants = db.hw.variants_hwe_pop.find({"all_pop_an_ratio_pass": True,
											"alt_af": { "$lt": 0.001 },
											"tandem_repeat": False,
											"seg_dup": False,
											"max_pop_af": { "$lte": RARE_HET_EXCESS_MAX_AF } })
	total_lines = variants.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for variant in variants:
		variant_id = variant['variant_id']
		pop = variant['pop']
		af = variant['af']
		gene_name = variant['gdit_gene_name']
		# FILTER: DO NOT INCLUDE RARE VARIANTS which cannot be potentially heterozygote advantageous
		if af > pop_het_excess_sign_thresholds[pop]: #  and af <= RARE_HET_EXCESS_MAX_AF
			if variant_ab_mid_ratios[variant_id] > ALLELE_BALANCE_MID_RATIO_THRESHOLD:
				if variant_id not in het_excess_unique_variants:
					pop_not_het_excess_variants[pop] += 1
					if variant_id not in not_het_excess_unique_variants:
						csq = variant['csq']
						if csq in not_het_excess_csqs:
							not_het_excess_csqs[csq] += 1
						else:
							not_het_excess_csqs[csq] = 1

					not_het_excess_unique_variants.add(variant_id)
				if gene_name not in pop_het_excess_genes:
					pop_not_het_excess_genes.add(gene_name)

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	not_het_excess_gene_stats = calculate_gene_group_stats(db, pop_not_het_excess_genes)
	het_excess_gene_stats = calculate_gene_group_stats(db, pop_het_excess_genes)

	f3_data = OrderedDict()
	f3_data['_id'] = 'f3_data'
	f3_data['pop_not_het_excess_variants'] = pop_not_het_excess_variants
	f3_data['pop_het_excess_variants'] = pop_het_excess_variants
	f3_data['not_het_excess_gene_stats'] = not_het_excess_gene_stats
	f3_data['het_excess_gene_stats'] = het_excess_gene_stats
	f3_data['pop_het_excess_clinvar'] = pop_het_excess_clinvar
	f3_data['not_het_excess_unique_variants'] = len(not_het_excess_unique_variants)
	f3_data['het_excess_unique_variants'] = len(het_excess_unique_variants)
	f3_data['not_het_excess_csqs'] = not_het_excess_csqs
	f3_data['het_excess_csqs'] = het_excess_csqs

	db.hw.stats.delete_one({'_id': f3_data['_id']})
	db.hw.stats.insert(f3_data)


def get_disease_gene_sets(db):
	ad_genes = set([])
	ar_genes = set([])
	ar_ad_genes = set([])

	genes = db.hw.gdit_genes.find({})
	for gene in genes:
		gene_name = gene['gene']
		omim_inheritance = gene['Inheritance_pattern']
		if omim_inheritance == 'AD':
			ad_genes.add(gene_name)
		if omim_inheritance == 'AR':
			ar_genes.add(gene_name)
		if omim_inheritance == 'AR,AD':
			ar_ad_genes.add(gene_name)
	return ad_genes, ar_genes, ar_ad_genes


def calculate_gene_group_stats(db, gene_set):
	ad_genes, ar_genes, ar_ad_genes = get_disease_gene_sets(db)

	ad_num = len(gene_set & ad_genes)
	ar_num = len(gene_set & ar_genes)
	ar_ad_num = len(gene_set & ar_ad_genes)
	total_num = len(gene_set)
	ad_percent = 100 * ad_num / float(total_num)
	ar_percent = 100 * ar_num / float(total_num)
	ar_ad_percent = 100 * ar_ad_num / float(total_num)

	gene_stats = OrderedDict()
	gene_stats['ad_num'] = ad_num
	gene_stats['ar_num'] = ar_num
	gene_stats['ar_ad_num'] = ar_ad_num
	gene_stats['total_num'] = total_num
	gene_stats['ad_percent'] = ad_percent
	gene_stats['ar_percent'] = ar_percent
	gene_stats['ar_ad_percent'] = ar_ad_percent

	return gene_stats


F3_SUBPLOT_LABEL_FONTSIZE = 12
F3_NORMAL_FONTSIZE = 8
F3_AXIS_FONTSIZE = 8
F3_LEGEND_FONTSIZE = 8

def draw_f3(db):
	print '--- Drawing Figure 3... ---'

	fig = plt.figure(1)
	f3_data = db.hw.stats.find_one({'_id': 'f3_data'})
	subplot_num = 1
	draw_f3_variant_pop_distribution(subplot_num, f3_data)
	subplot_num += 1
	draw_f3_variant_csqs(subplot_num, f3_data)
	subplot_num += 1
	draw_f3_clin_var(subplot_num, f3_data)
	subplot_num += 1
	draw_f3_ad_ar_stats(subplot_num, f3_data)

	report_f3_pop_stats(f3_data)
	report_f3_variants_csq_stats(f3_data)

	plt.figure(figsize = (2,2))
	gs1 = gridspec.GridSpec(2, 2)
	gs1.update(wspace=0.05, hspace=0.05)
	fig = plt.figure(1)
	fig.set_size_inches(7, 6)
	plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.98])
	plt.savefig(FIGURES_FOLDER + 'F3.png', format='png', dpi=300)
	plt.close(fig)


def report_f3_pop_stats(f3_data):
	he_unique_vars = f3_data['het_excess_unique_variants']
	nhe_unique_vars = f3_data['not_het_excess_unique_variants']
	print 'HetExc Unique', he_unique_vars
	print 'HetExc 7 pop', sum(f3_data['pop_het_excess_variants'].values())
	print 'HetExc genes', f3_data['het_excess_gene_stats']['total_num']
	print 'HetExc- Unique', nhe_unique_vars
	print 'HetExc- 7 pop', sum(f3_data['pop_not_het_excess_variants'].values())
	print 'HetExc- genes', f3_data['not_het_excess_gene_stats']['total_num']

	print 'Pop enrichemnt tests (Fisher Exact)'
	for pop, he_vars in f3_data['pop_het_excess_variants'].iteritems():
		he_data = [he_vars, he_unique_vars]
		nhe_vars = f3_data['pop_not_het_excess_variants'][pop]
		nhe_data = [nhe_vars, nhe_unique_vars]
		fe, p_value = fisher_exact([he_data, nhe_data])
		print '{} {} out of {} in HetExc vs {} out of {} in HetExc-'.format(pop, he_vars, he_unique_vars, nhe_vars, nhe_unique_vars)
		print '{} FE:{} p-value: {}'.format(pop, fe, p_value)
		table_name = 'f3_a_' + pop + '.csv'
		title = 'Comparison of proportions of variants deviating and not deviating from HWE due to excess of heterozygotes (HetExc and HetExc- respectively) in 7 major gnomAD populations.'
		group_names = ['HetExc', 'HetExc-']
		variable_names = [pop, 'all populations']
		export_fisher_test_results(table_name, title, group_names, variable_names, he_data, nhe_data, fe, p_value)


def draw_f3_clin_var(subplot_num, f3_data, title='', annotation_type='', xticks=[], x_offset=0):
	ax = plt.subplot(2,2,subplot_num)
	ax.text(-0.15, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')

	clin_var_stats = f3_data['pop_het_excess_clinvar']

	xs = clin_var_stats.values()
	colours = [C_RED, C_ORANGE, C_GREEN, C_LIGHT_GRAY]

	ys = ['Pathogenic (P)',
			 'Conflicting\nInterpretations\nof Pathogenicity (CIP)',
			 'Benign or\nLikely Benign (BLB)',
			 'Unknown (U)']

	width = 0.8 # the width of the bars 
	ind = np.arange(len(ys))  # the y locations for the groups
	ax.bar(ind, xs, width, color=colours)

	y_text_offset = max(xs) * 0.02
	for i, v in enumerate(xs):
		
		ax.text(i-0.05, v+y_text_offset, " {:,}".format(v), va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE)

	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)

	ax.set_xticks(ind)
	ax.set_xticklabels(['P', 'CIP', 'BLB', 'U'], minor=False, fontsize=F3_AXIS_FONTSIZE) # , rotation=15
	ax.set_ylim(ymax=max(xs) + x_offset)

	ax.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=True) # labels along the bottom edge are off

	if xticks:
		ax.set_xticks(xticks)

	if annotation_type == 'int':
		ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
	
	ax.set_title(title, fontsize=12)
	ax.set_ylabel('Number of variants (HetExc)', fontsize=F3_AXIS_FONTSIZE)
	ax.tick_params(axis='y', which='major', labelsize=F3_AXIS_FONTSIZE)

	# Adds CFTR and HBB examples to the plot
	#ax.text(0, 40, "CFTR\np.Phe508del", va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE, rotation=25)
	#ax.text(0, 75, "HBB\np.Glu7Val", va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE, rotation=25)

	clinvar_legend = OrderedDict()
	for x in range(0, len(ys)):
		clinvar_legend[ys[x]] = mpatches.Patch(color=colours[x])
	plt.legend(clinvar_legend.values(), clinvar_legend.keys(), loc='upper left', frameon=False, ncol=1, fontsize=F3_LEGEND_FONTSIZE)


def draw_f3_ad_ar_stats(subplot_num, f3_data, title='', annotation_type='', xticks=[], x_offset=0, combined_ar=True):
	ax = plt.subplot(2,2,subplot_num)
	ax.text(-0.16, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')

	nhe_gene_stats = f3_data['not_het_excess_gene_stats']
	he_gene_stats = f3_data['het_excess_gene_stats']
	nhe_total = nhe_gene_stats['total_num']
	he_total = he_gene_stats['total_num']

	if combined_ar:
		x_labels = ['AD', 'AR or AR,AD']
		nhe_num_ys = [nhe_gene_stats['ad_num'], nhe_gene_stats['ar_num'] + nhe_gene_stats['ar_ad_num']]
		nhe_percent_ys = [nhe_gene_stats['ad_percent'], nhe_gene_stats['ar_percent'] + nhe_gene_stats['ar_ad_percent']]
		he_num_ys = [he_gene_stats['ad_num'], he_gene_stats['ar_num'] + he_gene_stats['ar_ad_num']]
		he_percent_ys = [he_gene_stats['ad_percent'], he_gene_stats['ar_percent'] + he_gene_stats['ar_ad_percent']]
	else:
		x_labels = ['AD', 'AR', 'AR,AD']
		nhe_num_ys = [nhe_gene_stats['ad_num'], nhe_gene_stats['ar_num'], nhe_gene_stats['ar_ad_num']]
		nhe_percent_ys = [nhe_gene_stats['ad_percent'], nhe_gene_stats['ar_percent'], nhe_gene_stats['ar_ad_percent']]
		he_num_ys = [he_gene_stats['ad_num'], he_gene_stats['ar_num'], he_gene_stats['ar_ad_num']]
		he_percent_ys = [he_gene_stats['ad_percent'], he_gene_stats['ar_percent'], he_gene_stats['ar_ad_percent']]
	
	
	xs = np.arange(0, len(x_labels))
	bar_width, x_paddings = bar_chart_get_bar_width_and_x_paddings(2)

	ax.bar(xs + x_paddings[0], nhe_percent_ys, width=bar_width, label='HetExc-')#label='HetExcess-\n({:,})'.format(nhe_total))
	ax.bar(xs + x_paddings[1], he_percent_ys, width=bar_width, label='HetExc')#label='HetExcess\n({:,})'.format(he_total))

	# Percentages
	y_text_offset = max(nhe_percent_ys + he_percent_ys) * 0.02
	for i, v in enumerate(nhe_percent_ys):
		ax.text(i + x_paddings[0], v+y_text_offset, "{0:.1f}%".format(v), va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE)

		ax.text(i + x_paddings[0], v / 2.0, "{:,}".format(nhe_num_ys[i]), va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE, color='white')

	for i, v in enumerate(he_percent_ys):
		ax.text(i + x_paddings[1], v+y_text_offset, "{0:.1f}%".format(v), va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE)

		ax.text(i + x_paddings[1], v / 2.0, "{:,}".format(he_num_ys[i]), va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE)

	if combined_ar:
		nhe_ar_stats = [nhe_gene_stats['ar_num'] + nhe_gene_stats['ar_ad_num'], nhe_total]
		he_ar_stats = [he_gene_stats['ar_num'] + he_gene_stats['ar_ad_num'], he_total]
	else:	
		nhe_ar_stats = [nhe_gene_stats['ar_num'], nhe_total]
		he_ar_stats = [he_gene_stats['ar_num'], he_total]

	ar_fe, ar_p_value = fisher_exact([he_ar_stats, nhe_ar_stats])
	print 'AR HE/NHE Fihser:', ar_fe, ar_p_value, he_ar_stats, nhe_ar_stats
	title = 'Comparison of proportions of AR or AR,AD and all genes with at least one variant in HetExc and HetExc- datasets'
	export_fisher_test_results('f3_d_1.csv', title, ['HetExc', 'HetExc-'], ['AR or AR,AD', 'All'], 
							   he_ar_stats, nhe_ar_stats, ar_fe, ar_p_value)

	nhe_ad_stats = [nhe_gene_stats['ad_num'], nhe_total]
	he_ad_stats = [he_gene_stats['ad_num'], he_total]

	ad_fe, ad_p_value = fisher_exact([he_ad_stats, nhe_ad_stats])
	print 'AD HE/NHE Fihser:', ad_fe, ad_p_value, he_ad_stats, nhe_ad_stats
	title = 'Comparison of proportions of AD and all genes with at least one variant in HetExc and HetExc- datasets'
	export_fisher_test_results('f3_d_2.csv', title, ['HetExc', 'HetExc-'], ['AD', 'All'], 
							   he_ad_stats, nhe_ad_stats, ad_fe, ad_p_value)
	center = []
	height = []

	for x in range(0, len(xs)):
		center.append(xs[x] + x_paddings[0])
		center.append(xs[x] + x_paddings[1])
		height.append(nhe_percent_ys[x])
		height.append(he_percent_ys[x])


	barplot_annotate_brackets(0, 1, ad_p_value, center, height, dh=.06, barh=.03, fs=F3_NORMAL_FONTSIZE)
	barplot_annotate_brackets(2, 3, ar_p_value, center, height, dh=.06, barh=.03, fs=F3_NORMAL_FONTSIZE)
	if not combined_ar:
		nhe_ar_ad_stats = [nhe_gene_stats['ar_ad_num'], nhe_total]
		he_ar_ad_stats = [he_gene_stats['ar_ad_num'], he_total]
		ar_ad_fe, ar_ad_p_value = fisher_exact([he_ar_ad_stats, nhe_ar_ad_stats])
		barplot_annotate_brackets(4, 5, ar_ad_p_value, center, height, dh=.06, barh=.03, fs=F3_NORMAL_FONTSIZE)	
	
	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)

	if combined_ar:
		plt.legend(loc='upper left', frameon=False, ncol=1, fontsize=F3_LEGEND_FONTSIZE)
	else:
		plt.legend(loc=[0.01,0.7], frameon=False, ncol=1, fontsize=F3_LEGEND_FONTSIZE)

	ax.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=True) # labels along the bottom edge are off

	ax.set_ylabel('Genes (%)', fontsize=F3_AXIS_FONTSIZE)
	ax.tick_params(axis='y', which='major', labelsize=F3_AXIS_FONTSIZE)
	ax.set_xticks(xs)
	ax.set_xticklabels(x_labels, minor=False, fontsize=F3_AXIS_FONTSIZE)


def draw_f3_variant_pop_distribution(subplot_num, f3_data, title='', annotation_type='', xticks=[], x_offset=0):
	ax = plt.subplot(2,2,subplot_num)
	ax.text(-0.15, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')

	groups = ['HetExc-', 'HetExc']
	pop_nhe_variants = f3_data['pop_not_het_excess_variants']
	pop_he_variants = f3_data['pop_het_excess_variants']
	nhe_unique_variants = f3_data['not_het_excess_unique_variants']
	he_unique_variants = f3_data['het_excess_unique_variants']

	nhe_pop_variants_total = sum(pop_nhe_variants.values())
	he_pop_variants_total = sum(pop_he_variants.values())

	pop_variants_proportions = OrderedDict()
	pop_variants_nums = OrderedDict() 

	reverse_pop_order = POP_ORDER[:]
	reverse_pop_order.reverse()

	chi_sqr_nhe_variants = []
	chi_sqr_he_variants = []
	for pop in reverse_pop_order:
		pop_variants_proportions[pop] = []
		pop_variants_proportions[pop].append(pop_nhe_variants[pop] * 100 / float(nhe_pop_variants_total))
		pop_variants_proportions[pop].append(pop_he_variants[pop] * 100 / float(he_pop_variants_total))
		pop_variants_proportions[pop].append(0) # extra for space for figure legend

		pop_variants_nums[pop] = []
		pop_variants_nums[pop].append(pop_nhe_variants[pop])
		pop_variants_nums[pop].append(pop_he_variants[pop])
		chi_sqr_nhe_variants.append(pop_nhe_variants[pop])
		chi_sqr_he_variants.append(pop_he_variants[pop])

	chi_sqr_df = len(chi_sqr_nhe_variants) - 1
	bottom = [0] * len(groups) + [0] # extra for space for figure legend
	xs = range(0, len(groups) + 1) # extra for space for figure legend

	bar_width = 0.85

	for pop, group_proportions in pop_variants_proportions.iteritems():
		plt.bar(xs, group_proportions, bottom=bottom, color=POP_COLOURS[pop], width=bar_width)
		pop_variants = pop_variants_nums[pop]
		for x in range(0, len(bottom)):
			pop_proportion = group_proportions[x]
			if pop_proportion > 3:
				ax.text(x, (bottom[x] + pop_proportion / 2.3), "{:,}".format(pop_variants[x]), va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE, color='white')

			bottom[x] += pop_proportion

	ax.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=True) # labels along the bottom edge are off

	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)

	#ax.set_xlim(2)

	pop_legend = OrderedDict()
	for pop in POP_ORDER:
		pop_legend[pop] = mpatches.Patch(color=POP_COLOURS[pop])
	
	plt.legend(pop_legend.values(), pop_legend.keys(), loc='upper right', frameon=False, ncol=1, fontsize=F3_LEGEND_FONTSIZE)

	ax.set_ylabel('Variants (%)', fontsize=F3_AXIS_FONTSIZE)
	ax.tick_params(axis='y', which='major', labelsize=F3_AXIS_FONTSIZE)
	ax.set_xticklabels(['', 'HetExc-', 'HetExc'], fontsize=F3_AXIS_FONTSIZE)


def count_f3_variant_csq_groups(csq_dict):
	lof_csqs = set(['splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant'])

	lof = 0
	miss = 0
	syn = 0
	oth = 0
	all_variants = 0.0
	for csq, variant_num in csq_dict.iteritems():
		all_variants += variant_num
		if csq in lof_csqs:
			lof += variant_num
		elif csq == 'missense_variant':
			miss += variant_num
		elif csq == 'synonymous_variant':
			syn += variant_num
		else:
			oth += variant_num

	csq_group_num = [miss, syn, oth + lof]
	csq_group_percent = []
	for n in csq_group_num:
		csq_group_percent.append(n * 100 / all_variants)
	return csq_group_num, csq_group_percent


def draw_f3_variant_csqs(subplot_num, f3_data):
	ax = plt.subplot(2,2,subplot_num)
	ax.text(-0.16, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')

	x_labels = ['Missense', 'Synonymous', 'Other'] # 'Loss-of-Function', 
	nhe_num_ys, nhe_percent_ys = count_f3_variant_csq_groups(f3_data['not_het_excess_csqs'])
	he_num_ys, he_percent_ys = count_f3_variant_csq_groups(f3_data['het_excess_csqs'])

	xs = np.arange(0, len(x_labels))
	bar_width, x_paddings = bar_chart_get_bar_width_and_x_paddings(2, base_width=0.9)

	ax.bar(xs + x_paddings[0], nhe_percent_ys, width=bar_width, label='HetExc-')
	ax.bar(xs + x_paddings[1], he_percent_ys, width=bar_width, label='HetExc')

	# Percentages
	y_text_offset = max(nhe_percent_ys + he_percent_ys) * 0.02
	for i, v in enumerate(nhe_percent_ys):
		ax.text(i + x_paddings[0], v+y_text_offset, "{0:.1f}%".format(v), va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE)
		ax.text(i + x_paddings[0], v / 2.0, "{:,}".format(nhe_num_ys[i]), va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE, color='white')

	for i, v in enumerate(he_percent_ys):
		ax.text(i + x_paddings[1], v+y_text_offset, "{0:.1f}%".format(v), va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE)
		ax.text(i + x_paddings[1], v / 2.0, "{:,}".format(he_num_ys[i]), va='center', ha='center', fontsize=F3_NORMAL_FONTSIZE)


	center = []
	height = []

	for x in range(0, len(xs)):
		center.append(xs[x] + x_paddings[0])
		center.append(xs[x] + x_paddings[1])
		height.append(nhe_percent_ys[x])
		height.append(he_percent_ys[x])

	total_nhe = sum(nhe_num_ys)
	total_he = sum(he_num_ys)

	bar_1_index = 0
	bar_2_index = 1
	for x in range(0, len(nhe_num_ys)):
		csq_type = x_labels[x]
		nhe_csq_group = [nhe_num_ys[x], total_nhe]
		he_csq_group = [he_num_ys[x], total_he]
		csq_group_fe, csq_group_p_value = fisher_exact([nhe_csq_group, he_csq_group])
		title = 'Comparison of ' + csq_type + ' variants proportions in HetExc and HetExc- datasets'
		export_fisher_test_results('f3_b_' + str(x) + '.csv', title, ['HetExc', 'HetExc-'], [csq_type, 'All'], 
								   he_csq_group, nhe_csq_group, csq_group_fe, csq_group_p_value)

		barplot_annotate_brackets(bar_1_index, bar_2_index, csq_group_p_value, center, height, dh=.06, barh=.03, fs=F3_NORMAL_FONTSIZE)
		bar_1_index += 2
		bar_2_index += 2

	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)

	plt.legend(loc='upper right', frameon=False, ncol=1, fontsize=F3_LEGEND_FONTSIZE)

	ax.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=True) # labels along the bottom edge are off

	ax.set_ylabel('Variants (%)', fontsize=F3_AXIS_FONTSIZE)
	ax.tick_params(axis='y', which='major', labelsize=F3_AXIS_FONTSIZE)
	ax.set_xticks(xs)
	ax.set_xticklabels(x_labels, minor=False, fontsize=F3_AXIS_FONTSIZE)


def report_f3_variants_csq_stats(f3_data):
	he_miss = f3_data['het_excess_csqs']['missense_variant']
	he_syn = f3_data['het_excess_csqs']['synonymous_variant']
	he_all = f3_data['het_excess_unique_variants']

	nhe_miss = f3_data['not_het_excess_csqs']['missense_variant']
	nhe_syn = f3_data['not_het_excess_csqs']['synonymous_variant']
	nhe_all = f3_data['not_het_excess_unique_variants']

	print 'HetExc-'
	print 'Missense {} out of {}, {:.1f}%'.format(nhe_miss, nhe_all, float(nhe_miss) * 100 / nhe_all)
	print 'Synonymous {} out of {}, {:.1f}%'.format(nhe_syn, nhe_all, float(nhe_syn) * 100 / nhe_all)
	print 'HetExc'
	print 'Missense {} out of {}, {:.1f}%'.format(he_miss, he_all, float(he_miss) * 100 / he_all)
	print 'Synonymous {} out of {}, {:.1f}%'.format(he_syn, he_all, float(he_syn) * 100 / he_all)

################
### FIGURE 4 ###
################

def caculate_f4_data(db):
	print '--- Calculating data for Figure 4... ---'

	for pop in POP_ORDER:
		print pop
		variants = db.hw.variants_hwe_pop.find({ "pop": pop, "all_pop_an_ratio_pass": True, "af": { "$lte": RARE_HET_EXCESS_MAX_AF }, "alt_af": { "$lt": 0.001 }})

		het_exc_afs = []
		het_exc_hom_afs = []

		het_def_afs = []
		het_def_hom_afs = []

		normal_afs = []
		normal_hom_afs = []

		for variant in variants:
			af = variant['af']
			hom_af = float(variant['hom'] * 2) / variant['an']

			if variant['hwe_p_value'] < 0.05:
				if variant['balance'] == '-':
					het_exc_afs.append(af)
					het_exc_hom_afs.append(hom_af)
				else:
					het_def_afs.append(af)
					het_def_hom_afs.append(hom_af)
			else:
				normal_afs.append(af)
				normal_hom_afs.append(hom_af)

		# hwe 
		pop_ans = POP_SIZES[pop] * 2
		hwe_afs = []
		hwe_hom_afs = []
		for x in range(0, pop_ans):
			af = x / float(pop_ans)
			hom_af = af * af
			hwe_afs.append(af)
			hwe_hom_afs.append(hom_af)

			if af > RARE_HET_EXCESS_MAX_AF:
				break

		f4_pop_data = OrderedDict()
		f4_pop_data['_id'] = 'f4_data_' + pop
		f4_pop_data['het_exc_afs'] = het_exc_afs
		f4_pop_data['het_exc_hom_afs'] = het_exc_hom_afs
		f4_pop_data['het_def_afs'] = het_def_afs
		f4_pop_data['het_def_hom_afs'] = het_def_hom_afs
		f4_pop_data['normal_afs'] = normal_afs
		f4_pop_data['normal_hom_afs'] = normal_hom_afs
		f4_pop_data['hwe_afs'] = hwe_afs
		f4_pop_data['hwe_hom_afs'] = hwe_hom_afs
		db.hw.stats.delete_one({'_id': f4_pop_data['_id']})
		db.hw.stats.insert(f4_pop_data)


F4_SUBPLOT_LABEL_FONTSIZE = 12
F4_SCATTER_POINT_SIZE = 1
F4_NORMAL_FONTSIZE = 8
F4_LEGEND_FONTSIZE = 10
F4_AXIS_FONTSIZE = 10
F4_TITLE_FONTSIZE = 10

def draw_pop_hwe_scatter(db, subplot_num, pop):
	ax = plt.subplot(4,2,subplot_num)
	data = db.hw.stats.find_one({'_id': 'f4_data_' + pop})

	plt.scatter(data['het_def_afs'], data['het_def_hom_afs'], s=F4_SCATTER_POINT_SIZE, color=C_ORANGE)
	plt.scatter(data['het_exc_afs'], data['het_exc_hom_afs'], s=F4_SCATTER_POINT_SIZE, color=C_BLUE)
	plt.scatter(data['normal_afs'], data['normal_hom_afs'], s=F4_SCATTER_POINT_SIZE, color=C_LIGHT_GRAY)
	plt.scatter(data['hwe_afs'], data['hwe_hom_afs'], s=0.5, color='black')

	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	subplot_num_left_padding = -0.15
	if RARE_HET_EXCESS_MAX_AF == 0.1:
		ax.set_yticks([0.01, 0.02])
		ymax = 0.02
		ymin = -0.001
	elif RARE_HET_EXCESS_MAX_AF == 0.05:
		ax.set_yticks([0.0025, 0.005])
		ymax = 0.005
		ymin = -0.0003
		subplot_num_left_padding = -0.24
	else:
		ymax = RARE_HET_EXCESS_MAX_AF / 5.0
		ymin = -0.001

	ax.set_ylim(ymin=ymin, ymax=ymax)
	ax.text(0, ymax * 0.87, pop, fontsize=F4_TITLE_FONTSIZE, fontweight='bold')
	ax.text(subplot_num_left_padding, 1.15, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=F4_SUBPLOT_LABEL_FONTSIZE, fontweight='bold', va='top', ha='right')
	ax.tick_params(axis='both', which='major', labelsize=F4_NORMAL_FONTSIZE)


def draw_f4(db):
	print '--- Drawing Figure 4... ---'

	fig = plt.figure(1)
	subplot_num = 1
	for pop in POP_ORDER:
		draw_pop_hwe_scatter(db, subplot_num, pop)
		subplot_num += 1

	legend = OrderedDict()
	legend['Hardy-Weinberg Equilibrium'] = mpatches.Patch(color='black')
	legend['In Equilibrium'] = mpatches.Patch(color=C_LIGHT_GRAY)
	legend['Heterozygote Deficiency (HetDef)'] = mpatches.Patch(color=C_ORANGE)
	legend['Heterozygote Excess (HetExc)'] = mpatches.Patch(color=C_BLUE)
	fig.legend(legend.values(), legend.keys(), loc=[0.53,0.1], frameon=False, ncol=1, fontsize=F4_LEGEND_FONTSIZE)

	fig.text(0.5, 0.02, 'Allele Frequency (AF)', ha='center', va='center', fontsize=F4_AXIS_FONTSIZE)
	fig.text(0.03, 0.5, 'Homozygous Allele Frequency (Hom AF)', ha='center', va='center', rotation='vertical', fontsize=F4_AXIS_FONTSIZE)

	fig = plt.figure(1)
	fig.set_size_inches(7, 8)
	plt.tight_layout(rect=[0.06, 0.02, 0.99, 0.99], h_pad=1.7, w_pad=1.7)
	plt.savefig(FIGURES_FOLDER + 'F4.png', format='png', dpi=300)
	plt.close(fig)


##############
### TABLES ###
##############

def create_supplementary_table_het_exc_genes(db, all_genes=False):
	if all_genes:
		selected_genes = set()
		variants = db.hw.rare_het_excess_variants.find({ "ab_adjusted_hwe_het_exc": True })
		for variant in variants:
			selected_genes.add(variant['gdit_gene_name'])
	else:
		selected_genes = ['HBB',
						  'CFTR',
						  'PCSK9',
						  'LPL',
						  'CUBN',
						  'HSPG2',
						  'CETP',
						  'ABCA1',
						  'CHD4',
						  'CHD6',
						  'FRAS1',
						  'FREM2',
						 ]

	column_names = ['variant_id',
					'rsid',
					'csq',
					'canonical_transcript',
					'gdit_gene_name',
					'hgvsp',
					'hgvsc',
					'pop',	
					'ac',
					'an',
					'an_ratio',
					'af',
					'het',
					'exp_het',
					'hom',
					'exp_hom',
					'hwe_p_value',
					'ab_mid_ratio',
					'omim',	
					'Inheritance_pattern',	
					'phenotype',
					'clinvar',
					'clinvar_significance',
					'possible_hom_ab>09',
					'possible_hom_ab>08',

					'het_ac_proportion',
					'ab_adjusted_het',
					'ab_adjusted_hom',
					'ab_adjusted_hwe_p_value',
					'ab_adjusted_balance',
					]

	# Add gnomAD 3 data
	column_names += [
						'gnomad_v3_variant_id',
						'gnomad_v3_ac',
						'gnomad_v3_an',
						'gnomad_v3_an_ratio',
						'gnomad_v3_af',
						'gnomad_v3_het',
						'gnomad_v3_exp_het',
						'gnomad_v3_hom',
						'gnomad_v3_exp_hom',
						'gnomad_v3_hwe_p_value',
						'gnomad_v3_balance',
						'gnomad_v3_filters',
						'gnomad_v3_ab_mid_ratio',
						'gnomad_v3_possible_hom_ab>09',
						'gnomad_v3_possible_hom_ab>08',
						'gnomad_v3_het_ac_proportion',
						'gnomad_v3_ab_adjusted_het',
						'gnomad_v3_ab_adjusted_hom',
						'gnomad_v3_ab_adjusted_hwe_p_value',
						'gnomad_v3_ab_adjusted_balance',
					]


	table = [column_names]
	for gene_name in selected_genes:
		variants = db.hw.rare_het_excess_variants.find({'gdit_gene_name': gene_name})
		for variant in variants:
			row = []
			for column_name in column_names:
				if column_name in ['hgvsp', 'hgvsc']:
					gnomad_variant = db.hw.gnomad_variants.find_one({'variant_id': variant['variant_id']})
					row.append(gnomad_variant[column_name])
				else:
					row.append(variant[column_name])
			table.append(row)
	if all_genes:
		table_name = 'st1.csv'
	else:
		table_name = 't1.csv'
	output_csv = OUTPUT_FOLDER + table_name
	write_table_to_csv(table, output_csv)


def create_supplementary_table_1000g(db):
	hbb_variant = '11-5248232-T-A' # HBB
	chd6_variant = '20-40040825-C-G' # CHD6

	all_afr_individuals = db.g1000.individuals.find({ "Super Population": "AFR" })
	afr_individuals = set([])
	asw_and_acb_individuals = set([])

	aa_pops = ['Gambian in Western Divisions in the Gambia (GWD)', 
			   'Mende in Sierra Leone (MSL)', 
			   'Yoruba in Ibadan, Nigeria (YRI)', 
			   'Luhya in Webuye, Kenya (LWK)', 
			   'Esan in Nigeria (ESN)']
	am_pops = ['Americans of African Ancestry in SW USA (ASW)', 'African Caribbeans in Barbados (ACB)']

	for individual in all_afr_individuals:
		if individual['Population'] == 'ASW' or individual['Population'] == 'ACB':
			asw_and_acb_individuals.add(individual['_id'])
		else:
			afr_individuals.add(individual['_id'])

	variant_ids = [hbb_variant, chd6_variant]
	# gnomad 3 survived AFR variants
	'''
	variant_ids = [
					'11-5248232-T-A',
					'20-40040825-C-G',
					'11-126075608-G-A',
					'22-42610871-G-A',
					'1-22166484-C-T',
					'14-92471202-T-C',
					'20-33862165-G-A',
				   ]
	'''
	headers = ['Variant ID', 
			   'African (AFR) populations', 
			   'African American (AFR-AMR) populations',
			   'AFR Allele Count (AC)',
			   'AFR Allele Number (AN)',
			   'AFR Allele Frequency (AF)',
			   'AFR-AMR AC',
			   'AFR-AMR AN',
			   'AFR-AMR AF',
			   'AFR/AFR-AMR fold enrichement',
			   'AFR/AFR-AMR p-value',
			   'AFR individual IDs (Het)',
			   'AFR individual IDs (Hom)',
			   'AFR-AMR individual IDs (Het)',
			   'AFR-AMR individual IDs (Hom)',
			  ]
	table = [headers]

	for variant_id in variant_ids:
		aa_het_individuals = []
		aa_hom_individuals = []
		am_het_individuals = []
		am_hom_individuals = []

		aa_ac = 0
		am_ac = 0

		var_id_1kg = variant_id.replace('-','_')
		var_1kg = db.g1000.coding_variants.find_one({"_id": var_id_1kg})

		for individual_id, genotype in var_1kg['INDIVIDUALS'].iteritems():
			ac_num = 1
			if genotype == 3: # if homozygous
				ac_num = 2
			if individual_id in afr_individuals:
				if genotype == 3:
					aa_hom_individuals.append(individual_id)
				else:
					aa_het_individuals.append(individual_id)
				aa_ac += ac_num
			elif individual_id in asw_and_acb_individuals:
				if genotype == 3:
					am_hom_individuals.append(individual_id)
				else:
					am_het_individuals.append(individual_id)
				am_ac += ac_num

		aa_an = len(afr_individuals) * 2
		am_an = len(asw_and_acb_individuals) * 2

		aa_af = float(aa_ac) / aa_an
		am_af = float(am_ac) / am_an

		fe, p_value = fisher_exact([[aa_ac, aa_an],[am_ac, am_an]])

		row = [variant_id,
			   ', '.join(aa_pops),
			   ', '.join(am_pops),
			   aa_ac,
			   aa_an,
			   aa_af,
			   am_ac,
			   am_an,
			   am_af,
			   fe,
			   p_value,
			   ', '.join(aa_het_individuals),
			   ', '.join(aa_hom_individuals),
			   ', '.join(am_het_individuals),
			   ', '.join(am_hom_individuals),
		]
		table.append(row)

	output_csv = OUTPUT_FOLDER + 'st2.csv'
	write_table_to_csv(table, output_csv)


###################
### PAPER STATS ###
###################

def print_initial_variant_dataset_stats(db):
	print 'Initial dataset stats'
	variants = db.hw.variants_hwe_pop.find({"all_pop_an_ratio_pass": True, "alt_af": { "$lt": 0.001 }})

	total_variants = variants.count()
	unique_variants = set([])
	unique_genes = set([])
	extreme_het_exc_variants = set([])
	extreme_het_exc_variants_in_repeat = set([])
	miss_variants = set([])
	syn_variants = set([])

	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for variant in variants:
		unique_variants.add(variant['variant_id'])
		unique_genes.add(variant['gdit_gene_name'])

		if variant['af'] >= 0.23 and variant['hom'] == 0:
			extreme_het_exc_variants.add(variant['variant_id'])
			if variant['tandem_repeat'] or variant['seg_dup']:
				extreme_het_exc_variants_in_repeat.add(variant['variant_id'])

		if variant['csq'] == 'missense_variant':
			miss_variants.add(variant['variant_id'])
		if variant['csq'] == 'synonymous_variant':
			syn_variants.add(variant['variant_id'])

		line_number += 1
		bar.update((line_number + 0.0) / total_variants)
	bar.finish()

	print 'Total variants:', total_variants
	print 'Unique variants:', len(unique_variants)
	print 'Missense {}, {:.1f}%'.format(len(miss_variants), float(len(miss_variants)) / len(unique_variants))
	print 'Synonymous {}, {:.1f}%'.format(len(syn_variants), float(len(syn_variants)) / len(unique_variants))
	print 'Unique genes:', len(unique_genes)
	print 'For comparison with Graffelman study:'
	print 'Extreme HetExc (AF>=0.23; hom=0)', len(extreme_het_exc_variants)
	print 'Extreme HetExc (AF>=0.23; hom=0; seg dup or tandem repeat)', len(extreme_het_exc_variants_in_repeat)


def report_initial_and_future_dataset_pop_limits(db, pop):
	gnomad_af_threshold = calculate_single_pop_significant_af_thresholds_data(POP_SIZES[pop], 0)['af']
	m1_af_threshold = calculate_single_pop_significant_af_thresholds_data(1000000, 0)['af']
	m5_af_threshold = calculate_single_pop_significant_af_thresholds_data(5000000, 0)['af']

	print pop
	print "GnomAD AF limit:", gnomad_af_threshold
	print "1 million AF limit:", m1_af_threshold
	print "5 million AF limit:", m5_af_threshold
	
	total_variants = set([])

	gnomad_detectable_variants = set([])
	m1_detectable_variants = set([])
	m5_detectable_variants = set([])
	variants = db.hw.variants_hwe_pop.find({"pop": pop, "all_pop_an_ratio_pass": True, "alt_af": { "$lt": 0.001 }, "af": {"$lte": RARE_HET_EXCESS_MAX_AF }}) # , "af": {"$lt": RARE_HET_EXCESS_MAX_AF }
	
	for variant in variants:
		variant_id = variant['variant_id']
		total_variants.add(variant_id)
		af = variant['af']
		if af >= gnomad_af_threshold:
			gnomad_detectable_variants.add(variant_id)
		if af >= m1_af_threshold:
			m1_detectable_variants.add(variant_id)
		if af >= m5_af_threshold:
			m5_detectable_variants.add(variant_id)

	print 'Total', len(total_variants)
	print 'gnomAD detectable {} ({:.2f}%)'.format(len(gnomad_detectable_variants), len(gnomad_detectable_variants) * 100 / float(len(total_variants)))
	print '1 Million detectable {} ({:.2f}%)'.format(len(m1_detectable_variants), len(m1_detectable_variants) * 100 / float(len(total_variants)))
	print '5 Million detectable {} ({:.2f}%)'.format(len(m5_detectable_variants), len(m5_detectable_variants) * 100 / float(len(total_variants)))


def print_final_dataset_stats(db, ab_adjusted_hwe_het_exc=False):
	if ab_adjusted_hwe_het_exc:
		print '### Final HetExc dataset stats ###'
		variants = db.hw.rare_het_excess_variants.find({'ab_adjusted_hwe_het_exc': True})
	else:
		print '### HetExc dataset stats before HWE recalculation with adjusted AB ###'
		variants = db.hw.rare_het_excess_variants.find({})

	total_variants = variants.count()
	unique_variants = set([])
	unique_genes = set([])
	miss_variants = set([])
	syn_variants = set([])
	afr_variants = set([])

	for variant in variants:
		variant_id = variant['variant_id']
		unique_variants.add(variant_id)
		unique_genes.add(variant['gdit_gene_name'])
		if variant['pop'] == 'AFR':
			afr_variants.add(variant_id)
		if variant['csq'] == 'missense_variant':
			miss_variants.add(variant_id)
		if variant['csq'] == 'synonymous_variant':
			syn_variants.add(variant_id)

	print 'Total variants:', total_variants
	print 'Unique variants:', len(unique_variants)
	print 'Unique genes:', len(unique_genes)
	print 'Missense {}, {:.2f}%'.format(len(miss_variants), float(len(miss_variants)) * 100 / len(unique_variants))
	print 'Synonymous {}, {:.2f}%'.format(len(syn_variants), float(len(syn_variants)) * 100 / len(unique_variants))
	print 'AFR variants {}, {:.2f}%'.format(len(afr_variants), float(len(afr_variants)) * 100 / len(unique_variants))


def print_graffelman_japanese_stats(db):
	print 'Variant stats required to claim statistical significant HetExc in Graffelman et al. 2017 study:'
	stats = calculate_single_pop_significant_af_thresholds_data(104, 0, p_value_threshold=0.001)
	for key, value in stats.iteritems():
		print key, value


def print_chen_stats(db):
	print 'ExAC rs61002819 stats'
	ac = 30929 	
	an = 60138 	
	rare_hom = 7609
	het = ac - rare_hom * 2
	common_hom = an / 2 - rare_hom - het
	p_value = hwe(het, rare_hom, common_hom, mid_p=True)
	print 'HWE p-value', p_value


def report_clinvar_statuses_and_csqs_of_het_exc_ar_variants(db):
	all_variants = set()
	all_clinvar = set()
	benign_variants = set()
	syn_variants = set()
	variants = db.hw.rare_het_excess_variants.find({"ab_adjusted_hwe_het_exc": True, "$or": [ { "Inheritance_pattern": "AR" }, { "Inheritance_pattern": "AR,AD" } ] })
	for variant in variants:
		variant_id = variant['variant_id']

		all_variants.add(variant_id)
		if variant['csq'] == 'synonymous_variant':
			syn_variants.add(variant_id)

		if variant['clinvar_significance'] in CLINVAR_BENIGN_STATUSES:
			benign_variants.add(variant_id)

		if variant['clinvar_significance'] != '':
			all_clinvar.add(variant_id)

	print 'AR HetExc variants'
	print 'All:', len(all_variants)
	print 'Benign:', len(benign_variants)
	print 'Any ClinVar:', len(all_clinvar)
	print 'Synonymous:', len(syn_variants)


def report_cftr_g2_vs_g3(db):
	variant = db.hw.rare_het_excess_variants.find_one({'variant_id': '7-117199644-ATCT-A'})
	g2_het = variant['het']
	g2_hom = variant['hom']
	g3_het = variant['het_g3']
	g3_hom = variant['hom_g3']
	print 'gnomAD v2.1 het:', g2_het, 'hom:', g2_hom
	print 'gnomAD v3 het:', g3_het, 'hom:', g3_hom
	fe, p_value = fisher_exact([[g2_het, g2_hom],[g3_het, g3_hom]])
	print 'fold enrichement', fe, 'p-value', p_value
	print 'with extra Hom'
	fe, p_value = fisher_exact([[g2_het, g2_hom],[g3_het - 2, g3_hom + 1]])
	print 'fold enrichement', fe, 'p-value', p_value


def report_number_of_unique_het_exc_variants_with_skewed_allele_balance(db, ab_group):
	print ab_group
	variants = db.hw.rare_het_excess_variants.find({ ab_group: { "$gt": 0 } })
	unique_variants = set()
	for variant in variants:
		unique_variants.add(variant['variant_id'])
	print 'Unique HetExc variants with skewed allele balance:', len(unique_variants)


def report_number_of_more_significant_het_exc_variants(db, p_value_threshold):
	unique_variants = set()
	variants = db.hw.rare_het_excess_variants.find({"hwe_p_value": { "$lte": p_value_threshold }})
	for variant in variants:
		unique_variants.add(variant['variant_id'])
	print 'HetExc variants with p-value <=' + str(p_value_threshold) + ':', len(unique_variants)


def main():
	db = MongoDB()
	# Uncomment the functions to recreate/report Hardy-Weinberg Equilibrium (HWE) analysis results.
	# Original figures and supplementary tables are stored in "figures" and "tables" folders, respectively.
	# Re-running the functions that draw figures or export supplementary tables will rewrite original versions.

	############################################################
	### Functions to create figures and supplementary tables ###
	############################################################

	# For the first run, functions to calculate data for the figures has to be uncommented.
	# After that, the figures source data is stored in the database and
	# only "draw" functions can be used to recreate the figures.

	#caculate_f1_data(db)
	#draw_f1(db) # Figure 1 in the manuscript

	#calculate_f2_data(db)
	#draw_f2(db) # Figure 3 in the manuscript

	#caculate_f3_data(db)
	#draw_f3(db) # Figure 4 in the manuscript

	#caculate_f4_data(db)
	#draw_f4(db) # Figure 2 in the manuscript

	# Supplementary Table 1: HetExc variant dataset
	#create_supplementary_table_het_exc_genes(db, all_genes=True)

	'''
	Recreation of Supplementary Table 2 requires 1000 Genomes data.
	We used custom version of 1000 genomes database which was created from the VCF files.
	This database contained 2504 individuals, which were analysed in the original 1000 Genomes Phase 3 study:
	"An integrated map of structural variation in 2,504 human genomes"
	https://www.nature.com/articles/nature15394
	The custom version of the database is not shared due to its large size. 
	However, we provided all individual ids used in this analysis, so these results can be checked
	in the 1000 genomes variant browser:
	http://grch37.ensembl.org/Homo_sapiens/Info/Index
	OR VCF files:
	ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
	'''
	#create_supplementary_table_1000g(db)

	###############################################################
	### Functions to calculate stats reported in the manuscript ###
	###############################################################

	# Uncomment the functions to see the stats reported in the terminal.

	# Report stats of the initial dataset.
	#print_initial_variant_dataset_stats(db)

	# Report stats of HetExc dataset.
	#print_final_dataset_stats(db, ab_adjusted_hwe_het_exc=False)
	#print_final_dataset_stats(db, ab_adjusted_hwe_het_exc=True)
	
	# Report stats used in discussion of the previous work.
	#print_chen_stats(db)
	#print_graffelman_japanese_stats(db)
	
	# Report proportion of synonymous variants among HetExc variants observed in known AR genes.
	#report_clinvar_statuses_and_csqs_of_het_exc_ar_variants(db)
	
	# Report number of unique HetExc variants with skewed allele balance (>0.9)
	# which might have more homozygoues and therefore their HWE statistics might be inacurate 
	#report_number_of_unique_het_exc_variants_with_skewed_allele_balance(db, 'het_with_skewed_allele_balance_>09')
	#report_number_of_unique_het_exc_variants_with_skewed_allele_balance(db, 'het_with_skewed_allele_balance_>08')

	# Calculate required population size to detect statistical significant HetExc of c.448G>C (rs1800546) variant in ALDOB based on its AF
	#print calculate_minimum_population_number_required_to_achive_significance_for_af(0.0049, p_value_threshold=HWE_P_VALUE_THRESHOLD)
	
	# Report current AF statistical significance thresholds of the largest gnomAD population (Non-Finnish European) and 
	# how it will change when the population will increate to 1 and 5 millions. 
	#report_initial_and_future_dataset_pop_limits(db, "NFE") # gnomAD pop, 1 million, 5 millions

	# Report increase of individuals with homozygous Cystic fibrosis causing variants in gnomAD v3 in comparison to v2.1.1.
	#report_cftr_g2_vs_g3(db)

	# Temporary function to check ClinVar statuses of HetExc variants, used to manage hard coded values in Figure 3.
	#report_clinvar_statuses_in_rare_het_excess_genes(db)

	# Report number of HetExc variants with p-value<=X
	#report_number_of_more_significant_het_exc_variants(db, 0.01)

if __name__ == "__main__":
	sys.exit(main())