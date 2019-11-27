import os
import sys
import progressbar
import pymongo
import csv
import time
from collections import OrderedDict
from common import MongoDB, write_table_to_csv, get_xpos, file_len
from csv_reader import CsvReader
from gnomad_api import get_gene_variants_by_gene_name, get_transcript_exons

# !!! IMPORTANT !!! unzip "source_data.zip" to "hwe" folder

SOURCE_DATA_FOLDER = './source_data/'
GENE_DISCOVERY_INFORMATION_TOOLKIT_LEK = SOURCE_DATA_FOLDER + 'gene_discovery_informatics_toolkit_lek.csv'
USCS_SEGMENTAL_DUP_TSV = SOURCE_DATA_FOLDER + 'segmental_dups.txt'
USCS_TANDEM_REPEATS_TSV = SOURCE_DATA_FOLDER + 'simple_tandem_repeats.txt'


def import_gdit_genes(db):
	gdit_genes = CsvReader(GENE_DISCOVERY_INFORMATION_TOOLKIT_LEK)
	gdit_genes.import_to_db(db.hw, 'gdit_genes')
	db.hw.gdit_genes.create_index([('gene', pymongo.ASCENDING)], name='gene_1')

############
### USCS ###
############

def import_seg_dups(db):
	db.hw.seg_dups.drop()
	segmental_dup = CsvReader(USCS_SEGMENTAL_DUP_TSV, delimiter='\t')

	filtered_data = []
	for document in segmental_dup.data:
		document['chrom'] = document['chrom'][3:]
		document['otherChrom'] = document['otherChrom'][3:]

		if len(document['chrom']) < 3: 
			xstart = get_xpos(document['chrom'], document['chromStart'])
			xend = get_xpos(document['chrom'], document['chromEnd'])
			document['xstart'] = xstart
			document['xend'] = xend
			filtered_data.append(document)

	segmental_dup.data = filtered_data
	segmental_dup.import_to_db(db.hw, 'seg_dups')

	db.hw.seg_dups.create_index([('chrom', pymongo.ASCENDING)], name='chrom_1')
	db.hw.seg_dups.create_index([('chromStart', pymongo.ASCENDING)], name='chromStart_1')
	db.hw.seg_dups.create_index([('chromEnd', pymongo.ASCENDING)], name='chromEnd_1')
	db.hw.seg_dups.create_index([('xstart', pymongo.ASCENDING)], name='xstart_1')
	db.hw.seg_dups.create_index([('xend', pymongo.ASCENDING)], name='xend_1')


def import_tandem_repeats(db):
	db.hw.tandem_repeats.drop()
	tandem_repeats = CsvReader(USCS_TANDEM_REPEATS_TSV, delimiter='\t')

	pos_repeats = set([])

	total_lines = len(tandem_repeats.data)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for document in tandem_repeats.data:
		# added for patch chroms like chr19_gl000209_random
		if len(document['chrom']) > 5:
			continue

		chrom = document['chrom'][3:]
		start = document['chromStart']
		end = document['chromEnd']
		
		xstart = get_xpos(chrom, start)
		xend = get_xpos(chrom, end)

		for xpos in range(xstart, xend + 1):
			pos_repeats.add(xpos)

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	total_lines = len(pos_repeats)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()

	bulk = db.hw.tandem_repeats.initialize_unordered_bulk_op()
	counter = 0
	for xpos in pos_repeats:
		bulk.insert({ '_id': xpos })
		counter += 1

		if (counter % 500 == 0 and counter > 0):
			bulk.execute()
			bulk = db.hw.tandem_repeats.initialize_ordered_bulk_op()

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)

	if (counter % 500 != 0):
		bulk.execute()
	bar.finish()


##################
### GNOMAD API ###
##################

def import_gnomad_variants(db, retry_on_missed_genes=False):
	if not retry_on_missed_genes:
		db.hw.gnomad_variants.drop()
		db.hw.clinvar_variants.drop()

	gene_names = []
	if not retry_on_missed_genes:	
		genes = db.hw.gdit_genes.find({})
		for gene in genes:
			gene_names.append(gene['gene'])
	else:
		genes = db.hw.stats.find_one({'_id': 'gnomad_api_errors'})
		gene_names = genes['gene_names']


	not_found_genes = []
	error_genes = []
	no_variants_in_gnomad_genes = []

	total_genes = len(gene_names)
	gene_num = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for gene_name in gene_names:
		try:
			gnomad_gene_variants = get_gene_variants_by_gene_name(gene_name, timeout=120)
			if gnomad_gene_variants:
				clinvar_variants = gnomad_gene_variants['clinvar_variants']
				if clinvar_variants:
					num = len(clinvar_variants)
					db.hw.clinvar_variants.insert({'_id': gene_name, 'num': num, 'variants': clinvar_variants})
				else:
					db.hw.clinvar_variants.insert({'_id': gene_name, 'num': 0, 'variants': []})

				variants_data = gnomad_gene_variants['variants']
				variants = []

				for variant_data in variants_data:
					variant = OrderedDict()
					variant['gene_name'] = gnomad_gene_variants['gene_name']
					variant['other_names'] = gnomad_gene_variants['other_names']
					variant['gdit_gene_name'] = gene_name
					variant['gene_id'] = gnomad_gene_variants['gene_id']
					variant['canonical_transcript'] = gnomad_gene_variants['canonical_transcript']
					variant['chrom'] = gnomad_gene_variants['chrom']

					for key, value in variant_data.iteritems():
						variant[key] = value

					variants.append(variant)

				if len(variants) > 0:
					db.hw.gnomad_variants.insert_many(variants)
				else:
					no_variants_in_gnomad_genes.append(gene_name)
			else:
				not_found_genes.append(gene_name)
		except:
			print gene_name
			error_genes.append(gene_name)

		time.sleep(2) # some time to wait between API calls

		gene_num += 1
		bar.update((gene_num + 0.0) / total_genes)
	bar.finish()

	db.hw.stats.remove({'_id': 'gnomad_api_errors'})
	db.hw.stats.insert({'_id': 'gnomad_api_errors', 'gene_names': error_genes})

	if not retry_on_missed_genes:
		db.hw.stats.remove({'_id': 'not_found_gdit_genes_in_gnomad'})
		db.hw.stats.insert({'_id': 'not_found_gdit_genes_in_gnomad', 'gene_names': not_found_genes})

		db.hw.stats.remove({'_id': 'no_variants_in_gnomad'})
		db.hw.stats.insert({'_id': 'no_variants_in_gnomad', 'gene_names': no_variants_in_gnomad_genes})

	db.hw.gnomad_variants.create_index([('variant_id', pymongo.ASCENDING)], name='variant_id_1')
	db.hw.gnomad_variants.create_index([('xpos', pymongo.ASCENDING)], name='xpos_1')


class Exon():
	def __init__(self):
		self.transcript_id = ''
		self.gene_id = ''
		self.feature_type = ''
		self.chrom = ''
		self.strand = '+'
		self.start = 0
		self.stop = 0

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['transcript_id'] = self.transcript_id
		dictionary['gene_id'] = self.gene_id
		dictionary['feature_type'] = self.feature_type
		dictionary['chrom'] = self.chrom
		dictionary['strand'] = self.strand
		dictionary['start'] = self.start
		dictionary['stop'] = self.stop
		dictionary['xstart'] = get_xpos(self.chrom, self.start)
		dictionary['xstop'] = get_xpos(self.chrom, self.stop)
		return dictionary


def import_gnomad_exons(db, retry_on_missed_transcripts=False):
	
	error_transcripts = []

	if not retry_on_missed_transcripts:
		db.hw.exons.drop()
		transcripts = db.hw.stats.find_one({'_id': 'transcript_ids_to_gene_names'})
		transcripts = transcripts['transcript_ids_to_gene_names'].keys()
	else:
		transcripts = db.hw.stats.find_one({'_id': 'gnomad_api_exons_errors'})
		transcripts = transcripts['transcript_ids']


	total_transcripts = len(transcripts)
	transcript_num = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for transcript_id in transcripts:
		transcript_num += 1
		try:
			transcript_data = get_transcript_exons(transcript_id)
			exons = []
			for exon_data in transcript_data['exons']:
				exon = Exon()
				exon.transcript_id = transcript_data['transcript_id']
				exon.gene_id = transcript_data['gene_id']
				exon.feature_type = exon_data['feature_type']
				exon.chrom = transcript_data['chrom']
				exon.strand = exon_data['strand']
				exon.start = exon_data['start']
				exon.stop = exon_data['stop']

				exons.append(exon.get_dictionary())

			db.hw.exons.insert_many(exons)
		except:
			print transcript_id
			error_transcripts.append(transcript_id)

		#time.sleep(1) # some time to wait between API calls

		bar.update((transcript_num + 0.0) / total_transcripts)
	bar.finish()

	db.hw.stats.remove({'_id': 'gnomad_api_exons_errors'})
	db.hw.stats.insert({'_id': 'gnomad_api_exons_errors', 'transcript_ids': error_transcripts})

	db.hw.exons.create_index([('transcript_id', pymongo.ASCENDING)], name='transcript_id_1')
	db.hw.exons.create_index([('feature_type', pymongo.ASCENDING)], name='feature_type_1')


def main():
	db = MongoDB()
	# For the first run, please carefully read the comments and uncomment all of the following functions:

	# Gene Discovery Informatics Toolkit
	# https://www.nature.com/articles/s41525-019-0081-z
	#import_gdit_genes(db)

	# USCS browser (table) simple tandem repeats and segmental duplications. (25-01-2018)
	# same references as in Graffelman Japanese paper.
	#import_seg_dups(db)
	#import_tandem_repeats(db)

	# !!! IMPORTANT !!! 
	# This function works with gnomAD API (i.e. requires internet connection) and takes ~15 hours to run
	# gnomAD API sometimes (rarely) fail to return data for valid genes
	# reruning analysis on those genes fixes the problem (might require several reruns)
	#import_gnomad_variants(db)
	#import_gnomad_variants(db, retry_on_missed_genes=True)

	# import gnomAD
	# same as with variants API sometimes fail, have to run twice
	#import_gnomad_exons(db)
	#import_gnomad_exons(db, retry_on_missed_transcripts=True)

if __name__ == "__main__":
	sys.exit(main())