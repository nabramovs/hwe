import csv
import pymongo
import numpy as np
from collections import OrderedDict
from decimal import Decimal


#################
### CONSTANTS ###
#################
DB_HOST = 'localhost'
DB_PORT = 27017
DB_NAME_HW = 'hw'
DB_NAME_1000G = '1000g'
DB_NAME_EXAC = 'exac'


class MongoDB():
	"""Database Client."""
	def __init__(self):
		client = pymongo.MongoClient(host=DB_HOST, port=DB_PORT, document_class=OrderedDict)
		self.hw = client[DB_NAME_HW]
		self.g1000 = client[DB_NAME_1000G]
		self.exac = client[DB_NAME_EXAC]


def write_table_to_csv(table, output_csv, delimiter=','):
	"""Write table (list of lists) to csv."""
	output_file = open(output_csv,'w+')
	writer = csv.writer(output_file, delimiter=delimiter)

	for row in table:
		writer.writerow(row)

	output_file.close()


def is_float(x):
	"""Check if value (e.g. string) can be converted to float."""
	try:
		a = float(x)
	except ValueError:
		return False
	else:
		return True


def is_int(x):
	"""Check if value (e.g. string) can be converted to integer."""
	try:
		a = float(x)
		b = int(a)
	except ValueError:
		return False
	else:
		return a == b

def get_confidence_intervals(stats, stats_name='', alpha=0.95):
	p = ((1.0-alpha)/2.0) * 100
	lower = max(0.0, np.percentile(stats, p))
	p = (alpha+((1.0-alpha)/2.0)) * 100
	upper = min(100, np.percentile(stats, p))
	print('%s %.2f confidence interval %.1f%% and %.1f%%' % (stats_name, alpha, lower, upper))
	return lower, upper


# Calculate length of a file
def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

# FROM GNOMAD BROWSER:

# Added to support xpos 
CHROMOSOMES = ['%s' % x for x in range(1, 23)]
CHROMOSOMES.extend(['X', 'Y', 'MT'])
CHROMOSOME_TO_CODE = { item: i+1 for i, item in enumerate(CHROMOSOMES) }

def get_single_location(chrom, pos):
	"""
	Gets a single location from chromosome and position
	chr must be actual chromosme code (chrY) and pos must be integer
	Borrowed from xbrowse
	"""
	return CHROMOSOME_TO_CODE[chrom] * int(1e9) + pos


def get_xpos(chrom, pos):
	"""
	Borrowed from xbrowse
	"""
	return get_single_location(chrom, int(pos))