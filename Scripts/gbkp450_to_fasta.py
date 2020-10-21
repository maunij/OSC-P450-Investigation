#Takes a .gbk and p450_list.csv as input, outputs a FASTA file of gene_id and translation

import sys
import csv
import Bio
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

fileGBK = sys.argv[1]
fileP450 = sys.argv[2]

def get_gene_id_CO(feature):
	"Get the gene ID from locus_tag, gene name or protein id etc. Order in function defines order in which ID will be found"
	if 'gene' in feature.qualifiers and 'protein_id' in feature.qualifiers and 'locus_tag' in feature.qualifiers:
		return ('locus_tag|' + str(feature.qualifiers['locus_tag'][0]) + '|gene|' + str(feature.qualifiers['gene'][0]) + '|protein|' + str(feature.qualifiers['protein_id'][0]) + '|').replace(' ', '')
	if 'gene' in feature.qualifiers and 'protein_id' in feature.qualifiers:
		return ('gene|' + str(feature.qualifiers['gene'][0]) + '|protein|' + str(feature.qualifiers['protein_id'][0]) + '|').replace(' ', '')
	if 'locus_tag' in feature.qualifiers and 'protein_id' in feature.qualifiers:
		return ('gene|' + str(feature.qualifiers['locus_tag'][0]) + '|protein|' + str(feature.qualifiers['protein_id'][0]) + '|').replace(' ', '')
	if 'protein_id' in feature.qualifiers:
		return ('protein|' + feature.qualifiers['protein_id'][0] + '|').replace(' ', '')
	if 'gene' in feature.qualifiers and 'pacid' in feature.qualifiers:
		return ('gene|' + feature.qualifiers['gene'][0] + '|pacid|' + str(feature.qualifiers['pacid'][0]) + '|').replace(' ', '')
	if 'gene' in feature.qualifiers and 'ID' in feature.qualifiers:
		return ('gene|' + feature.qualifiers['gene'][0] + '|ID|' + str(feature.qualifiers['ID'][0]) + '|').replace(' ', '')
	if 'gene' in feature.qualifiers and 'translation' in feature.qualifiers:
		return ('gene|' + feature.qualifiers['gene'][0] + '|proteinmd5|' + hashlib.md5(str(feature.qualifiers['translation'][0]).encode('utf-8')).hexdigest() + '|').replace(' ', '')
	if 'pacid' in feature.qualifiers and 'translation' in feature.qualifiers:
		return ('pacid|' + feature.qualifiers['pacid'][0] + '|proteinmd5|' + hashlib.md5(str(feature.qualifiers['translation'][0]).encode('utf-8')).hexdigest() + '|').replace(' ', '')
	if 'locus_tag' in feature.qualifiers and 'translation' in feature.qualifiers:
		return ('locus_tag|' + feature.qualifiers['locus_tag'] + '|proteinmd5|' + hashlib.md5(str(feature.qualifiers['translation'][0]).encode('utf-8')).hexdigest() + '|').replace(' ', '')
	if 'gene' in feature.qualifiers:
		return ('gene|' + feature.qualifiers['gene'][0] + '|').replace(' ', '')
	if 'pacid' in feature.qualifiers:
		return ('pacid|' + feature.qualifiers['pacid'][0] + '|').replace(' ', '')
	if 'locus_tag' in feature.qualifiers:
		return ('locus_tag|' + feature.qualifiers['locus_tag'][0] + '|').replace(' ', '')
	return "no_tag_found"
	
def get_gene_id(feature):
	"Get the gene ID from locus_tag, gene name or protein id, in that order"
	if 'locus_tag' in feature.qualifiers:
		return feature.qualifiers['locus_tag'][0]
	if 'gene' in feature.qualifiers:
		return feature.qualifiers['gene'][0]
	if 'protein_id' in feature.qualifiers:
		return feature.qualifiers['protein_id'][0]
	return "no_gene_ID_found"


loci = []
anno = []


for record in SeqIO.parse(fileGBK, "genbank"):
	for feature in record.features:
		if feature.type == "CDS":
			gene_id = get_gene_id_CO(feature)
			with open(fileP450, 'r') as f:
				reader = csv.reader(f)
				for row in reader:
					for field in row:
						if field == gene_id:
							if gene_id not in loci:
								gene_ann = feature.qualifiers["translation"][0]
								anno.append(gene_ann)
								loci.append(gene_id)


dict_ann_id = dict(zip(anno, loci))

species = os.path.basename(os.path.normpath(os.path.abspath(os.path.join(fileGBK, os.pardir))))

with open(species + "_flanking_p450s_and_translations.fasta", "w+") as ofile:
	for i in range(len(anno)):
		ofile.write(">" +loci[i] + "\n" +anno[i] + "\n")

