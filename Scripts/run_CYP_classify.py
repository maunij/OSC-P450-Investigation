#!/usr/bin/env python

from __future__ import division
import random
import subprocess
import re
import csv
import logging
from Bio import SeqIO
import sys

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


def execute(commands, input=None):
	"Execute commands in a system-independent manner"
	if input is not None:
		stdin_redir = subprocess.PIPE
	else:
		stdin_redir = None
	try:
		proc = subprocess.Popen(" ".join(commands), stdin=stdin_redir, shell=True,
								stdout=subprocess.PIPE,
								stderr=subprocess.PIPE)
		out, err = proc.communicate(input=input)
		retcode = proc.returncode
		return out, err, retcode
	except OSError as e:
		logging.debug("%r %r returned %r", commands, input[:40] if input is not None else None, e)
		raise

def parse_cdhit_file(file_path):
	"Parsing CD-HIT result file into array"
	clusters = []
	gene_to_cluster = {}
	
	cdhit_file = open(file_path, "r")
	re_content = re.compile(".+\>(?P<gene_id>.+)(\.\.\.) (\*|at (?P<identity_val>[.*0-9]+)).*")
	parsing_content = False
	for line in cdhit_file.readlines():
		if line[0] == '>':
			new_cluster = True
		else:
			match = re_content.match(line)
			if match and (new_cluster or len(clusters) > 0):
				gene_id = match.group("gene_id")
				identity_val = -1.00
				if new_cluster:
					clusters.append({"rep_gene" : "", "genes" : {}})
					new_cluster = False
				if match.group("identity_val") == None:
					clusters[len(clusters) - 1]["rep_gene"] = gene_id
					identity_val = 100.00
				else:
					identity_val = float(match.group("identity_val"))
				clusters[len(clusters) - 1]["genes"][gene_id] = identity_val
				gene_to_cluster[gene_id] = len(clusters) - 1
	cdhit_file.close()
	
	return clusters, gene_to_cluster

def keywithmaxval(d):
     """ a) create a list of the dict's keys and values; 
         b) return the key with the max value"""  
     v=list(d.values())
     k=list(d.keys())
     return k[v.index(max(v))]



def main (argv):
	def printhelp():
		print("\nusage:\nCYP_classify.py CYPs.fasta full_CYP_database.clstr CYP_database_representatives.fasta outdir")
	
	if len(sys.argv)<5:
		print('\nNot enough arguments.')
		printhelp()
		sys.exit()


	################################################################## PARSE INPUT
	CYPs_fasta	=sys.argv[1]
	CYP_database_clstr =sys.argv[2]
	CYP_database_repres =sys.argv[3]
	outdir	=sys.argv[4]
	
	
	print('Reading CYP database...')
	# Read database
	training_clusters, training_gene_to_cluster = parse_cdhit_file(CYP_database_clstr)
	cluster_to_fam = {}
	training_genes_set = set([])
	cluster_count = -1
	for cluster in training_clusters:
		cluster_count +=1
		fam_set = set([])
		fam_list = []
		total_count = 0
		for gene_id in cluster['genes'].keys():
			training_genes_set.update([gene_id])
			total_count +=1
			fam = gene_id.split("_sid")[0]
			fam_set.update([fam])
			fam_list.append(fam)
		
		fam_counts = {}
		for fam in fam_set:
			fam_count = fam_list.count(fam)
			fam_counts[fam] = fam_count
		
		cluster_to_fam[cluster_count] = fam_counts
	
	database_records = []
	for seq_record in SeqIO.parse(CYP_database_repres, 'fasta'):
		database_records.append(seq_record)
	print('Reading input CYPs...')
	
	test_records = []
	for seq_record in SeqIO.parse(CYPs_fasta, 'fasta'):
		test_records.append(seq_record)
		if seq_record.name[:19] in training_genes_set:
			print('Warning, the name of sequence %s clashes with a sequence in the database! Please check!' %(seq_record.name))
	
	print('Combining CYPs with database...')
	all_records = database_records + test_records
	SeqIO.write(all_records, str(outdir) + '/' + str(CYPs_fasta) + '_combined_sequences.fasta', 'fasta')
	
	print('Clustering database...')
	command = ["cd-hit", "-i", str(outdir) + '/' + str(CYPs_fasta) + '_combined_sequences.fasta', "-o", str(outdir) + '/' + str(CYPs_fasta) + '_combined_sequences_CDHIT_OUT', "-c", "0.50", "-n", "3", "-T", "0"]
	error = False
	retcode = 0
	out, err, retcode = execute(command)
	
	print('Reading CD-HIT output...')
	test_clusters, test_gene_to_cluster = parse_cdhit_file(str(outdir) + '/' + str(CYPs_fasta) + '_combined_sequences_CDHIT_OUT.clstr')
	
	print('Classifying CYPs...')
	gene_to_classification = {}
	cluster_count = -1
	for cluster in test_clusters:
		fam_scores = []
		cluster_count += 1
		for gene_id in cluster['genes'].keys():
			if gene_id in training_genes_set:
				cluster_id = training_gene_to_cluster[gene_id]
				fam_scores.append(cluster_to_fam[cluster_id])
		if len(fam_scores) > 1:
			new_fam_scores = {}
			for score_count in fam_scores:
				for fam in score_count.keys():
					if fam not in new_fam_scores.keys():
						new_fam_scores[fam] = score_count[fam]
					else:
						new_fam_scores[fam] = new_fam_scores[fam] + score_count[fam]
			fam_scores = [new_fam_scores]
		for gene_id in cluster['genes'].keys():
			if gene_id not in training_genes_set:
				gene_to_classification[gene_id] = fam_scores
	
	print('Writing CYP classifications...')
	with open(str(outdir)+'/'+str(CYPs_fasta)+'_classified_CYPs_summary.csv', 'w') as fo:
		writer = csv.writer(fo)
		for gene_id, fam_scores in gene_to_classification.iteritems():
			if len(fam_scores) == 0:
				predicted_fam = 'Uncertain'
			elif len(fam_scores[0]) == 1:
				predicted_fam = fam_scores[0].keys()[0]
			else:
				total_count = sum(fam_scores[0].itervalues())
				for fam, score_count in fam_scores[0].iteritems():
					fam_scores[0][fam] = score_count/float(total_count)
	 			# If all of the ratios in the dictionary are equal, pick a random score
				if len(set(fam_scores[0].itervalues())) <= 1:
					predicted_fam = [key for key in fam_scores[0].keys()]
				else:
					predicted_fam = keywithmaxval(fam_scores[0])
			writer.writerow([gene_id, predicted_fam])


if __name__ == "__main__":
	main(sys.argv[1:])