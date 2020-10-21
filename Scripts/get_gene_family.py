#!/usr/bin/env python
"""
Usage:
    get_gene_family.py <genome> <hmm>
"""

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
import logging
import os
import errno
import pprint
import subprocess
from collections import defaultdict
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

import warnings
# Don't display the SearchIO experimental warning, we know this.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	from Bio import SearchIO

import pprint
from operator import itemgetter
import operator
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, ExactPosition
import csv
import re
import hashlib
import sys

import numpy as np
import scipy
# from scipy import stats

# if sys.version_info[:2] < (3, 3):
#     old_print = print
#     def print(*args, **kwargs):
#         flush = kwargs.pop('flush', False)
#         old_print(*args, **kwargs)
#         file = kwargs.get('file', sys.stdout)
#         if flush and file is not None:
#             file.flush()



class Signature(object):
	"""Secondary metabolite signature"""
	def __init__(self, name, _type, description, cutoff, path):
		self.name = name
		self.type = _type
		self.description = description
		self.cutoff = cutoff
		self.path = path

class HmmSignature(Signature):
	"""HMM signature"""
	def __init__(self, name, description, cutoff, hmm_file):
		self.hmm_file = (hmm_file)
		self.name = name
		super(HmmSignature, self).__init__(name, 'model',
			  description, cutoff, hmm_file)


profile_hmm = sys.argv[2]

_signature_profiles = [HmmSignature(profile_hmm, profile_hmm, -1, profile_hmm)]


def detect_signature_genes_withnofilter_CO(seq_record, enabled_clustertypes):
	"Function to be executed by module"
	logging.info('Detecting gene clusters using HMM library')
	feature_by_id = get_feature_dict(seq_record)
	full_fasta = get_multifasta_CO(seq_record)
	results = []
	sig_by_name = {}
	results_by_id = {}
	for sig in _signature_profiles:
		runresults = run_hmmsearch_CO(sig.path, full_fasta)
		for runresult in runresults:
			#Store result if it is above cut-off
			for hsp in runresult.hsps:
				if hsp.bitscore > sig.cutoff:
					results.append(hsp)
					if not results_by_id.has_key(hsp.hit_id):
						results_by_id[hsp.hit_id] = [hsp]
					else:
						results_by_id[hsp.hit_id].append(hsp)
		sig_by_name[sig.name] = sig
	results_nofilter, results_by_id_nofilter = results[:], results_by_id.copy()
	# Get overlap tables (for overlap filtering etc)
	overlaps = get_overlaps_table(seq_record)
	# Filter results of overlapping genes
	results, results_by_id = filter_result_overlapping_genes_CO_nohmm(results, results_by_id, overlaps, feature_by_id)
# 	# Filter multiple results of the same model in one gene
	results, results_by_id = filter_result_multiple(results, results_by_id)
	return results, results_by_id, results_nofilter, results_by_id_nofilter


def detect_signature_genes_nofilter_CO(seq_record, enabled_clustertypes):
	"Function to be executed by module"
	logging.info('Detecting gene clusters using HMM library')
	feature_by_id = get_feature_dict(seq_record)
	full_fasta = get_multifasta_CO(seq_record)
	results = []
	sig_by_name = {}
	results_by_id = {}
	for sig in _signature_profiles:
		runresults = run_hmmsearch_CO(sig.path, full_fasta)
		for runresult in runresults:
			#Store result if it is above cut-off
			for hsp in runresult.hsps:
				if hsp.bitscore > sig.cutoff:
					results.append(hsp)
					if not results_by_id.has_key(hsp.hit_id):
						results_by_id[hsp.hit_id] = [hsp]
					else:
						results_by_id[hsp.hit_id].append(hsp)
		sig_by_name[sig.name] = sig
	return results, results_by_id
	

def get_feature_dict(seq_record):
	"""Get a dictionary mapping features to their IDs"""
	features = get_cds_features(seq_record)
	feature_by_id = {}
	for feature in features:
		gene_id = get_gene_id_CO(feature)
		feature_by_id[gene_id] = feature
	return feature_by_id

def get_cds_features(seq_record):
	"Return all CDS features for a seq_record"
	return get_all_features_of_type(seq_record, "CDS")

def get_all_features_of_type(seq_record, types):
	"Return all features of the specified types for a seq_record"
	if isinstance(types, str):
		# force into a tuple
		types = (types, )
	features = []
	for f in seq_record.features:
		if f.type in types:
			features.append(f)
	return features

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


def get_gene_annotation(feature):
	"Get the gene annotation from the product qualifier"
	if 'product' in feature.qualifiers:
		return feature.qualifiers['product'][0]
	return "unannotated gene"


def run_hmmsearch_CO(query_hmmfile, target_sequence):
	"Run hmmsearch"
	command = ["hmmsearch",
			   query_hmmfile, '-']
	try:
		out, err, retcode = execute(command, input=target_sequence)
	except OSError:
		return []
	if retcode != 0:
		logging.debug('hmmsearch returned %d: %r while searching %r', retcode,
						err, query_hmmfile)
		return []
	res_stream = StringIO(out)
	results = list(SearchIO.parse(res_stream, 'hmmer3-text'))
	return results


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

def get_overlaps_table(seq_record):
	"""Given seq_record, returns array of overlapping genes
	and the corresponding gene_id-to-indexes"""
	overlaps = []
	overlap_by_id = {}
	features = get_cds_features(seq_record)
	if len(features) < 1:
		return overlaps, overlap_by_id
	features.sort(key=lambda feature: feature.location.start)
	i = 0
	j = i + 1
	cds_queue = []
	while i < len(features):
		if j >= len(features):
			break
		cds = features[i]
		ncds = features[j]
		if (cds.location.end <= ncds.location.start + 1):
			overlaps.append([])
			cds_queue.append(cds)
			for cds in cds_queue:
				overlap_by_id[get_gene_id_CO(cds)] = len(overlaps) - 1
				overlaps[-1].append(cds)
			cds_queue = []
			i = j
		else:
			if (cds.location.end < ncds.location.end):
				cds_queue.append(cds)
				i = j
			else:
				cds_queue.append(ncds)
		j += 1
	overlaps.append([])
	cds_queue.append(features[i])
	for cds in cds_queue:
		overlap_by_id[get_gene_id_CO(cds)] = len(overlaps) - 1
		overlaps[-1].append(cds)
	return overlaps, overlap_by_id

def filter_result_overlapping_genes_CO_nohmm(results, results_by_id, overlaps, feature_by_id):
	# filter results of overlapping genes (only gene with the best score can retain its result, unless it is a predicted gene from selenoprofiles)
	overlap_id_with_result = {}
	for cds in results_by_id.keys():
		if overlaps[1][cds] not in overlap_id_with_result.keys():
			overlap_id_with_result[overlaps[1][cds]] = [cds]
		elif cds not in overlap_id_with_result[overlaps[1][cds]]:
			overlap_id_with_result[overlaps[1][cds]].append(cds)
#	 print overlap_id_with_result
	for overlap_id in overlap_id_with_result.keys():
		best_hit_scores = {}
		equal_hits_by_cds = {}
		for cds in overlap_id_with_result[overlap_id]:
			for hit in results_by_id[cds]:
				feature = feature_by_id[hit.hit_id]
				if 'selenoprofiles' in hit.hit_id:
					hit.bitscore = 0
				if (hit.query_id not in best_hit_scores) or (best_hit_scores[hit.query_id] < hit.bitscore):
					best_hit_scores[hit.query_id] = hit.bitscore	
		for cds in overlap_id_with_result[overlap_id]:
			to_delete = []
			for hit in results_by_id[cds]:
				feature = feature_by_id[hit.hit_id]
				if (hit.bitscore < best_hit_scores[hit.query_id]):
					to_delete.append(hit)
			for hit in to_delete:
				del results[results.index(hit)]
				del results_by_id[cds][results_by_id[cds].index(hit)]
				if len(results_by_id[cds]) < 1:
					del results_by_id[cds]
		# If range is identical, choose the hit with the earliest alphabetical annotation
		for cds in overlap_id_with_result[overlap_id]:
			to_delete = []
			if cds in results_by_id.keys():
				for hit in results_by_id[cds]:
					feature = feature_by_id[hit.hit_id]
					if (hit.bitscore == best_hit_scores[hit.query_id]):
						equal_hits_by_cds[cds] = hit
		sorted_equal_hits_by_cds = sorted(equal_hits_by_cds.items(), key=operator.itemgetter(1))
		del equal_hits_by_cds[sorted_equal_hits_by_cds[0][0]]
		for cds in equal_hits_by_cds.keys():
			for hit in results_by_id[cds]:
				del results[results.index(hit)]
				del results_by_id[cds][results_by_id[cds].index(hit)]
				if len(results_by_id[cds]) < 1:
					del results_by_id[cds]
	return results, results_by_id

def filter_result_multiple(results, results_by_id):
	#Filter multiple results of the same model within a gene
	for cds in results_by_id.keys():
		best_hit_scores = {}
		to_delete = []
		for hit in results_by_id[cds]:
			if (hit.query_id not in best_hit_scores) or (best_hit_scores[hit.query_id] < hit.bitscore):
				best_hit_scores[hit.query_id] = hit.bitscore
		for hit in results_by_id[cds]:
			if (hit.bitscore < best_hit_scores[hit.query_id]):
				to_delete.append(hit)
		for hit in to_delete:
			del results[results.index(hit)]
			del results_by_id[cds][results_by_id[cds].index(hit)]
			if len(results_by_id[cds]) < 1:
				del results_by_id[cds]
	return results, results_by_id



def get_multifasta_CO(seq_record):
	"""Extract multi-protein FASTA from all CDS features in sequence record"""
	features = get_cds_features(seq_record)
	all_fastas = []
	seen_gene_id_counter = {}
	for feature in features:
		gene_id = get_gene_id_CO(feature)
		description = get_gene_annotation(feature)
		# Check to see if there is a translation entry
		if 'translation' not in feature.qualifiers:
			logging.debug("No translation entry for %s, skipping.." %(gene_id))
			continue
		fasta_seq = feature.qualifiers['translation'][0]
		if "-" in str(fasta_seq):
			fasta_seq = Seq(str(fasta_seq).replace("-",""), generic_protein)
		# Never write empty fasta entries
		if len(fasta_seq) == 0:
			logging.debug("No translation sequence for %s, skipping" %(gene_id))
			continue
		gene_id = str(gene_id) + ' ' + str(description)
		all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
	full_fasta = "\n".join(all_fastas)
	
	return full_fasta


def check_duplicate_gene_ids_CO(seq_record):
    "Fix duplicate locus tags so that they are different"
    NO_TAG = "no_tag_found"
    high_water_mark = 0
    all_ids = defaultdict(lambda: False)
    seq_ids = get_cds_features(seq_record)
    for cdsfeature in seq_ids:
    	gene_id = get_gene_id_CO(cdsfeature)
#     	cdsfeature.qualifiers['product'][0].replace(' ', '')
#     	cdsfeature.qualifiers['locus_tag'][0].replace(' ', '')
#     	if 'gene' in cdsfeature.qualifiers:
#     		cdsfeature.qualifiers['gene'][0].replace(' ', '')
        if not all_ids[gene_id]:
            all_ids[gene_id] = True
        else:
            if gene_id == NO_TAG:
                x = high_water_mark + 1
            else:
                x = 1
            id_str = "%s_%s" % ( gene_id[:8], x)
            while all_ids[id_str]:
                x += 1
                id_str = "%s_%s" % ( gene_id[:8], x)
            print("generated id %r,", id_str)
            cdsfeature.qualifiers['product'] = [id_str]
            cdsfeature.qualifiers['locus_tag'] = [id_str]
            all_ids[id_str] = True
            if gene_id == NO_TAG:
                high_water_mark = x


def main(genome, profile_hmm, output):
	species = os.path.basename(os.path.normpath(os.path.abspath(os.path.join(genome, os.pardir))))
	print('Looking for %s families in %s ...' %(profile_hmm, species))
	feature_hits = []
	gene_id_hits = []
	for seq_record in SeqIO.parse(genome, "gb"):
		check_duplicate_gene_ids_CO(seq_record)
		results, results_by_id, results_nofilter, results_by_id_nofilter = detect_signature_genes_withnofilter_CO(seq_record, 'plant')
		for result in results:
			gene_id_hits.append(result.hit.id)
		features = get_cds_features(seq_record)
		for feature in features:
			gene_id = get_gene_id_CO(feature)
			if gene_id in gene_id_hits:
				record = SeqRecord(Seq(str(feature.qualifiers['translation'][0])), id=gene_id)
				feature_hits.append(record)
	SeqIO.write(feature_hits, output, "fasta")

if __name__ == "__main__":
    main(*sys.argv[1:])

