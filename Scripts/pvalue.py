import scipy
from scipy import stats
import csv
import sys

fileA = sys.argv[1]

def main(genome):
	with open(fileA) as fi:
		reader = csv.reader(fi)
		labels = reader.next()[1:]
		flank_counts = reader.next()[1:]
		all_counts = reader.next()[1:]
	
	label_to_p = {}
	
	for lab in labels:
		index = labels.index(lab)
		lab_flanking = int(flank_counts[index])
		lab_remaining = int(all_counts[index]) - lab_flanking
		if lab_remaining < 0:
			lab_remaining = 0
		# TODO: bug with CYP_classify script not counting duplicates properly
		other_flanking = sum(map(int, [x for i,x in enumerate(flank_counts) if i!=index]))
		other_remaining = sum(map(int, [x for i,x in enumerate(all_counts) if i!=index]))
		
		LOD, pvalue = scipy.stats.fisher_exact([[lab_flanking,lab_remaining],[other_flanking,other_remaining]], alternative='greater')
		label_to_p[lab] = pvalue
		print(str(pvalue))
		print(str(lab_flanking))
	
	with open(str(fileA) + '_pvalue.csv', "a+") as pf:
		writer = csv.writer(pf)
		for lab in labels:
			writer.writerow([str(lab), str(label_to_p[lab])])
	
if __name__ == "__main__":
    main(*sys.argv[1:])