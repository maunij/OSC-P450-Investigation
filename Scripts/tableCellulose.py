# opens {spe}_flanking_total_counts.csv_pvalue and grabs the pvalue for Cellulose Synthase

import os
import csv
import sys

fileCSV = sys.argv[1]

x = fileCSV.split("_flanking")
i = False

with open(fileCSV, "r") as f:
	with open("tableCellulose.csv", "a+") as ofile:
		reader = csv.reader(f)
		writer = csv.writer(ofile)
		for row in reader:
			if row[0] == "Cellulose synthase":
				i = True
				pval = row[1]
				writer.writerow([x[0], pval])
				break
		if i == False:
			pval = "n/a"
			writer.writerow([x[0], pval])