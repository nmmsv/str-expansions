# This script draw's a heat map for 
#		coverage vs. Delta (diff between ref and alt allele)
# To be called by $REPO/test/heat_map_wrapper.py

import os, sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
from scipy.stats import linregress

parser = argparse.ArgumentParser('Estimate STR length, calculate error, and draw insert size plot.')
parser.add_argument('--in-dir', 	type = str, required = True)
parser.add_argument('--plot-path', 	type = str, required = True)
parser.add_argument('--list0', 		type = str, required = True)
parser.add_argument('--list1', 		type = str, required = True)
parser.add_argument('--capt0',		type = str, required = True)
parser.add_argument('--capt1',		type = str, required = True)
parser.add_argument('--exp-name', 	type = str, required = True)


args = parser.parse_args()

in_dir = args.in_dir
plot_path = args.plot_path
list0 = [int(i) for i in args.list0.split()]
list1 = [int(i) for i in args.list1.split()]
capt0 = args.capt0
capt1 = args.capt1
exp_name = args.exp_name
limit = 1

heat_array = np.zeros((len(list0), len(list1)))

i = 0
determ_coeff_r2 = []
for it0 in list0:
	j = 0
	alt_allele_estm = []
	alt_allele_true = list1
	corr_plot_path = os.path.dirname(plot_path) + '/' + capt0 + str(it0) + '.pdf'
	for it1 in list1:
		in_path = in_dir + '/' + capt0 + str(it0) +\
						'_' + capt1 + str(it1) + '.txt'
		with open (in_path, 'r') as file_handle:
			for row in file_handle:
				if row[0] == '#' or row[0] == '0' or row.split() == []:
					pass
				else:
					cols = row.split()
					alt_allele_estm.append(float(cols[5]))
					error = float(cols[6])
					if error < limit:
						heat_array[i][j] = error
					else:
						heat_array[i][j] = limit
		j = j + 1

	slope, intercept, r_value, p_value, std_err = \
			linregress(alt_allele_true, alt_allele_estm)

	determ_coeff_r2.append(r_value ** 2)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.scatter(alt_allele_true, alt_allele_estm)
	ax.set_xlabel('True Alternate Allele')
	ax.set_ylabel('Estimated Alternate Allele')
	line_estm = [slope * v + intercept for v in alt_allele_true]
	ax.plot(alt_allele_true, line_estm, 'r')
	fig.savefig(corr_plot_path)
	i = i + 1


r2_plot_path = os.path.dirname(plot_path) + '/A_r2.pdf'
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(list0, determ_coeff_r2)
ax.set_title('Coefficient of determination vs. coverage')
ax.set_xlabel(capt0)
ax.set_ylabel('Coefficient of Determination (r^2)')
ax.set_xticks(list0)
fig.savefig(r2_plot_path)

fig2 = plt.figure()
ax = fig2.add_subplot(111)
cax = ax.imshow(heat_array, cmap = 'hot_r', interpolation = 'nearest')
ax.set_title('Squared Difference Error \
	(#inferred_Alt - #true_Alt)^2   d=')

ax.set_ylabel(capt0)
plt.yticks(np.arange(len(list0)))
labels = [item.get_text() for item in ax.get_yticklabels()]
labels = [str(item) for item in list0]
ax.set_yticklabels(labels)

ax.set_xlabel(capt1)
plt.xticks(np.arange(len(list1)))
labels = [item.get_text() for item in ax.get_xticklabels()]
labels = [str(item) for item in list1]
ax.set_xticklabels(labels)

fig2.colorbar(cax)
fig2.savefig(plot_path)

