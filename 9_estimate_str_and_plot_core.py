import os, sys
import csv
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import argparse

parser = argparse.ArgumentParser('Estimate STR length, calculate error, and draw insert size plot.')
parser.add_argument('--in-path', 	type = str, required = True)
parser.add_argument('--estm-path', 	type = str, required = True)
parser.add_argument('--plot-path', 	type = str, required = True)
parser.add_argument('--temp-dir', 	type = str, required = True)
parser.add_argument('--ref-allele', type = int, required = True)
parser.add_argument('--alt-allele', type = int, required = True)
parser.add_argument('--motif', 		type = str, required = True)
parser.add_argument('--dist-mean', 	type = int, required = True)


args = parser.parse_args()

in_path = args.in_path
estm_path = args.estm_path
plot_path = args.plot_path
temp_dir = args.temp_dir
ref_allele = args.ref_allele
alt_allele = args.alt_allele
motif = args.motif
dist_mean = args.dist_mean

print 'Plotting and STR estimation for' + in_path
estm_file_handle = open(estm_path, 'w')
# # StackOverflow
# from scipy.optimize import curve_fit
def gauss(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)


os.system('samtools view ' + in_path + ' | cut -f 9 | sort | uniq -c | \
	tr -s \' \' \'\\t\' | sed -e \'s/^[ \\t]*//\' > ' \
	+ temp_dir + 'tempHist.txt')

Counts = []
data = {}
with open(temp_dir + 'tempHist.txt', 'r') as f:
	for line in csv.reader(f, dialect='excel-tab'):
		x = np.abs(int(line[1]))
		if x > 0 and x < 10000:
			if x not in data:
				data[x] = int(line[0])
			else:
				data[x] = data[x] + int(line[0])

x = np.array(data.keys()).reshape(-1, 1)
y = np.array(data.values())

# http://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html#sklearn.mixture.GaussianMixture.fit
deq = np.array(x.repeat(y)).reshape(-1, 1)
model = GaussianMixture(n_components = 2,\
						max_iter = 100 )
model.fit(deq)

param0 = [model.means_[0 , 0], np.sqrt(model.covariances_[0, 0, 0]), model.weights_[0]]
param1 = [model.means_[1 , 0], np.sqrt(model.covariances_[1, 0, 0]), model.weights_[1]]

if np.abs(dist_mean - param1[0]) < np.abs(dist_mean - param0[0]):		# swap if 1 is ref allele
	temp = param0
	param0 = param1
	param1 = temp 		# Swap

delta0 = np.abs(dist_mean - param0[0]) / len(motif)
delta1 = np.abs(dist_mean - param1[0]) / len(motif)

est_allele = np.abs(param1[0] - param0[0]) / len(motif) + ref_allele

error = ((est_allele - alt_allele) / alt_allele) ** 2


estm_file_handle.write('#\t\tReference Allele:\t' + str(ref_allele) + '\n')
estm_file_handle.write('#(True)\t\tAlternate Allele:\t' + str(alt_allele) + '\n')
estm_file_handle.write('#(Estimate)\t\tAlternate Allele:\t' + str(est_allele) + '\n' + '\n')
estm_file_handle.write('#peak\tmean\tstd_dev\tWeight\tDelta(*motif)\tSqDiff\n')
estm_file_handle.write('0\t{0:.2f}'.format(param0[0]) + \
		'\t{0:.2f}'.format(param0[1]) + '\t{0:.2f}'.format(param0[2]) + \
		'\t{0:.2f}'.format(delta0) +'\t0.0\n')
estm_file_handle.write('1\t{0:.2f}'.format(param1[0]) + \
		'\t{0:.2f}'.format(param1[1]) + '\t{0:.2f}'.format(param1[2]) + \
		'\t{0:.2f}'.format(delta1) + '\t{0:.3f}'.format(error) + '\n')

# the histogram of the data'
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(x, y)
ax.set_xlabel('Mapped Mate Pair Outer Distance')
ax.set_ylabel('Count')


# ax.plot(x,bimodal(x,*params),color='red',lw=3,label='model')
pred  = model.predict(x)
x0 = x[pred == 0]
x1 = x[pred == 1]

cent = (param1[0] + param0[0]) / 2
stdd = max(param0[1], param1[1])
xplot = np.linspace(cent - 10 * stdd, cent + 10 * stdd, 100)
ax.plot(xplot, max(y) * gauss(xplot, param0[0], param0[1], param0[2]), \
	color = 'red', lw = 3, label = 'model0')
ax.plot(xplot, max(y) * gauss(xplot, param1[0], param1[1], param1[2]), \
	color = 'blue', lw = 3, label = 'model1')
ax.axvline(x = dist_mean, color = 'k', linestyle = '--')
fig.savefig(plot_path)


estm_file_handle.close()

