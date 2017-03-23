import os, sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn.mixture import GaussianMixture

import errno
# StackOverflow
def mkdir_p(path):
    try:
        os.makedirs(path)
        print 'Created path: ', path
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

# # StackOverflow
# from scipy.optimize import curve_fit
def gauss(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)



if len(sys.argv) == 1:
	print '### Usage python vis_insert_size.py expName repoDir nCopyCount nCopy1 nCopy2 ... nCopyn extension'
	sys.exit()
elif sys.argv[1] == "default":
	expName = 'HTT2'
	nCopyList = [21, 36]
	repoDir = '/storage/nmmsv/str-expansions/'
	ext = '.pdf'
else:
	expName = sys.argv[1]
	repoDir = sys.argv[2]
	nCopyCount = int(sys.argv[3])
	nCopyList = []
	for i in range(4, 4 + nCopyCount):
		nCopyList.append(int(sys.argv[i]))
	ext = sys.argv[4 + nCopyCount]

print '## Aligning Reads from Simulated Data: ##'
print '# Experiment Name: ', expName
print '# Repository Directory: ', repoDir
print '# List of nCopies for Generated Genomes: ', nCopyList
print '## Visualization parameters: ##'
print '# Output extension: ', ext

print 


mkdir_p(repoDir + expName + '/plots')
for nc in nCopyList:
	nCopy = str(nc)

	outPath = repoDir + expName + '/plots/' + expName + '_' + nCopy + ext

	os.system('samtools view ' + repoDir + expName +\
		'/aligned_read/' + expName + '_' + nCopy +\
		'/' + expName + '_' + nCopy + '_flt.bam | cut -f 9 | sort | uniq -c | \
		tr -s \' \' \'\\t\' | sed -e \'s/^[ \\t]*//\' > ' \
		+ repoDir + expName + '/temp/tempHist.txt')


	Counts = []
	data = {}
	with open(repoDir + expName + '/temp/tempHist.txt', 'r') as f:
		for line in csv.reader(f, dialect='excel-tab'):
			x = np.abs(int(line[1]))
			if x > 0 and x < 10000:
				if x not in data:
					data[x] = int(line[0])
				else:
					data[x] = data[x] + int(line[0])

	x = np.array(data.keys()).reshape(-1, 1)
	y = np.array(data.values())
	# StackOverflow
	# params,cov=curve_fit(bimodal,x,y, p0=[1500, 100, 100, 1500, 100, 100], maxfev=30000)
	# # sigma=np.sqrt(np.diag(cov))
	# print '1) mean: ', params[0], '\tstd_dev: ', params[1], '\tA: ', params[2]
	# print '2) mean: ', params[3], '\tstd_dev: ', params[4], '\tA: ', params[5]
	# print
	##

	# http://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html#sklearn.mixture.GaussianMixture.fit
	deq = np.array(x.repeat(y)).reshape(-1, 1)
	model = GaussianMixture(n_components = 2,\
							max_iter = 100 )

	model.fit(deq)

	param0 = [model.means_[0 , 0], np.sqrt(model.covariances_[0, 0, 0]), model.weights_[0]]
	param1 = [model.means_[1 , 0], np.sqrt(model.covariances_[1, 0, 0]), model.weights_[1]]
	print '######'
	print 'nCopy = ', nCopy
	print 'peak1 = ', param0[0], '\tDelta = ', (dist - param0[0]) / len(motif)
	print param0
	print param1
	print model.converged_

	# the histogram of the data'
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.scatter(x, y)
	ax.set_xlabel('Insert Size')
	#ax.set_xlim(left=0, right=100)
	ax.set_ylabel('Count')
	#ax.set_ylim(bottom=0, top=3000)


	# ax.plot(x,bimodal(x,*params),color='red',lw=3,label='model')
	pred  = model.predict(x)
	x0 = x[pred == 0]
	x1 = x[pred == 1]

	cent = (param1[0] + param0[0]) / 2
	stdd = max(param0[1], param1[1])
	xplot = np.linspace(cent - 10 * stdd, cent + 10 * stdd, 100)
	ax.plot(xplot, max(y) * gauss(xplot, param0[0], param0[1], param0[2]), color = 'red', lw = 3, label = 'model0')
	ax.plot(xplot, max(y) * gauss(xplot, param1[0], param1[1], param1[2]), color = 'blue', lw = 3, label = 'model1')
	fig.savefig(outPath)




