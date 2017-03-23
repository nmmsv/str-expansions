import os, sys
import csv
import numpy as np
import matplotlib
import pickle
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



# if len(sys.argv) == 1:
# 	print '### Usage python vis_insert_size.py expName repoDir nCopyCount nCopy1 nCopy2 ... nCopyn extension'
# 	sys.exit()
# elif sys.argv[1] == "default":
# 	expName = 'HTT2'
# 	nCopyList = [21, 36]
# 	repoDir = '/storage/nmmsv/str-expansions/'
# 	ext = '.pdf'
# else:
# 	expName = sys.argv[1]
# 	repoDir = sys.argv[2]
# 	nCopyCount = int(sys.argv[3])
# 	nCopyList = []
# 	for i in range(4, 4 + nCopyCount):
# 		nCopyList.append(int(sys.argv[i]))
# 	ext = sys.argv[4 + nCopyCount]


expName = sys.argv[1]
repoDir = sys.argv[2]
ext = '.' + sys.argv[3]

expCaseDir = repoDir + 'experiments/' + expName

try:
	with open(expCaseDir + '/profile.txt', 'r') as f:
		arg_dict = pickle.load(f)
except:
	raise

refGenomeDir = arg_dict['ref_gen_dir']
nCopyList = arg_dict['num_copy']
bamFilter = arg_dict['bam_filter']
numThreads = arg_dict['num_threads']
numReadsList = arg_dict['num_reads']
locus = 'loci/' + arg_dict['locus'] + '.bed'
l1 = arg_dict['read_len']
flankLength = arg_dict['flank_len']
dist = arg_dict['read_ins_mean']
motif = arg_dict['motif']
refAllele = arg_dict['ref_allele_count']
heatLimit = arg_dict['heat_map_limit']
print '## Aligning Reads from Simulated Data: ##'
print '# Experiment Name: ', expName
print '# Repository Directory: ', repoDir
print '# List of nCopies for Generated Genomes: ', nCopyList
print '# Number of read pairs: ', numReadsList
print '## Visualization parameters: ##'
print '# Output extension: ', ext

print 

print '# Creating new directories..'
expDir = repoDir + 'experiments/'
alignDir = expDir + expName + '/aligned_read'
plotDir = expDir + expName + '/plots'
mkdir_p(plotDir)

heatMap = np.zeros((len(nCopyList), len(numReadsList)))
i = 0
for nc in nCopyList:
	nCopy = str(nc)
	j = 0
	inCaseDir =  alignDir + '/' + expName + '_nc' + nCopy
	for numReads in numReadsList:
		nRead = str(numReads)
		if len(numReadsList) > 1:
			subCaseName = 'nReads_' + nRead
			inSubCaseDir = inCaseDir + '/' + subCaseName
		else:
			inSubCaseDir = inCaseDir

		inPath = inSubCaseDir + '/' + expName + '_nc' + nCopy + '_nr' + nRead
		outPath = plotDir + '/' + expName + '_nc' + nCopy + '_nr' + nRead + ext

		os.system('samtools view ' + inPath + '_flt.bam | cut -f 9 | sort | uniq -c | \
			tr -s \' \' \'\\t\' | sed -e \'s/^[ \\t]*//\' > ' \
			+ expDir + expName + '/temp/tempHist.txt')


		Counts = []
		data = {}
		with open(expDir + expName + '/temp/tempHist.txt', 'r') as f:
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
		print 'nCopy = ', nCopy, '\tnRead=', nRead

		if np.abs(dist - param1[0]) < np.abs(dist - param0[0]):		# swap if 1 is ref allele
			temp = param0
			param0 = param1
			param1 = temp 		# Swap

		delta0 = np.abs(dist - param0[0]) / len(motif)
		delta1 = np.abs(dist - param1[0]) / len(motif)

		altAllele = np.abs(param1[0] - param0[0]) / len(motif) + refAllele
		err = (altAllele - nc) ** 2
		print 'peak0=', '{0:.5f}'.format(param0[0]), '\tDelta=', '{0:.4f}'.format(delta0), '*', motif
		print 'peak1=', '{0:.5f}'.format(param1[0]), '\tDelta=', '{0:.4f}'.format(delta1), '*', motif, '\tSqDiff=', '{0:.4f}'.format(err)
		print param0
		print param1
		if err > heatLimit:
			heatMap[i][j] = heatLimit
		else:
			heatMap[i][j] = err
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

		j = j + 1
	i = i + 1


fig2 = plt.figure()
ax = fig2.add_subplot(111)
cax = ax.imshow(heatMap, cmap = 'hot_r', interpolation = 'nearest')
ax.set_title('Squared Difference Error (#inferred_Alt - #true_Alt)^2\td=' + str(dist))
ax.set_xlabel('Average Coverage')
plt.xticks(np.arange(len(numReadsList)))
labels = [item.get_text() for item in ax.get_xticklabels()]
labels = [str(item * l1 / flankLength) for item in numReadsList]
ax.set_xticklabels(labels)

ax.set_ylabel('Number of Copies (Ref = 10)')
plt.yticks(np.arange(len(nCopyList)))
labels = [item.get_text() for item in ax.get_yticklabels()]
labels = [str(item) for item in nCopyList]
ax.set_yticklabels(labels)

fig2.colorbar(cax)
fig2.savefig(plotDir + '/' + expName + '_HeatMap.pdf')

