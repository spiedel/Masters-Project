from topfitter.analysis import AnalysisFrame, PredictionArray
from topfitter.fitting import FitHandler
import matplotlib.pyplot as plt
import numpy as np

#plot stepped bar chart from bin edges and weights
def plotHist(ax, bins, values, label="", colour='b', linestyle='-'):
	ax.plot([x[0] for x in bins] + [bins[-1][1]],
		np.append(values, values[-1]),
		ds='steps-post', label=label, ls=linestyle, color=colour)

#form prediction array
af = AnalysisFrame.from_hdf('../../HDF/analyses/ATLAS_1407.0371/ATLAS_1407.0371.h5')
pa = PredictionArray(af)

#setup fit handler
fh = FitHandler(pa)

print(fh.wilcos)

#predictions
sm = fh.predict()
wc = fh.predict(qq3_i33i=1)
wcName = "qq3\_i33i"

#iterate through different observables, plotting each
for entry in sm.index.levels[0][:-1]:
	
	fig=plt.figure()
	ax = plt.subplot(1,1,1)

	print(entry)
	
	#define data to be plotted for each observable
	smData = sm.loc[entry]
	smBins = smData.index
	wcBins = wc.loc[entry].index
	
	smValues = sm.loc[(entry, 'value')]
	smErr = sm.loc[(entry, 'error')]
	wcValues = wc.loc[(entry,'value')]
	wcErr = wc.loc[(entry, 'error')]

	print(smBins)
	print(smValues)

	#plot standard model against one wc set to one
	plotHist(ax, smBins, smValues.values, "SM\_nlo", "b")
	plotHist(ax, wcBins, wcValues.values, wcName, "r")

	#create legend
	ax.legend()

	#save plot under entry name
	name = entry.replace('/', '_')
	print(name)
	plt.show()	
	fig.savefig(name, format='pdf')
