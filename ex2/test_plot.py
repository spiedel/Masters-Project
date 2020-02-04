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
af = AnalysisFrame.from_hdf('../../HDF/CMS_2018_I1662081.h5')
pa = PredictionArray(af)

#setup fit handler
#fh = FitHandler(pa)

print(pa.wilcos)

#predictions
data = pa.reference
wc = pa.atpoint(ud1_33ii=1, qu8_ii33=1, qq1_ii33=1, qd8_33ii=1, ud8_33ii=1, uG_33=1, uu_ii33=1, uu_i33i=1, qq3_ii33=1, qu8_33ii=1, qq1_i33i=1, qu1_33ii=1, qu1_ii33=1, qd1_33ii=1, qq3_i33i=1, G=1)
#wc = fh.predict(G=1)
wcName = "All Wilson Coefficients set to 1"

print(pa.reference.columns)
obs = set([x[0] for x in pa.reference.index])

#obs = pa.reference.index
print(obs)

#iterate through different observables, plotting each
for entry in obs:
	print(entry)
	
	fig=plt.figure()
	ax = plt.subplot(1,1,1)

	print(entry)
	
	#define data to be plotted for each observable
	measuredData = data.loc[entry]
	measuredBins = measuredData.index
	wcBins = wc.loc[entry].index
	
	measuredValues = data.loc[(entry, 'value')]
	measuredErr = data.loc[(entry, 'error')]
	wcValues = wc.loc[(entry,'value')]
	wcErr = wc.loc[(entry, 'error')]

	print(measuredBins)
	print(measuredValues)

	#plot standard model against one wc set to one
	plotHist(ax, measuredBins, measuredValues.values, "data", "b")
	plotHist(ax, measuredBins, measuredValues.values + measuredErr.values, "error up", "b", "--")
	plotHist(ax, measuredBins, measuredValues.values - measuredErr.values, "error down", "b", "--")
	plotHist(ax, wcBins, wcValues.values, wcName, "r")

	#create legend
	ax.legend()

	#save plot under entry name
	name = "comp_plots/" + entry.replace('/', '_') + ".pdf"
	print(name)
	plt.show()	
	fig.savefig(name, format='pdf')
