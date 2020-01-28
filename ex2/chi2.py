from topfitter.analysis.frames import PredictionArray, AnalysisFrame
from topfitter.fitting import FitHandler
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare, chi2
from scipy import stats

def chi2_local(observed, expected):
	#function to calculate chi2 from two arrays
	if len(observed) != len(expected):
		print("Inputted data doesn't match length")
		return 0.

	chi2_sum = 0.
	for i in range(len(observed)):
	# calculate each chi2 part
		chi2_part = (observed[i]-expected[i])**2/expected[i]
		chi2_sum = chi2_sum + chi2_part

	return chi2_sum


def getPred(pa, wcValue1, wcValue2):
	#funtion for predicting the theory values for a given pair of wilson coefficients
	#returns
	prediction = pa.atpoint(G=wcValue1, uG_33=wcValue2)#, qq3_i33i=wcValue2)
	return [x[0] for x in prediction.values]


#generate prediction array from data
#af = AnalysisFrame.from_hdf('/nfs/topfitter/sbrown/HDF/analyses/ATLAS_2017_I1604029/ATLAS_2017_I1604029.h5')
af = AnalysisFrame.from_hdf('../../HDF/CMS_2018_I1662081.h5')
pa = PredictionArray(af.xr)

#get data values and error
obs = [x[0] for x in pa.reference.values]
err = [x[1] for x in pa.reference.values]

print(pa.wilcos)

print(pa.reference)
#print(pa.atpoint(G=1, qq3_i33i=1).values)
#print(chisquare(obs, [x[0] for x in pa.atpoint(G=3, qq3_i33i=3).values]))

#split by observable? use pa.reference.loc[level, 'value'] to get list instead

 
#plotting noValues*noValues with given x and y ranges
noValues = 30
xBins = np.linspace(-10, 5, noValues)
yBins = np.linspace(-10, 5, noValues)

i = 0
j = 0
chi2Array = np.zeros((noValues, noValues))
for x in xBins:
	for y in yBins:		
		pred = getPred(pa, x, y)
		chi2Value = chisquare(obs, f_exp=pred)[0]
		chi2Array[i][j] = chi2Value
		j+=1
	j=0
	i+=1

print(chi2Array)
#Calculate delta chi2

chi2Array = chi2Array - chi2Array.min()

#calculate p from the survival function (1-cdf)
pArray = chi2.sf(chi2Array, 2)

print(chi2Array)
print(pArray)

#plot
fig, ax = plt.subplots()

#cont = ax.contourf(xBins, yBins, chi2Array)
cont = ax.contourf(xBins, yBins, pArray, [0.65, 0.95, 0.99, 1.0]) 
fig.colorbar(cont, ax=ax)

plt.show()
fig.savefig("contour_plots/contours_p_{}.png".format(noValues))


