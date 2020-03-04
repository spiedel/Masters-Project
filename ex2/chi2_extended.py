from topfitter.analysis.frames import PredictionArray, AnalysisFrame
from topfitter.fitting import FitHandler
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare, chi2
from scipy import stats

def chi2_local(observed, expected, error):
	#function to calculate chi2 from two arrays
	if len(observed) != len(expected):
		print("Inputted data doesn't match length")
		return 0.

	chi2_sum = 0.
	for i in range(len(observed)):
	# calculate each chi2 part
		chi2_part = (observed[i]-expected[i])**2/error[i]
		chi2_sum = chi2_sum + chi2_part

	return chi2_sum


def getPred(fh, wcValue1, wcValue2):
	#funtion for predicting the theory values for a given pair of wilson coefficients
	#{'qq3_i33i', 'G', 'qq1_ii33', 'uu_i33i', 'qu1_33ii', 'qu8_33ii', 'qq1_i33i', 'qu8_ii33', 'qq3_ii33', 'uG_33', 'ud1_33ii', 'qd8_33ii', 'qu1_ii33', 'uu_ii33', 'qd1_33ii', 'ud8_33ii'}
	#returns
	prediction = fh.predict(qq3_i33i=wcValue1, qq1_i33i=wcValue2)#iqq3_i33i=wcValue2)#, uG_33=wcValue2)#
	return [x[0] for x in prediction.values]


wilco1 = "qq3\_i33i"
wilco2 = "qq1\_i33i"

#generate prediction array from data
#af = AnalysisFrame.from_hdf('/nfs/topfitter/sbrown/HDF/analyses/ATLAS_2017_I1604029/ATLAS_2017_I1604029.h5')
af = AnalysisFrame.from_hdf('../../HDF/CMS_2018_I1662081.h5')
pa = PredictionArray(af.xr)
fh = FitHandler(pa)

#get data values and error
obs = [x[0] for x in fh.reference.values]
err = [x[1] for x in fh.reference.values]
print(err)

print(pa.wilcos)
#print(pa.reference)
#print(fh.reference)
#print(pa.atpoint(G=1, qq3_i33i=1))
#print(chisquare(obs, [x[0] for x in pa.atpoint(G=3, qq3_i33i=3).values]))

#split by observable? use pa.reference.loc[level, 'value'] to get list instead

 
#plotting noValues*noValues with given x and y ranges
noValues = 40
xBins = np.linspace(-3, 3, noValues)
yBins = np.linspace(-3, 3, noValues)

i = 0
j = 0
chi2Array = np.zeros((noValues, noValues))
for x in xBins:
	for y in yBins:		
		pred = getPred(fh, x, y)
		chi2Value = chi2_local(obs, pred, err)
		chi2Array[i][j] = chi2Value
		j+=1
	j=0
	i+=1
#print(pred)

#print(chi2Array)
#Calculate delta chi2

chi2Array = chi2Array - chi2Array.min()

#calculate p
pArray = chi2.sf(chi2Array, 2)




#plot
fig, ax = plt.subplots()
#cont = ax.contourf(xBins, yBins, chi2Array)
cont = ax.contourf(xBins, yBins, pArray, [0.05]) 
print(cont.collections[0].get_paths()[0])
print(cont.allsegs[0])

#fig.colorbar(cont, ax=ax)
#plt.xlabel(wilco1)
#plt.ylabel(wilco2)
#plt.show()
#fig.savefig("contour_plots/contours_p_{}_{}.png".format(wilco1.replace("\\",""), wilco2.replace("\\","")))
