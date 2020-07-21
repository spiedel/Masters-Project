from topfitter.analysis.frames import PredictionArray, AnalysisFrame
from topfitter.fitting import FitHandler
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare, chi2
from scipy import stats

def chi2_local(observed, expected, error):
	#function to calculate chi2 from two arrays
	if observed.size != expected.size:
		print("Inputted data doesn't match length")
		return 0.
	chi2_sum = 0.
	for i in range(observed.size):
	# calculate each chi2 part
		if error[i] == 0.:
			print("Check this ", observed[i], " at point", i)
		else:
			diff = observed[i] - expected[i]
			chi2_part = (diff*diff)/(error[i]*error[i])
			chi2_sum = chi2_sum + chi2_part
	return chi2_sum



def getPred(pa, wcdict_tmp):
	#funtion for predicting the theory values for a given pair of wilson coefficients
	#returns array of returned prediction values
	prediction = pa.atpoint(**wcdict_tmp)
	prediction = dropLastValue(prediction)

	return [x[0] for x in prediction.values]



def getPPlot(pa, wcdict, noValues, xBins, yBins, obs, err, wilco1 = "qq3_i33i", wilco2 = "qq1_i33i"):
	#function to generated p values and output them as a table and as a plot

	wcdict_tmp = wcdict.copy()

	i = 0
	j = 0
	chi2Array = np.zeros((noValues, noValues))
	for i, x in enumerate(xBins):
		for j, y in enumerate(yBins):
			wcdict_tmp[wilco1] = x
			wcdict_tmp[wilco2] = y
			pred = np.array(getPred(pa, wcdict_tmp))
			chi2Value = chi2_local(obs, pred, err)
			chi2Array[j,i] = chi2Value

	chi2Array = chi2Array - chi2Array.min()

	#calculate p
	pArray = chi2.cdf(chi2Array, 2)

	np.savetxt("analysis_tables/pArray_{}_{}.txt".format(wilco1, wilco2), pArray)
	#print(chi2Array)
	#print(pArray)

	#plot
	fig, ax = plt.subplots()

	#cont = ax.contourf(xBins, yBins, chi2Array)
	cont = ax.contourf(xBins, yBins, pArray, [0.0, 0.68, 0.95, 0.99], cmap='viridis_r') 
	fig.colorbar(cont, ax=ax)
	plt.xlabel(wilco1.replace("_", r"\_"))
	plt.ylabel(wilco2.replace("_", r"\_"))
	#plt.show()
	fig.savefig("contour_plots2/contours_p_{}_{}.png".format(wilco1, wilco2))
	plt.close(fig)

	return pArray

def dropLastValue(df):
	#function to drop last value of each observable for closure
	last_index = df.index[0]
	prev_index = ["",""]
	toDrop = []

 
	for index in df.index:
		if index[0] != last_index[0]:
			toDrop += [prev_index]
			last_index = index
		prev_index = index

	for index in toDrop:
		df = df.drop(index)

	return df
		


def setupChi2(afPath = '../../HDF/CMS_2018_I1662081.h5', numvalues=25, range=(-5,5)):
	#load in from hdf and setup general variables
	af = AnalysisFrame.from_hdf(afPath)
	pa = PredictionArray(af.xr)
	#fh = FitHandler(pa)

	
	#get data values and error
	pseudoDat = pa.atpoint()
	pseudoDat = dropLastValue(pseudoDat)
	obs = np.array([x[0] for x in pseudoDat.values])
	err = np.array([x[1] for x in pseudoDat.values])

	#initialise a dictionary of wilco values
	wcdict_tmp = {}
	wcnames = pa.wilcos
	for name in wcnames:
		wcdict_tmp[name] = 0.
	
	#plotting noValues*noValues with given x and y ranges
	noValues = numvalues
	xBins = np.linspace(range[0], range[1], noValues)
	yBins = np.linspace(range[0], range[1], noValues)

	return	pa, wcdict_tmp, noValues, xBins, yBins, obs, err

if __name__ == "__main__":

	pa, wcdict, noValues, xBins, yBins, obs, err = setupChi2("../../HDF/ATLAS_2019_I1707015.h5", numvalues=5, range=(-3,3))
	#pa, wcdict, noValues, xBins, yBins, obs, err = setupChi2(numvalues=40, range=(-3,3))


	getPPlot(pa, wcdict, noValues, xBins, yBins, obs, err)
