import numpy as np
import scipy as sp
from topfitter.analysis.frames import PredictionArray, AnalysisFrame
from topfitter.fitting import FitHandler
from numpy import corrcoef
from sklearn.feature_selection import mutual_info_regression
import matplotlib.pyplot as plt
import yoda


from chi2 import chi2_local, getPred

yoda.plotting.setup_mpl()

def analyse(wilco1="uu_i33i", wilco2="uu_ii33", afPath='../../HDF/CMS_2018_I1662081.h5'):

    #initialise topfitter objects
    af = AnalysisFrame.from_hdf(afPath)
    pa = PredictionArray(af.xr)
    fh = FitHandler(pa)

    #get data values and error
    obs = np.array([x[0] for x in fh.reference.values])
    #err = np.array([x[1] for x in fh.reference.values])
    print(fh.reference)

    #initialise dict for wilson coefficients
    wcdict = {}
    wcnames = pa.wilcos
    for name in wcnames:
        wcdict[name] = 0.

    #load in analysis from file
    loadArr = np.loadtxt("int_vals_maj.txt")

    #initialise loop variables
    #pred = np.zeros((loadArr[0].size, len(obs)))
    i = 0

    mutInfs = np.zeros((4,len(obs)))
    for lens, wc1Values, wc2Values in zip(loadArr[::3], loadArr[1::3], loadArr[2::3]):
        j = 0
        
        pred = np.zeros((wc1Values.size, len(obs)))

        # loop through x y values and get bin predictions there
        for wc1,wc2 in zip(wc1Values,wc2Values):
            wcdict[wilco1] = wc1
            wcdict[wilco2] = wc2

            binValues = np.array(getPred(fh, wcdict))
            
            pred[j] = binValues
             
            j += 1

        mutInfs[i] = mutual_info_regression(pred, lens)

        i += 1


    return mutInfs , fh.reference


ret, ref = analyse()

length = len(ret[0])
h = yoda.Histo1D(length, 0, length-1, "foo")

for i in range(length):
    h.fillBin(i, ret[0][i])
#yoda.plotting.mk_figaxes_1d()
fig, ax = yoda.plotting.plot_hist_1d(h)
plt.xticks(range(length), ref.index, rotation='vertical')
#plt.margins(0.5)
#plt.subplots_adjust(bottom=0.5)
plt.savefig('mutual_info')
plt.show()


#fig, ax = plt.subplots()

#plt.show()