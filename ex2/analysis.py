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

    af = AnalysisFrame.from_hdf(afPath)
    pa = PredictionArray(af.xr)
    fh = FitHandler(pa)

    #get data values and error
    obs = np.array([x[0] for x in fh.reference.values])
    err = np.array([x[1] for x in fh.reference.values])
    print(fh.reference)

    wcdict = {}
    wcnames = pa.wilcos
    for name in wcnames:
        wcdict[name] = 0.


    loadArr = np.loadtxt("int_vals_maj.txt")
    pred = np.zeros((4, loadArr[0].size, len(obs)))
    lenAlongAx = []
    i = 0

    ret = []
    for wc1Values, wc2Values in zip(loadArr[::2], loadArr[1::2]):
        j = 0
        #pred = np.zeros((len(wc1Values), len(obs)))
        #chitot = np.zeros(len(wc1Values))
        
        for wc1,wc2 in zip(wc1Values,wc2Values):
            wcdict[wilco1] = wc1
            wcdict[wilco2] = wc2

            binValues = np.array(getPred(fh, wcdict))
            
            pred[i][j] = binValues    
            j += 1

        i += 1

    
    ret = mutual_info_regression(pred, )

    return ret , fh.reference


ret, ref = analyse()

length = len(ret)
# h = yoda.Histo1D(len(ret), 0, len(ret)-1, "foo")

# for i in range(len(ret)):
#     h.fillBin(i, ret[i])
# #yoda.plotting.mk_figaxes_1d()
# fig, ax = yoda.plotting.plot_hist_1d(h)
plt.xticks(range(len(ret)), ref.index, rotation='vertical')
#plt.margins(0.5)
#plt.subplots_adjust(bottom=0.5)
plt.savefig('mutual info')
plt.show()


#fig, ax = plt.subplots()

#plt.show()