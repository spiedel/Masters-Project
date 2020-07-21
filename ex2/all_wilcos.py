from chi2_initial_values import setupChi2, getPPlot
from chi2_extended import getCoordsFromEllipse
from analysis import analyseWilcos
import matplotlib.pyplot as plt

afPath = "../../HDF/ATLAS_2019_I1707015.h5"
fh, wcdict, noValues, xBins, yBins, obs, err = setupChi2(afPath, range=(-1,1), numvalues=20)

#uncomment below to check specific wilcos, with their value entered
#wcdict = {"G":0, "qq3_i33i":0}

checkedWilcos = []
for wilco1 in wcdict:
    checkedWilcos += [wilco1]

    for wilco2 in wcdict:
        if wilco2 not in checkedWilcos:
            print ("({},{})".format(wilco1, wilco2))
            pArray = getPPlot(fh, wcdict, noValues, xBins, yBins, obs, err, wilco1, wilco2)
            try:
                getCoordsFromEllipse(xBins, yBins, pArray, wilco1, wilco2)
                analyseWilcos(fh, wcdict, obs, err, wilco1, wilco2)
            except:
                print("Not elliptical corr")
                plt.close()
            
            del pArray

        
