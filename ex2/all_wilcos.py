from chi2_initial_values import setupChi2, getPPlot
from chi2_extended import getCoordsFromEllipse
from analysis import analyseWilcos
import matplotlib.pyplot as plt

fh, wcdict, noValues, xBins, yBins, obs, err = setupChi2()

checkedWilcos = []
for wilco1 in wcdict:
    checkedWilcos += [wilco1]

    for wilco2 in wcdict:
        if wilco2 not in checkedWilcos:
            print ("({},{})".format(wilco1, wilco2))
            pArray = getPPlot(fh, wcdict, noValues, xBins, yBins, obs, err, wilco1, wilco2)
            try:
                getCoordsFromEllipse(xBins, yBins, pArray, wilco1, wilco2)
                analyseWilcos(wilco1, wilco2)
            except:
                print("Not elliptical corr")
                plt.close()
            
            del pArray

        