

from chi2 import setupChi2, getPPlot
from chi2_extended import getCoordsFromEllipse
from analysis import analyseWilcos

fh, wcdict, noValues, xBins, yBins, obs, err = setupChi2()

checkedWilcos = []
for wilco1 in wcdict:
    #remove this key from dict as no longer needs to be compared to 
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
