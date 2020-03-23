import numpy as np
import yoda
from numpy import corrcoef
from sklearn.feature_selection import mutual_info_regression
import matplotlib.pyplot as plt
from topfitter.analysis.frames import PredictionArray, AnalysisFrame
from topfitter.fitting import FitHandler
from chi2 import getPred

def plotHist(ax, bins, values, label="", colour='b', linestyle='-'):
    #function to plot a stepped histogram given bins and bin values
	ax.plot(bins,
		np.append(values, values[-1]),
		ds='steps-post', label=label, ls=linestyle, color=colour)

def plotMut(mutInfReg, indexRef, uniqRef, name, w1, w2):
    #function to plot the mutual info score with appropriate labelling
    length = len(mutInfReg[0])
    bins = range(0, length+1, 1)
    fig, ax = plt.subplots()
    plotHist(ax, bins, mutInfReg[0])
    plotHist(ax,bins,mutInfReg[1], colour='r')
    plt.xticks(indexRef, uniqRef, fontsize=10, rotation=0, ha='left')
    #plt.subplots_adjust(bottom=0.2)
    plt.xlabel("Bins of observables")
    plt.ylabel("Mutual info regression score")
    plt.title("Mutual info for the {} axis of ellipse and for wilcos {} and {}".format(name, w1, w2),
               fontsize=12)
    plt.savefig('mutual_info_{}'.format(name))
    plt.show()
    plt.close(fig)

def analyse(wilco1="uu_i33i", wilco2="uu_ii33", afPath='../../HDF/CMS_2018_I1662081.h5'):
    #get the mutual info score given the list of points in the file

    #initialise topfitter objects
    af = AnalysisFrame.from_hdf(afPath)
    pa = PredictionArray(af.xr)
    fh = FitHandler(pa)

    #get data values and erroraf 
    obs = np.array([x[0] for x in fh.reference.values])
    err = np.array([x[1] for x in fh.reference.values])

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

    #array to store mutual info values for each pos/neg axis
    mutInfs = np.zeros((4,len(obs)))
    for lens, wc1Values, wc2Values in zip(loadArr[::3], loadArr[1::3], loadArr[2::3]):
        
        pred = np.zeros((wc1Values.size, len(obs)))

        # loop through x y values and get bin predictions there
        for j, (wc1,wc2) in enumerate(zip(wc1Values,wc2Values)):
            wcdict[wilco1] = wc1
            wcdict[wilco2] = wc2

            binValues = np.array(getPred(fh, wcdict))

            pred[j] = np.random.normal(binValues, err)



        #fig, ax = plt.subplots()
        #print (pred.shape)
        #plotHist(ax, range(len(wc1Values)+1),pred[:,13])
        #plt.show()
        #ax.cla()
        #fig, ax = plt.subplots()
        #plotHist(ax, range(len(obs)+1), pred[13])
        #plt.show()
        #plt.close(fig)

        mutInfs[i] = mutual_info_regression(pred, lens)

        i += 1


    return mutInfs , fh.reference.index

def main():
    wilco1 = "qq3_i33i"
    wilco2 = "qq1_i33i"
    mutInfReg, ref = analyse(wilco1=wilco1, wilco2=wilco2)

    #get unique copies of each index
    uniqRef = ref.levels[0]
    #get index of first appearance of each
    refFull = ref.get_level_values(0).tolist()

    indexRef = [refFull.index(refPart) for refPart in uniqRef]

    #format values in uniqRef for display
    #   replaces _ with \_, removes all the values before the variable name,
    #   and drops the normed from the end
    uniqRef = [r'\_'.join(x.replace('_', r'\_').split('/')[1].split(r'\_')[:-1]) for x in uniqRef]

    #plot on histogram
    plotMut(mutInfReg[0:2], indexRef, uniqRef, "maj", 
            wilco1.replace('_', r'\_'), wilco2.replace('_', r'\_'))
    plotMut(mutInfReg[2:4], indexRef, uniqRef, "min", 
            wilco1.replace('_', r'\_'), wilco2.replace('_', r'\_'))


    #fig, ax = plt.subplots()

    #plt.show()
    return 0

if __name__ == "__main__":
    main()