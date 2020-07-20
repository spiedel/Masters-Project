import numpy as np
import yoda
from numpy import corrcoef
from sklearn.feature_selection import mutual_info_regression
import matplotlib.pyplot as plt
import matplotlib
from topfitter.analysis.frames import PredictionArray, AnalysisFrame
from topfitter.fitting import FitHandler
from chi2_initial_values import getPred
matplotlib.rcParams.update({'errorbar.capsize': 2})

def plotHist(ax, bins, values, label="", colour='b', linestyle='-'):
    #function to plot a stepped histogram given bins and bin values
	ax.plot(bins,
		np.append(values, values[-1]),
		ds='steps-post', label=label, ls=linestyle, color=colour)

def plotMut(mutInfReg, err, indexRef, uniqRef, name, w1, w2):
    #function to plot the mutual info score with appropriate labelling
    length = len(mutInfReg[0])
    bins = np.array(range(0, length+1, 1))
    fig, ax = plt.subplots()
    #positive
    plotHist(ax, bins, mutInfReg[0], label='Positive values')
    binCent = (bins[1:] + bins[:-1])/2.0
    ax.errorbar(binCent, mutInfReg[0], yerr=err[0], fmt='none', color='b')


    #negative
    plotHist(ax,bins,mutInfReg[1], label = 'Negative Values', colour='r')
    ax.errorbar(binCent, mutInfReg[1], yerr=err[1], fmt='none', color='r')
    ax.set_ylim(bottom=0)
    plt.xticks(indexRef, uniqRef, fontsize=10, rotation=0, ha='left')
    plt.xlabel("Bins of observables")
    plt.ylabel("Mutual info regression score")
    plt.legend()
    plt.savefig('analysis_plots/handdone/mutual_info_{}_{}_{}.png'.format(name, w1, w2))
    #plt.savefig('ref_{}.png'.format(name))
    #plt.show()
    plt.close(fig)

def plotDiff(mutInf1, mutInf2, err1, err2, indexRef, uniqRef, name, w1, w2):
    #function to plot the mutual info score with appropriate labelling
    length = len(mutInf1[0])
    bins = np.array(range(0, length+1, 1))
    fig, ax = plt.subplots()
    #positive
    plotHist(ax, bins, mutInf1[0]-mutInf2[0], label='Positive values')
    binCent = (bins[1:] + bins[:-1])/2.0
    ax.errorbar(binCent, mutInf1[0]-mutInf2[0], yerr=np.sqrt(err1[0]*err1[0] + err2[0]*err2[0]), fmt='none', color='b')

    #negative
    plotHist(ax,bins,mutInf1[1]-mutInf2[1], label = 'Negative Values', colour='r')
    ax.errorbar(binCent, mutInf1[1]-mutInf2[1], yerr=np.sqrt(err1[1]*err1[1] + err2[1]*err2[1]), fmt='none', color='r')
    plt.xticks(indexRef, uniqRef, fontsize=10, rotation=0, ha='left')
    plt.xlabel("Bins of observables")
    plt.ylabel("Mutual info regression score")
    plt.legend()
    plt.savefig('analysis_plots/handdone/mutual_info_diff_{}_{}_{}.png'.format(name, w1, w2))
    plt.close(fig)

def analyse(wcdict, obs, err, wilco1="qq3_i33i", wilco2="qq1_i33i"):
    #get the mutual info score given the list of points in the file

    #load in analysis from file
    loadArr = np.loadtxt("analysis_tables/int_vals_maj_{}_{}.txt".format(wilco1, wilco2))

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

            binValues = np.full(len(obs), j)#np.array(getPred(pa, wcdict))#
            

            pred[j] = np.random.normal(binValues, err)



        # fig, ax = plt.subplots()
        # print (pred.shape)
        # plotHist(ax, np.append(lens, lens[-1]),pred[:,0])
        # plt.xlabel("Length along the contour ellipse axis")
        # plt.ylabel("Bin Value")
        # plt.show()
        # ax.cla()
        # fig, ax = plt.subplots()
        # plotHist(ax, range(len(obs)+1), pred[13])
        # plt.show()
        # plt.close(fig)

        mutInfs[i] = mutual_info_regression(pred, lens)

        i += 1


    return mutInfs

def analyseWilcos(fh, wcdict, obs, err, wilco1 = "qq3_i33i", wilco2 = "qq1_i33i"):
    
    mutInfs = analyse(wcdict, obs, err, wilco1=wilco1, wilco2=wilco2)
    ref = fh.index
    mutInfs = [mutInfs]
    for i in range(15):
        mutInfs = np.append(mutInfs, [analyse(wilco1=wilco1, wilco2=wilco2)[0]], axis=0)

    
    mutInfReg = np.mean(mutInfs, axis=0)
    mutInfErr = np.std(mutInfs, axis=0)
    
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
    plotMut(mutInfReg[0:2], mutInfErr[0:2], indexRef, uniqRef, "maj", 
            wilco1, wilco2)
    plotMut(mutInfReg[2:4], mutInfErr[2:4], indexRef, uniqRef, "min", 
            wilco1, wilco2)

    #plt.show()
    return mutInfReg, mutInfErr, indexRef, uniqRef

if __name__ == "__main__":
    #TODO fix so can be run as main again
    mut1, err1, index, ref = analyseWilcos("uu_i33i", "uu_ii33")
    mut2, err2, index, ref = analyseWilcos()

    plotDiff(mut1, mut2, err1, err2, index, ref, "", "uu_i33i_uu_ii33", "qq3_i33i_qq1_i33i")
