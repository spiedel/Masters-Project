import yoda
import matplotlib.pyplot as plt
import numpy as np


def plotHist(ax, bins, values, label="", colour='b', linestyle='-'):
    #function to plot a stepped histogram given bins and bin values
	ax.plot(bins, np.append(values, values[-1]),
		ds='steps-post', label=label, ls=linestyle, color=colour)


def readFile(path):
    
    inList = yoda.read(path, asdict=False)
    return inList

def comp(data, sim, labelsDict):
    for dataplot, simplot in zip(data, sim):
        fig, ax = plt.subplots()
        #set the plot specific parameters
        labels = labelsDict[dataplot.name()]
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        if ( labels[2] == True ):
            ax.set_yscale('log')

        #prepare arrays for plotting simulated values
        simBins = np.append(simplot.xMins(), simplot.xMax())
        simValues = simplot.yVals()
        simErr = simplot.yErrs()

        #plot simulated values
        plotHist(ax, simBins, simValues, label="Simulated result from Rivet")
        ax.errorbar(simplot.xMids(), simplot.yVals(), yerr=simErr, fmt='none', color='b')

        #find errors for data
        absErrs = dataplot.yErrs().T
        data = dataplot.yVals()
        errUp = absErrs[0]-data
        errDown = data-absErrs[1]

        #plot data
        ax.scatter(dataplot.xVals(), dataplot.yVals(), color='r', label='Data')
        ax.errorbar(dataplot.xVals(), data, yerr=[errUp, errDown], color='r', fmt='.')

        plt.legend()
        #save plot then close it for next run
        plt.savefig('check_plots/{}.png'.format(dataplot.name()))
        plt.close(fig)


if __name__ == "__main__":
    dataList = readFile("data.yoda")
    simList = readFile("out.yoda")

    h = yoda.Histo1D

    labelsDict = { 
        "d01-x01-y01":[r"$\mathrm{p}_\mathrm{T}(\gamma)$ [GeV]", 
            r"$\frac{1}{\sigma}\frac{\mathrm{d}\sigma}{\mathrm{d}_\mathrm{T}(\gamma)}$ [GeV]",
            True], 
        "d02-x01-y01":[r"$|\eta(\gamma)|$", 
                r"$\frac{1}{\sigma}\frac{\mathrm{d}\sigma}{\mathrm{d}|\eta(\gamma)|}$",
                False], 
        "d03-x01-y01":[r"$\Delta\mathrm{R}\left(\gamma,\mathrm{l}\right)$", 
            r"$\frac{1}{\sigma}\frac{\mathrm{d}\sigma}{\mathrm{d}\eta}$", 
            False], 
        "d04-x01-y01":[r"$\mathrm{p}_\mathrm{T}(\gamma)$ [GeV]",
            r"$\frac{1}{\sigma}\frac{\mathrm{d}\sigma}{\mathrm{d}_\mathrm{T}(\gamma)}$ [GeV]",
            True],
        "d05-x01-y01":[r"$|\eta(\gamma)|$", 
            r"$\frac{1}{\sigma}\frac{\mathrm{d}\sigma}{\mathrm{d}|\eta(\gamma)|}$",
            False], 
        "d06-x01-y01":[r"$\Delta\mathrm{R}\left(\gamma,\mathrm{l}\right)_{\mathrm{min}}$", 
            r"$\frac{1}{\sigma}\frac{\mathrm{d}\sigma}{\mathrm{d}\eta}$", 
            False],
        "d07-x01-y01":[r"$|\Delta\eta\left(\mathrm{l,l}\right)|$",
            r"$\frac{1}{\sigma}\frac{\mathrm{d}\sigma}{\mathrm{d}|\Delta\eta\left(\mathrm{l,l}\right)|}$",
            False], 
        "d08-x01-y01":[r"$\Delta\phi\left(\mathrm{l,l}\right)$",
            r"$\frac{1}{\sigma}\frac{\mathrm{d}\sigma}{\mathrm{d}\Delta\phi\left(\mathrm{l,l}\right)}$",
            False]
    }

    comp(dataList, simList, labelsDict)


    