import yoda
import matplotlib.pyplot as plt


def plotHist(ax, bins, values, label="", colour='b', linestyle='-'):
    #function to plot a stepped histogram given bins and bin values
	ax.plot(bins,
		np.append(values, values[-1]),
		ds='steps-post', label=label, ls=linestyle, color=colour)


def readFile(path):
    
    inList = yoda.read(path, asdict=False)

    for name in inList:
        print( name )

    return inList

def comp(data, sim):
    for dataplot, simplot in zip(data, sim):
        fig, ax = plt.subplot(1, 1)
        plotHist(ax, simplot.bins, simplot.values, label="Simulated result from Rivet")
        fig.plot()


if __name__ == "__main__":
    dataList = readFile("data.yoda")
    dataSL = dataList[:3]
    dataDL = dataList[3:]
    simDL = readFile("plots_DL.yoda")[3:]
    simSL = readFile("plots_SL.yoda")[:3]