#imports
#tf
from topfitter.analysis.frames import PredictionArray, AnalysisFrame
from topfitter.fitting import FitHandler
#numpy
import numpy as np
#matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
#scipy
from scipy.stats import chisquare, chi2
from scipy import stats

def chi2_local(observed, expected, error):
	#function to calculate chi2 from two arrays
	if len(observed) != len(expected):
		print("Inputted data doesn't match length")
		return 0.

	chi2_sum = 0.
	for i in range(len(observed)):
	# calculate each chi2 part
		chi2_part = (observed[i]-expected[i])**2/error[i]
		chi2_sum = chi2_sum + chi2_part

	return chi2_sum


def getPred(fh, wcdict):
	#funtion for predicting the theory values for a given pair of wilson coefficients
	#{'qq3_i33i', 'G', 'qq1_ii33', 'uu_i33i', 'qu1_33ii', 'qu8_33ii', 'qq1_i33i', 'qu8_ii33', 'qq3_ii33', 'uG_33', 'ud1_33ii', 'qd8_33ii', 'qu1_ii33', 'uu_ii33', 'qd1_33ii', 'ud8_33ii'}
	#returns
	prediction = fh.predict(**wcdict)#iqq3_i33i=wcValue2)#, uG_33=wcValue2)#
	return [x[0] for x in prediction.values]

def rotatePoints(x, y, angle):
	"""Function to rotate a point in x,y space clockwise around the axis by angle radians"""
	xnew = x * np.cos(angle) - y * np.sin(angle)
	ynew = x * np.sin(angle) + y * np.cos(angle)
	return xnew, ynew

def fitEllipse(x,y):
	"""Function that takes two points and fits an equation of an ellipse to it. 
	Source http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html """
	x = x[:,np.newaxis]
	y = y[:,np.newaxis]
	D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
	S = np.dot(D.T,D)
	C = np.zeros([6,6])
	C[0,2] = C[2,0] = 2; C[1,1] = -1
	E, V =  np.linalg.eig(np.dot(np.linalg.inv(S), C))
	n = np.argmax(np.abs(E))
	a = V[:,n]
	return a

def getAngle(ellipseConsts):
	""" Return rotation angle in radians of the fitted ellipse given the ellipse constants from fitting
	Rotation gives anticlockwise rotation from positive x axis to semi major axis
	Assumes a != c - ie it is not a circle
	 """
	a, b, c = ellipseConsts[0], ellipseConsts[1]/2, ellipseConsts[2]
	if b == 0:
		if a < c:
			return 0
		else:
			return 0.5 * np.pi
	else:
		if a < c:
			return 0.5 * np.arctan(2*b/(a-c))
		else:
			return 0.5 * ( np.pi + np.arctan(2*b/(a-c)) )


def getAxLen(ellipseConsts):
	""" Returns lenth of the two axes of the ellipse given the ellipse constants from fitting
	Assumes a != c - ie it is not a circle
	 """
	a, b, c = ellipseConsts[0], ellipseConsts[1]/2, ellipseConsts[2]
	d, f, g = ellipseConsts[3]/2, ellipseConsts[4]/2, ellipseConsts[5]

	#get equation parts
	denom = 2 * (a*f*f + c*d*d + g*b*b - 2*b*d*f - a*c*g)
	majNum = (b*b - a*c) * ( np.sqrt((a-c)**2 + 4*b*b) - (a+c) )
	minNum = (b*b - a*c) * ( -1 * np.sqrt((a-c)**2 + 4*b*b) - (a+c) )

	majLen = np.sqrt(denom/majNum)
	minLen = np.sqrt(denom/minNum)

	if majLen < minLen:
		return minLen, majLen
	else:
		return majLen, minLen

def getAxCen(ellipseConsts):
	""" Returns lenth of the two axes of the ellipse given the ellipse constans from fitting
	Assumes a != c - ie it is not a circle
	 """
	a, b, c = ellipseConsts[0], ellipseConsts[1]/2, ellipseConsts[2]
	d, f = ellipseConsts[3]/2, ellipseConsts[4]/2

	x = (c*d - b*f) / (b*b - a*c)
	y = (a*f - b*d) / (b*b - a*c)

	return x, y

def getCoordsInt(lenNarrow, lenWide, xCen, yCen, angle):
	posInt = np.linspace(lenNarrow, lenWide, 50)
	negInt = np.linspace(-1*lenWide, -1*lenNarrow, 50)
	print(xCen)
	
	#yPos = np.full(len(posInt), yCen)
	yPos = np.zeros(len(posInt))
	yNeg = np.full(len(negInt), -1 * yCen)

	xRotPos, yRotPos = rotatePoints(posInt, yPos, angle)
	xRotNeg, yRotNeg = rotatePoints(negInt, yNeg, angle)


	return xRotPos+xCen, yRotPos+yCen, xRotNeg, yRotNeg


wilco1 = "qq3_i33i"
wilco2 = "qq1_i33i"

#generate prediction array from data
#af = AnalysisFrame.from_hdf('/nfs/topfitter/sbrown/HDF/analyses/ATLAS_2017_I1604029/ATLAS_2017_I1604029.h5')
af = AnalysisFrame.from_hdf('../../HDF/CMS_2018_I1662081.h5')
pa = PredictionArray(af.xr)
fh = FitHandler(pa)

#get data values and error
obs = [x[0] for x in fh.reference.values]
err = [x[1] for x in fh.reference.values]
print(err)

wcdict = {}
wcnames = pa.wilcos
for name in wcnames:
	wcdict[name] = 0.


#print(pa.reference)
#print(fh.reference)
#print(pa.atpoint(G=1, qq3_i33i=1))
#print(chisquare(obs, [x[0] for x in pa.atpoint(G=3, qq3_i33i=3).values]))

#split by observable? use pa.reference.loc[level, 'value'] to get list instead

 
#plotting noValues*noValues with given x and y ranges
noValues = 40
xBins = np.linspace(-3, 3, noValues)
yBins = np.linspace(-3, 3, noValues)

i = 0
j = 0
chi2Array = np.zeros((noValues, noValues))
for x in xBins:
	for y in yBins:
		wcdict[wilco1] = x
		wcdict[wilco2] = y
		pred = getPred(fh, wcdict)
		chi2Value = chi2_local(obs, pred, err)
		chi2Array[i][j] = chi2Value
		j+=1
	j=0
	i+=1
#print(pred)

#print(chi2Array)
#Calculate delta chi2

chi2Array = chi2Array - chi2Array.min()

#calculate p
pArray = chi2.sf(chi2Array, 2)

#plot
fig, ax = plt.subplots()
#cont = ax.contourf(xBins, yBins, chi2Array)
cont = ax.contour(xBins, yBins, pArray, [0.05, 0.35]) 
#print(cont.collections[0].get_paths()[0].vertices)

wideCont = cont.allsegs[0][0]
narrowCont = cont.allsegs[1][0]
#print(wideCont)
#print(narrowCont)

xWide = np.array([x[0] for x in wideCont])
yWide = np.array([x[1] for x in wideCont])
xNarrow = np.array([x[0] for x in narrowCont])
yNarrow = np.array([x[1] for x in narrowCont])

aWide = fitEllipse(xWide, yWide)
aNarrow = fitEllipse(xNarrow, yNarrow)
theta = getAngle(aNarrow)
lenWideMaj, lenWideMin = getAxLen(aWide)
lenNarrowMaj, lenNarrowMin = getAxLen(aNarrow)
xCen, yCen = getAxCen(aNarrow)


print(lenWideMaj)
print(lenNarrowMaj)

#get positive interpolation xValues
xPos, yPos, xNeg, yNeg = getCoordsInt( lenNarrowMaj, lenWideMaj, xCen, yCen, theta)


fig.colorbar(cont, ax=ax)
ax.scatter(xPos, yPos, c = 'g', marker='o')
plt.xlabel(wilco1.replace("_",r"\_"))
plt.ylabel(wilco2.replace("_",r"\_"))
plt.show()
fig.savefig("contour_plots/rotationtest_{}_{}".format(wilco1.replace("_",""), wilco2.replace("_","")))
