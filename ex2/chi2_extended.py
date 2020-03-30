#imports
#numpy
import numpy as np
#matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

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
	Rotation gives anticlockwise rotation from positive x axis to major axis
	Assumes a != c - ie it is not a circle
	 """
	a, b, c = ellipseConsts[0], ellipseConsts[1]/2, ellipseConsts[2]
	if b == 0:
		if a < c:
			angle = 0
		else:
			angle = 0.5 * np.pi
	else:
		if a < c:
			angle = 0.5 * np.arctan(2*b/(a-c))
		else:
			angle = 0.5 * ( np.pi + np.arctan(2*b/(a-c)) )
	
	return angle


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
		return minLen, majLen, True
	else:
		return majLen, minLen, False

def getAxCen(ellipseConsts):
	""" Returns lenth of the two axes of the ellipse given the ellipse constans from fitting
	Assumes a != c - ie it is not a circle
	 """
	a, b, c = ellipseConsts[0], ellipseConsts[1]/2, ellipseConsts[2]
	d, f = ellipseConsts[3]/2, ellipseConsts[4]/2

	x = (c*d - b*f) / (b*b - a*c)
	y = (a*f - b*d) / (b*b - a*c)

	return x, y

def getCoordsInt(number, lenNarrow, lenWide, xCen, yCen, angle):
	posInt = np.linspace(lenNarrow, lenWide, number)
	negInt = np.linspace(-1*lenWide, -1*lenNarrow, number)
	
	#yPos = np.full(len(posInt), yCen)
	yPos = np.zeros(len(posInt))
	yNeg = np.full(len(negInt), -1 * yCen)

	xRotPos, yRotPos = rotatePoints(posInt, yPos, angle)
	xRotNeg, yRotNeg = rotatePoints(negInt, yNeg, angle)


	return xRotPos+xCen, yRotPos+yCen, xRotNeg+xCen, yRotNeg+yCen


def getCoordsFromEllipse(xBins, yBins, pArray, wilco1="qq3_i33i", wilco2="qq1_i33i"):
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
	lenWideMaj, lenWideMin, swapFlag = getAxLen(aWide)

	lenNarrowMaj, lenNarrowMin, swapFlag = getAxLen(aNarrow)
	if swapFlag == True:
		theta = theta + np.pi/2

	xCen, yCen = getAxCen(aNarrow)

	#get positive interpolation xValues
	numPoints = 50
	xPos, yPos, xNeg, yNeg = getCoordsInt( 
		numPoints, lenNarrowMaj, lenWideMaj, xCen, yCen, theta )


	majPoints = np.linspace(lenNarrowMaj, lenWideMaj, numPoints)

	#get for minor axis
	xPosMin, yPosMin, xNegMin, yNegMin = getCoordsInt( 
		numPoints, lenNarrowMin, lenWideMin, xCen, yCen, theta-np.pi/2 )
	minPoints = np.linspace(lenNarrowMin, lenWideMin, numPoints)


	fig.colorbar(cont, ax=ax)
	ax.scatter(xPos, yPos, c = 'g', marker='x')
	ax.scatter(xNeg, yNeg, c = 'g', marker='x')
	ax.scatter(xPosMin, yPosMin, c = 'r', marker='x')
	ax.scatter(xNegMin, yNegMin, c = 'r', marker='x')
	plt.xlabel(wilco1.replace("_",""))
	plt.ylabel(wilco2.replace("_",""))
	#plt.show()
	fig.savefig("contour_plots/rotationtest_{}_{}".format(wilco1.replace("_",""), wilco2.replace("_","")))
	plt.close()

	np.savetxt( "analysis_tables/int_vals_maj_{}_{}.txt".format(wilco1, wilco2), np.array( [ 
			majPoints, xPos, yPos, majPoints, xNeg, yNeg, minPoints, 
			xPosMin, xNegMin, minPoints, yPosMin, yNegMin ] ) )

if __name__ == "__main__":
	#calculate p
	pArray = np.loadtxt("pArray.txt")
	wilco1 = "qq3_i33i"
	wilco2 = "qq1_i33i"

	#plotting noValues*noValues with given x and y ranges
	noValues = 40
	xBins = np.linspace(-3, 3, noValues)
	yBins = np.linspace(-3, 3, noValues)

	getCoordsFromEllipse(xBins, yBins, pArray, wilco1, wilco2)

