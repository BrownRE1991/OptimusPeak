#!/usr/bin/env /usr/local/bin/python3

import sys
import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import warnings
warnings.filterwarnings("ignore")

#uses nmrGlue to extract spectra.
def readNMR( filename ):
    dic,data = ng.pipe.read(filename)
    data = np.flip(data, axis=0)
    data = np.flipud(data)
    return dic, data

#takes the spectrum, the upper and lower thresholds, and the value to collapse to ("zero"), and returns the array where all values between the upper and lower thresholds are converted to "zero". This function does the grunt works after the threshold has been determined. You can use it with any threshold you want.
def applyFixedThreshold(specData, threshold1, threshold2, zero):
    th1 = 0
    th2 = 0
    if(threshold1 >= threshold2):
        th1 = threshold1
        th2 = threshold2
    if(threshold1 < threshold2):
        th1 = threshold2
        th2 = threshold1
    if(threshold1 == threshold2):
        print("Threshold Error")
    tempData = np.ravel(specData.copy())
    tempData[(tempData <= th1) & (tempData >= th2)] = zero
    if(specData.ndim == 2):
        tempData = tempData.reshape(specData.shape[0],specData.shape[1])
    if(specData.ndim == 3):
        tempData = tempData.reshape(specData.shape[0],specData.shape[1], specData.shape[2])
    return tempData
    
#def showHist(data, fileName, count, log):
#    hist = cv2.calcHist([data.ravel()], [0], None, [256], [0, 256])
#    plt.plot(hist, color='k')
#    plt.clf()
#    plt.plot(hist, color='k')
#    plt.xlim([0, 256])
#    plt.xlabel("Intensity")
#    plt.ylabel("Frequency")
#    plt.title('Histogram of Intensities' + str(count))
#    if log:
#        plt.yscale('log')
#    plt.clf()
#    plt.close()
    
#def showHist2(data, fileName, count, log):
#    hist = data.ravel()
#    plt.hist(hist, bins=5000)
#    plt.clf()
#    plt.hist(hist, bins=5000)
#    plt.xlabel("Intensity")
#    plt.ylabel("Frequency")
#    plt.title('Histogram of Intensities')
#    if log:
#        plt.yscale('log')
#    plt.clf()
#    plt.close()

#model used to fit the noise 
def gaussian(x, amp, cen, wid):
    bottom = 2*(wid**2)
    front = amp/(np.sqrt(2*math.pi) * wid)
    start = x - cen
    step2 = -(start)**2
    step3 = step2/bottom
    if(step3.size != 1):
        for x in range(step3.size):
            step3[x] = math.exp(step3[x])
    else:
        step3 = math.exp(step3)
    step5 = front*step3
    return (np.array(step5))

#calculates the thresholding by fitting a gaussian curve to the histogram of intensities, normalizing, and determining the threshold to be the center of the points whose intensiy is between 10^-6 and 10^-12. This process is repeated if there are less than 20 bins between the upper and lower thresholds.
def modelHistGauss2(data, row, count, man):
    name1 = "Threshold/Plots/" + "Gaus2" + row + "_HistPlusModel.png"
    name2 = "Threshold/Plots/" + "Gaus2" + row + "_HistPlusModelLog.png"
    workingData = data.ravel()
    m = np.mean(workingData)
    sig = np.sqrt(np.var(workingData))
    med = np.median(workingData)
    x = np.linspace(min(workingData), max(workingData), 5000)
    n, binEdge, patches = plt.hist(workingData, bins =5000)
    n = n/workingData.size
    plt.clf()
    maxH = max(n)
    params = [maxH, med, sig]
    pars, cov = curve_fit(f=gaussian, xdata=x, ydata=n, p0=params)
    model = gaussian(x, pars[0], pars[1], pars[2])
    minThresh = (10*max(n), 0)
    maxThresh = (10*max(n), 0)
    maxModel = max(model)
    index = np.where((model < 9.05763694396763e-06) & (model > 2.5471292322959e-12))
    if(len(index[0]) == 0):
        index = np.where(model > 0)
    
    maxStart = np.where((binEdge[index[0]] >= 0))
    if(len(maxStart[0]) > 0):
        if(len(maxStart[0])%2 > 0):
            minMaxStart = np.median(binEdge[index[0][maxStart[0]]])
        else:
            middle = int(len(binEdge[index[0][maxStart[0]]])/2)-1
            minMaxStart = binEdge[index[0][maxStart[0][middle]]]
        maxStart2 = np.where(binEdge == minMaxStart)
        maxThresh = (binEdge[maxStart2[0]], maxStart2[0])
    
    minStart = np.where((binEdge[index[0]] < 0))
    if(len(minStart[0]) > 0):
        if(len(minStart[0])%2 > 0):
            minMaxStart = np.median(binEdge[index[0][minStart[0]]])
        else:
            middle = int(len(binEdge[index[0][minStart[0]]])/2)-1
            minMaxStart = binEdge[index[0][middle]]
        minStart2 = np.where(binEdge == minMaxStart)
        minThresh = (binEdge[minStart2[0]], minStart2[0])
    else:
        minThresh = maxThresh
        
    if((len(maxStart[0]) == 0) & (len(minStart[0]) > 0)):
        maxThresh = minThresh
    
    if(maxThresh[1] - minThresh[1] > 99):
        if((maxThresh[1] != 0)):
            x = np.linspace(min(workingData), max(workingData), 5000)
            plt.hist(workingData, bins =5000)
            plt.plot(x, model*workingData.size)
            plt.axvline(x=binEdge[maxThresh[1]][0],color='red')
            plt.xlabel("Intensity")
            plt.ylabel("Frequency")
            plt.title('Histogram of Intensities' + str(count))
            plt.clf()
            plt.close()
            
            x = np.linspace(min(workingData), max(workingData), 5000)
            n, binEdge, patches = plt.hist(workingData, bins =5000)
            plt.yscale('log')
            plt.plot(x, model*workingData.size)
            plt.axvline(x=binEdge[maxThresh[1]][0],color='red')
            plt.ylim(10**(-0.5), 2*max(n))
            plt.xlabel("Intensity")
            plt.ylabel("Frequency")
            plt.title('Histogram of Intensities' + str(count))
            plt.clf()
            plt.close()
            newdata = applyFixedThreshold(data, binEdge[maxThresh[1]][0], binEdge[minThresh[1]-1][0], 0.0)
            return binEdge[maxThresh[1]][0], newdata
    
    if((maxThresh[1] != 0)):
        lowercut = minThresh[1]-40
        uppercut = maxThresh[1]+40
        dif = uppercut - lowercut
        newset = workingData[((workingData <= binEdge[uppercut]) & (workingData >= binEdge[lowercut]))]
        if(dif == 80):
            bins = 10000
        else:
            bins = 900
        m = np.mean(workingData)
        sig = np.sqrt(np.var(newset))
        med = np.median(newset)
        x = np.linspace(min(newset), max(newset), int(bins))
        n, binEdge, patches = plt.hist(newset, int(bins))
        n = n/workingData.size
        plt.clf()
        maxH = max(n)
        params = [maxH, med, sig]
        pars, cov = curve_fit(f=gaussian, xdata=x, ydata=n, p0=params)
        model = gaussian(x, pars[0], pars[1], pars[2])
        minThresh = (10*max(n), 0)
        maxThresh = (10*max(n), 0)
        maxModel = max(model)
        index = np.where((model < 9.05763694396763e-06) & (model > 2.5471292322959e-12))
        maxStart = np.where((binEdge[index[0]] > 0))
        if(len(maxStart[0]) > 0):
            if(len(maxStart[0])%2 > 0):
                minMaxStart = np.median(binEdge[index[0][maxStart[0]]])
            else:
                middle = int(len(binEdge[index[0][maxStart[0]]])/2)-1
                minMaxStart = binEdge[index[0][maxStart[0][middle]]]
            maxStart2 = np.where(binEdge == minMaxStart)
            maxThresh = (binEdge[maxStart2[0]], maxStart2[0])
    
        minStart = np.where((binEdge[index[0]] < 0))
        if(len(minStart[0]) > 0):
            if(len(minStart[0])%2 > 0):
                minMaxStart = np.median(binEdge[index[0][minStart[0]]])
            else:
                middle = int(len(binEdge[index[0][minStart[0]]])/2)-1
                minMaxStart = binEdge[index[0][middle]]
            minStart2 = np.where(binEdge == minMaxStart)
            minThresh = (binEdge[minStart2[0]], minStart2[0])
        else:
            minThresh = maxThresh
            
        if((len(maxStart[0]) == 0) & (len(minStart[0]) > 0)):
            maxThresh = minThresh

        if(maxThresh[1] != minThresh[1]):
            x = np.linspace(min(newset), max(newset), int(bins))
            n, binEdge, patches = plt.hist(newset, bins =int(bins))
            plt.plot(x, model*newset.size)
            plt.axvline(x=binEdge[maxThresh[1]],color='red')
            plt.xlabel("Intensity")
            plt.ylabel("Frequency")
            plt.title('Histogram of Intensities' + str(count))
            plt.clf()
            plt.close()
            x = np.linspace(min(newset), max(newset), int(bins))
            plt.hist(newset, bins =int(bins))
            plt.yscale('log')
            plt.plot(x, model*newset.size)
            plt.axvline(x=binEdge[maxThresh[1]],color='red')
            plt.ylim(10**(-0.5), 2*max(n))
            plt.xlabel("Intensity")
            plt.ylabel("Frequency")
            plt.title('Histogram of Intensities' + str(count))
            plt.clf()
            plt.close()
            newdata = applyFixedThreshold(data, binEdge[maxThresh[1]][0], binEdge[minThresh[1]-1][0], 0.0)
            return binEdge[maxThresh[1]][0], newdata
            
    x = np.linspace(min(workingData), max(workingData), 10000)
    n, binEdge, patches = plt.hist(workingData, bins =10000)
    n = n/newset.size
    plt.clf()
    maxH = max(n)
    params = [maxH, med, sig]
    pars, cov = curve_fit(f=gaussian, xdata=x, ydata=n, p0=params)
    model = gaussian(x, pars[0], pars[1], pars[2])
    minThresh = (10*max(n), 0)
    maxThresh = (10*max(n), 0)
    maxModel = max(model)
    index = np.where(model == maxModel)
    for x in range(model.size):
        if(x <= index[0][0]):
            if model[x] >= n[x]:
                if(model[x] > 1):
                    if(model[x] < minThresh[0]):
                        minThresh = (model[x], x)
        if(model.size-x-1 >= index[0][0]):
            if(model[model.size-x-1] >= n[x]):
                if(model[model.size-x-1] > 1):
                    if(model[model.size-x-1] < maxThresh[0]):
                        maxThresh = (model[model.size-x-1], (model.size-x-1))
    if(maxThresh[1] != minThresh[1]):
        newdata = applyFixedThreshold(data, binEdge[maxThresh[1]], binEdge[minThresh[1]], 0.0)
        return binEdge[maxThresh[1]], newdata
    else:
        print("Error")
    newdata = applyFixedThreshold(data, binEdge[maxThresh[1]], binEdge[maxThresh[1] - 1], 0.0)
    return binEdge[maxThresh[1]], newdata
    
#creates the final output. Can be run in binary or text output. It is used in this program as binary output to save space.
def printPrintList(printList, outfile, how, flag):
    openfile = open(outfile, how)
    if(flag == 0):
    	for x in printList:
    		line = ""
    		for y in x:
    			if(y != 0):
    				line = line + str(y) + " "
    		line = line +"\n"
    		openfile.write(line)
    if(flag == 1):
    	np.asarray(printList, dtype=np.dtype('d')).tofile(openfile)
    openfile.close()
    

#main algorithm. Reads in the spectrum using nmrGlue, extracts information, thresholds is, and outputs thresholded datapoints to binary file.Currently works for 2 and 3D spectra, but can be extended to more. The part that takes the longest here is printing to file. Python sucks at triple loops, but I am not quite sure how to get around the nmrGlue data structure. I might tackle it in future work.
if(len(sys.argv) < 2):
    print("ERROR. Not enough Arguments")
    exit();
file1 = sys.argv[1]

dicCC,dataCC_orig = readNMR(file1)
GausThresh2, GausData2 = modelHistGauss2(dataCC_orig, "test", 1, 0)

if(dataCC_orig.ndim == 2):
    uc_dim1 = ng.pipe.make_uc(dicCC, GausData2, dim=0)
    ppm_dim1 = uc_dim1.ppm_scale()
    ppm_dim1_0, ppm_dim1_1 = uc_dim1.ppm_limits()
    uc_dim2 = ng.pipe.make_uc(dicCC, GausData2, dim=1)
    ppm_dim2 = uc_dim2.ppm_scale()
    ppm_dim2_0, ppm_dim2_1 = uc_dim2.ppm_limits()
    
    printList = []
    for x in range(0, len(ppm_dim1)):
        for y in range(0, len(ppm_dim2)):
            if((ppm_dim1[x] >= ppm_dim1_1) & (ppm_dim1[x] <= ppm_dim1_0)):
                if((ppm_dim2[y] >= ppm_dim2_1) & (ppm_dim2[y] <= ppm_dim2_0)):
                    string1 = str(ppm_dim1[x])+ "ppm"
                    string2 = str(ppm_dim2[y])+ "ppm"
                    if(GausData2[uc_dim1(string1)][uc_dim2(string2)] != 0):
                        printList.append([ppm_dim1[x], ppm_dim2[y], GausData2[uc_dim1(string1)][uc_dim2(string2)], 0])
    printPrintList(printList, "Threshold/Outlist.bin", "wb", 1)
if(dataCC_orig.ndim == 3):
    uc_dim1 = ng.pipe.make_uc(dicCC, GausData2, dim=0)
    ppm_dim1 = uc_dim1.ppm_scale()
    ppm_dim1_0, ppm_dim1_1 = uc_dim1.ppm_limits()
    uc_dim2 = ng.pipe.make_uc(dicCC, GausData2, dim=1)
    ppm_dim2 = uc_dim2.ppm_scale()
    ppm_dim2_0, ppm_dim2_1 = uc_dim2.ppm_limits()
    uc_dim3 = ng.pipe.make_uc(dicCC, GausData2, dim=2)
    ppm_dim3 = uc_dim3.ppm_scale()
    ppm_dim3_0, ppm_dim3_1 = uc_dim3.ppm_limits()
    
    #print(ppm_dim3_0, ppm_dim3_1)
    #print(ppm_dim2_0, ppm_dim2_1)
    #print(ppm_dim1_0, ppm_dim1_1)
    
    printList = []
    #outfile1 = "Threshold/Outlist.txt"
    outfile2 = "Threshold/Outlist.bin"
    openfile = open(outfile1, "w")
    openfile.close()
    openfile = open(outfile2, "wb")
    openfile.close()
    for x in range(0, len(ppm_dim1)):
        for y in range(0, len(ppm_dim2)):
            for w in range(0, len(ppm_dim3)):
                if((ppm_dim1[x] >= ppm_dim1_1) & (ppm_dim1[x] <= ppm_dim1_0)):
                    if((ppm_dim2[y] >= ppm_dim2_1) & (ppm_dim2[y] <= ppm_dim2_0)):
                        if((ppm_dim3[w] >= ppm_dim3_1) & (ppm_dim3[w] <= ppm_dim3_0)):
                            string1 = str(ppm_dim1[x])+ "ppm"
                            string2 = str(ppm_dim2[y])+ "ppm"
                            string3 = str(ppm_dim3[w])+ "ppm"
                            if((GausData2[uc_dim1(string1)][uc_dim2(string2)][uc_dim3(string3)] != 0)):
                            	printList.append([ppm_dim1[x], ppm_dim2[y], ppm_dim3[w], GausData2[uc_dim1(string1)][uc_dim2(string2)][uc_dim3(string3)], 0])
                            	if(len(printList) == 10000):
                                	printPrintList(printList, outfile2, "ab", 1)
                                	while(len(printList) > 0):
                                		printList.pop()
    if(len(printList) > 0):
    	printPrintList(printList, outfile2, "ab", 1)
    	while(len(printList) > 0):
        	printList.pop()
print(GausThresh2)

