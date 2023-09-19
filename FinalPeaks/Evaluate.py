import sys
import os
import numpy
import math

#this only works for 2D. nvm i fixed it

dim = 2
peaksFlag = 0;
orderFlag = 0;
fitness = 0;
AAdims = []
AAdims.append('NA')
file1 = sys.argv[1]
if(len(sys.argv) > 2):
    peaksFlag = int(sys.argv[2])
if(len(sys.argv) > 3):
    orderFlag = int(sys.argv[3])
f = open(file1,"r")
stuff = f.readlines()
f.close()
TruthList = []
temp = []
temp.append(0.0)
for x in stuff:
    temp = []
    if(peaksFlag == 0):
        #print(((x.strip('\n'))).split("\t\t"))
        line = ((x.strip('\n'))).split("\t\t")
        #print(line)
        TruthList.append(line)
    if(peaksFlag == 1):
        line = ((x.strip('\n'))).split(" ")
        #print(line)
        for a in range(0,len(line)):
            if(line[a] != ''):
                temp.append(line[a])
        #print(temp)
        TruthList.append(temp)
        #TruthList.append(((x.strip('\n'))).split(" "))
        if(TruthList[len(TruthList)-1][0] == '#'):
            if(TruthList[len(TruthList)-1][1] == 'Number'):
                if(TruthList[len(TruthList)-1][2] == 'of'):
                    if(TruthList[len(TruthList)-1][3] == 'dimensions'):
                        dim = int(TruthList[len(TruthList)-1][4])
                        while(len(AAdims) < dim):
                            AAdims.append('NA')
        if(TruthList[len(TruthList)-1][0] == '#INAME'):
            if(TruthList[len(TruthList)-1][1] == '1'):
                        AAdims[0] = TruthList[len(TruthList)-1][2]
        if(len(AAdims) > 1):
            if(TruthList[len(TruthList)-1][0] == '#INAME'):
                if(TruthList[len(TruthList)-1][1] == '2'):
                        AAdims[1] = TruthList[len(TruthList)-1][2]
        if(len(AAdims) > 2):
            if(TruthList[len(TruthList)-1][0] == '#INAME'):
                if(TruthList[len(TruthList)-1][1] == '3'):
                        AAdims[2] = TruthList[len(TruthList)-1][2]
        if(len(AAdims) > 3):
            if(TruthList[len(TruthList)-1][0] == '#INAME'):
                if(TruthList[len(TruthList)-1][1] == '4'):
                        AAdims[3] = TruthList[len(TruthList)-1][2]
    if(peaksFlag == 2):
        #print(((x.strip('\n'))).split("\t\t"))
        line = ((x.strip('\n'))).split("  ")
        #print(line)
        if(len(line) > 1):
            while('' in line):
                line.pop(line.index(''))
        #print(line)
        if(len(line) > 1): 
        	if( line[0] != "Assignment"):
        		#print((line[1],line[2]))
        		TruthList.append((line[1],line[2]))
                        
#print(AAdims)
#print(dim)

#print(TruthList)

found = []

if(peaksFlag == 0):
    AAnum = 0
    aminoA = ""
    n_shift = 0
    h_shift = 0
    TruthList2 = []
    for x in TruthList:
        if(x[0][0] == ';'):
            TruthList2.append([AAnum, aminoA, h_shift, n_shift])
            AAnum = 0
            aminoA = ""
            n_shift = 0
            h_shift = 0
        if(x[0][0] != ';'):
            if(x[0][0] != '#'):
                if(x[0][0] != 'N'):
                    if(x[0][0] != 'H'):
                        if(x[0][0] != 'C'):
                            thing = x[0].split(" ")
                            AAnum = int(thing[0])
                            aminoA  = thing[1]
                if(x[0][0] == 'N'):
                    n_shift = float(x[1])
                if(x[0][0] == 'H'):
                    h_shift = float(x[1])

    files = os.listdir("./Models_Decon")
    starting = "./Models/" + files[0]
    models = []
    #found = []
    count = 0

    for x in files:
        if(x != ".DS_Store"):
            starting = "./Models_Decon/" + x
            #print(starting)
            f = open(starting, "r")
            stuff = numpy.array(f.readlines())
            f.close()
            lines = []
            for x in stuff:
                lines.append((x.strip('\n')).split(","))

            dataset = numpy.array(lines)
            sizeIncrease = 2
            for y in dataset:
                #print(y)
                if(dim == 2):
                    for w in TruthList2:
                        if(w[0] != 'sequence'):
                            distance = ((float(w[3]) - float(y[0]))/(sizeIncrease*float(y[2])))*((float(w[3]) - float(y[0]))/(sizeIncrease*float(y[2]))) + ((float(w[2]) - float(y[1]))/(sizeIncrease*float(y[3])))*((float(w[2]) - float(y[1]))/(sizeIncrease*float(y[3])))
                            #this might be backwards. dif1 might need to be compared to models[count][3] and dif2 might need to be compared to models[count][2]. The linewidth difference makes a large impact on possible candidates in the hydrogen dimension (cause its so small compared to the others)
                            if(distance < 1):
                                found.append((w, y, distance))

if(peaksFlag == 1):
    TruthList2 = []
    temp = []
    for x in range(0, dim):
        temp.append(0.0)
    for x in TruthList:
        temp = []
        for z in range(0, dim):
            temp.append(0.0)
        if((x[0] != '#') & (x[0] != '#INAME')):
            for y in range(0, dim):
                temp[y] = x[1+y]
            TruthList2.append(temp)
        
    totalTruth = len(TruthList2)
    
    files = os.listdir("./Models_Decon")
    starting = "./Models/" + files[0]
    models = []
    count = 0
    hillsWithOutMatches = 0
    models = 0
    flag = 0
    hillsFound = 0
    left = []
    tooMany = 0
    
    lines = []
    for x in files:
        if(x != ".DS_Store"):
            starting = "./Models_Decon/" + x
            f = open(starting, "r")
            stuff = numpy.array(f.readlines())
            f.close()
            for x in stuff:
                lines.append((x.strip('\n')).split(","))
            
    models = len(lines)
    hillsWithOutMatches = len(TruthList2)
    
    
    howmanyGrounthTruthFound = 0;
    
    fitness = 0
    dataset = lines
    for x in range(0, len(dataset)):
        left.append(x)
        fitness = fitness + float(dataset[x][len(dataset[x])-2])
    sizeIncrease = 0.5
    for w in TruthList2:
        flag = 0
        current = (w, dataset[0], 1000000)
        distance = 10000000
        distance2 = 100000
        for y in dataset:
            if(orderFlag == 0):
                if(dim == 2):
                    if((float(y[2]) != 0) & (float(y[3]) != 0)):
                        distance = ((float(w[1]) - float(y[0]))/(sizeIncrease*float(y[2])))*((float(w[1]) - float(y[0]))/(sizeIncrease*float(y[2]))) + ((float(w[0]) - float(y[1]))/(sizeIncrease*float(y[3])))*((float(w[0]) - float(y[1]))/(sizeIncrease*float(y[3])))
                        distance2 = (float(w[1]) - float(y[0]))*(float(w[1]) - float(y[0])) + (float(w[0]) - float(y[1]))*(float(w[0]) - float(y[1]))
                if(dim == 3):
                    if((float(y[3]) != 0) & (float(y[4]) != 0) & (float(y[5]) != 0)):
                        distance = ((float(w[2]) - float(y[0]))/(sizeIncrease*float(y[3])))*((float(w[2]) - float(y[0]))/(sizeIncrease*float(y[3]))) + ((float(w[1]) - float(y[1]))/(sizeIncrease*float(y[4])))*((float(w[1]) - float(y[1]))/(sizeIncrease*float(y[4]))) + ((float(w[0]) - float(y[2]))/(sizeIncrease*float(y[5])))*((float(w[0]) - float(y[2]))/(sizeIncrease*float(y[5])))
                        distance2 = (float(w[2]) - float(y[0]))*(float(w[2]) - float(y[0])) + (float(w[1]) - float(y[1]))*(float(w[1]) - float(y[1])) + (float(w[0]) - float(y[2]))*(float(w[0]) - float(y[2]))
            if(orderFlag == 1):
                if(dim == 2):
                    if((float(y[2]) != 0) & (float(y[3]) != 0)):
                        distance = ((float(w[0]) - float(y[0]))/(sizeIncrease*float(y[3])))*((float(w[0]) - float(y[0]))/(sizeIncrease*float(y[3]))) + ((float(w[1]) - float(y[1]))/(sizeIncrease*float(y[2])))*((float(w[1]) - float(y[1]))/(sizeIncrease*float(y[2])))
                        distance2 = (float(w[0]) - float(y[0]))*(float(w[0]) - float(y[0])) + (float(w[1]) - float(y[1]))*(float(w[1]) - float(y[1]))
                if(dim == 3):
                    if((float(y[3]) != 0) & (float(y[4]) != 0) & (float(y[5]) != 0)):
                        distance = ((float(w[0]) - float(y[0]))/(sizeIncrease*float(y[5])))*((float(w[0]) - float(y[0]))/(sizeIncrease*float(y[5]))) + ((float(w[1]) - float(y[1]))/(sizeIncrease*float(y[4])))*((float(w[1]) - float(y[1]))/(sizeIncrease*float(y[4]))) + ((float(w[2]) - float(y[2]))/(sizeIncrease*float(y[3])))*((float(w[2]) - float(y[2]))/(sizeIncrease*float(y[3])))
                        distance2 = (float(w[0]) - float(y[0]))*(float(w[0]) - float(y[0])) + (float(w[1]) - float(y[1]))*(float(w[1]) - float(y[1])) + (float(w[2]) - float(y[2]))*(float(w[2]) - float(y[2]))
            if((distance < 1) & (distance2 < 0.5)):
                if(distance < current[2]):
                    current = (w, y, distance)
        if(current[2] < 1):
            found.append(current)
            #left.pop(left.index(dataset.index(current[1])))
            flag = 1
            if(dataset.index(current[1]) in left):
                #found.append(current)
                howmanyGrounthTruthFound = howmanyGrounthTruthFound +1
                left.pop(left.index(dataset.index(current[1])))
                #flag = 1
            else:
                tooMany = tooMany + 1
        if(flag == 1):
            hillsWithOutMatches = hillsWithOutMatches - 1
            hillsFound = hillsFound  + 1
        # if(flag != 1):
#             if(dim == 2):
#                 found.append(((0,0), y, 1000000))
#             if(dim == 3):
#                 found.append(((0,0,0), y, 1000000))
    
    
    HillsMissing = len(left)
    
    #print(left)
    
    for b in left:
        if(dim == 2):
            found.append(((0,0), dataset[b], 1000000))
        if(dim == 3):
            found.append(((0,0,0), dataset[b], 1000000))
        
    
    #output = "Hills found: " + str(HillsFound) + " Total ground truth hills: " + str(totalTruth) + " Hills not found: " + str(HillsMissing) + " Hills without assignments: " + str(hillsWithOutMatches)
    output = "Models Assigned: " + str(hillsFound) + " Total ground truth hills: " + str(totalTruth) + " Hills not found: " + str(hillsWithOutMatches) + " Hills without assignments: " + str(HillsMissing)
    print(output)
    print("total number of models: " + str(models))
    print("GroundTruth Assigned: " + str(howmanyGrounthTruthFound))
    thing = hillsFound - (totalTruth - hillsWithOutMatches)
    otherthing = 0
    if(models > 0):
        otherthing = fitness/models
    print("extra assignments: " + str(thing))
    print("Average RMSE: " + str(otherthing))
    
if(peaksFlag == 2):
    TruthList2 = []
    if("w1" in TruthList[0]):
        if("w2" in TruthList[0]):
            if("w3" in TruthList[0]):
                dim = 3
            else:
                dim = 2
        else:
            dim = 1
    
    for x in TruthList:
        if((x[0] != "Assignment") & (len(TruthList[0]) == len(x))):
            TruthList2.append(x)
            
    #print(len(TruthList2))
    #print(TruthList2)
    totalTruth = len(TruthList2)
    
    files = os.listdir("./Models_Decon")
    starting = "./Models/" + files[0]
    models = []
    count = 0
    hillsWithOutMatches = 0
    models = 0
    flag = 0
    hillsFound = 0
    left = []
    tooMany = 0
    
    lines = []
    for x in files:
        if(x != ".DS_Store"):
            starting = "./Models_Decon/" + x
            f = open(starting, "r")
            stuff = numpy.array(f.readlines())
            f.close()
            for x in stuff:
                lines.append((x.strip('\n')).split(","))
            
    models = len(lines)
    hillsWithOutMatches = len(TruthList2)
    
    
    howmanyGrounthTruthFound = 0;
    
    fitness = 0
    dataset = lines
    for x in range(0, len(dataset)):
        left.append(x)
        fitness = fitness + float(dataset[x][len(dataset[x])-2])
    sizeIncrease = 0.5
    
    for w in TruthList2:
        flag = 0
        current = (w, dataset[0], 1000000)
        distance = 10000000
        distance2 = 100000
        for y in dataset:
            if(orderFlag == 0):
                if(dim == 2):
                    if((float(y[2]) != 0) & (float(y[3]) != 0)):
                        distance = ((float(w[1]) - float(y[0]))/(sizeIncrease*float(y[2])))*((float(w[1]) - float(y[0]))/(sizeIncrease*float(y[2]))) + ((float(w[0]) - float(y[1]))/(sizeIncrease*float(y[3])))*((float(w[0]) - float(y[1]))/(sizeIncrease*float(y[3])))
                        distance2 = (float(w[1]) - float(y[0]))*(float(w[1]) - float(y[0])) + (float(w[0]) - float(y[1]))*(float(w[0]) - float(y[1]))
                if(dim == 3):
                    if((float(y[3]) != 0) & (float(y[4]) != 0) & (float(y[5]) != 0)):
                        distance = ((float(w[2]) - float(y[0]))/(sizeIncrease*float(y[3])))*((float(w[2]) - float(y[0]))/(sizeIncrease*float(y[3]))) + ((float(w[0]) - float(y[1]))/(sizeIncrease*float(y[4])))*((float(w[0]) - float(y[1]))/(sizeIncrease*float(y[4]))) + ((float(w[1]) - float(y[2]))/(sizeIncrease*float(y[5])))*((float(w[1]) - float(y[2]))/(sizeIncrease*float(y[5])))
                        distance2 = (float(w[2]) - float(y[0]))*(float(w[2]) - float(y[0])) + (float(w[0]) - float(y[1]))*(float(w[0]) - float(y[1])) + (float(w[1]) - float(y[2]))*(float(w[1]) - float(y[2]))
            if(orderFlag == 1):
                if(dim == 2):
                    if((float(y[2]) != 0) & (float(y[3]) != 0)):
                        distance = ((float(w[0]) - float(y[0]))/(sizeIncrease*float(y[3])))*((float(w[0]) - float(y[0]))/(sizeIncrease*float(y[3]))) + ((float(w[1]) - float(y[1]))/(sizeIncrease*float(y[2])))*((float(w[1]) - float(y[1]))/(sizeIncrease*float(y[2])))
                        distance2 = (float(w[0]) - float(y[0]))*(float(w[0]) - float(y[0])) + (float(w[1]) - float(y[1]))*(float(w[1]) - float(y[1]))
                if(dim == 3):
                    if((float(y[3]) != 0) & (float(y[4]) != 0) & (float(y[5]) != 0)):
                        distance = ((float(w[0]) - float(y[0]))/(sizeIncrease*float(y[5])))*((float(w[0]) - float(y[0]))/(sizeIncrease*float(y[5]))) + ((float(w[1]) - float(y[1]))/(sizeIncrease*float(y[4])))*((float(w[1]) - float(y[1]))/(sizeIncrease*float(y[4]))) + ((float(w[2]) - float(y[2]))/(sizeIncrease*float(y[3])))*((float(w[2]) - float(y[2]))/(sizeIncrease*float(y[3])))
                        distance2 = (float(w[0]) - float(y[0]))*(float(w[0]) - float(y[0])) + (float(w[1]) - float(y[1]))*(float(w[1]) - float(y[1])) + (float(w[2]) - float(y[2]))*(float(w[2]) - float(y[2]))
            if((distance < 1) & (distance2 < 0.5)):
                if(distance < current[2]):
                    current = (w, y, distance)
        if(current[2] < 1):
            found.append(current)
            #left.pop(left.index(dataset.index(current[1])))
            flag = 1
            if(dataset.index(current[1]) in left):
                #found.append(current)
                howmanyGrounthTruthFound = howmanyGrounthTruthFound +1
                left.pop(left.index(dataset.index(current[1])))
                #flag = 1
            else:
                tooMany = tooMany + 1
        if(flag == 1):
            hillsWithOutMatches = hillsWithOutMatches - 1
            hillsFound = hillsFound  + 1
        # if(flag != 1):
        #     if(dim == 2):
        #         found.append(((0,0), y, 1000000))
        #     if(dim == 3):
        #         found.append(((0,0,0), y, 1000000))
    
    
    HillsMissing = len(left)
    
    #print(left)
    
    for b in left:
        if(dim == 2):
            found.append(((0,0), dataset[b], 1000000))
        if(dim == 3):
            found.append(((0,0,0), dataset[b], 1000000))
    
    #output = "Hills found: " + str(HillsFound) + " Total ground truth hills: " + str(totalTruth) + " Hills not found: " + str(HillsMissing) + " Hills without assignments: " + str(hillsWithOutMatches)
    output = "Models Assigned: " + str(hillsFound) + " Total ground truth hills: " + str(totalTruth) + " Hills not found: " + str(hillsWithOutMatches) + " Hills without assignments: " + str(HillsMissing)
    print(output)
    print("total number of models: " + str(models))
    print("GroundTruth Assigned: " + str(howmanyGrounthTruthFound))
    thing = hillsFound - (totalTruth - hillsWithOutMatches)
    otherthing = 0
    if(models > 0):
        otherthing = fitness/models
    print("extra assignments: " + str(thing))
    print("Average RMSE: " + str(otherthing))

openfile = open("Results.txt", "w")
for x in found:
    toPrint = ""
    for y in x[0]:
        toPrint = toPrint + str(y) + " "
    toPrint = toPrint + "  Model: "
    for y in x[1]:
        toPrint = toPrint + y + " "
    toPrint = toPrint + str(x[2]) + " "
    toPrint = toPrint + "\n"
    openfile.write(toPrint)
openfile.close()
