import sys
import os
import numpy
import math

def fitness(model):
    return(model[6])

dim = int(sys.argv[1])
a1 = sys.argv[2]
a2 = sys.argv[3]
#print(len(sys.argv))
if(dim > 2):
    if(len(sys.argv) > 4):
        a3 = sys.argv[4]
    else:
        print("Nucleus 3 missing")
if(dim > 3):
    if(len(sys.argv) > 5):
        a3 = sys.argv[5]
    else:
        print("Nucleus 4 missing")
        
files = os.listdir("./Models_Decon")
models = []
found = []
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
        for x in dataset:
            models.append(x)
        
        
        
openfile = open("Results.peaks", "w")
toPrint = "# Number of dimensions " + str(dim) + "\n"
toPrint = toPrint + "#FORMAT xeasy3D_LW\n"
if(a1 == "H"):
    toPrint = toPrint + "#INAME 1       1H\n"
if(a1 == "C"):
    toPrint = toPrint + "#INAME 1      13C\n"
if(a1 == "N"):
    toPrint = toPrint + "#INAME 1      15N\n"
if(a2 == "H"):
    toPrint = toPrint + "#INAME 2       1H\n"
if(a2 == "C"):
    toPrint = toPrint + "#INAME 2      13C\n"
if(a2 == "N"):
    toPrint = toPrint + "#INAME 2      15N\n"
if(dim > 2):
    if(a3 == "H"):
        toPrint = toPrint + "#INAME 3       1H\n"
    if(a3 == "C"):
        toPrint = toPrint + "#INAME 3      13C\n"
    if(a3 == "N"):
        toPrint = toPrint + "#INAME 3      15N\n"
if(dim > 3):
    if(a4 == "H"):
        toPrint = toPrint + "#INAME 4       1H\n"
    if(a4 == "C"):
        toPrint = toPrint + "#INAME 4      13C\n"
    if(a4 == "N"):
        toPrint = toPrint + "#INAME 4      15N\n"
#print(toPrint)
openfile.write(toPrint)

models.sort(key=fitness)
toPrint2 = ""
for x in range(1, len(models)+1):
    if(x < 10):
        toPrint2 = "     "
    if((x > 9) and (x < 100)):
        toPrint2 = "    "
    if(x > 99 | x < 1000):
        toPrint2 = "   "
    if(len(str(round(float(models[x-1][1]), 3))) == 3):
        toPrint2 = toPrint2 + str(x) + "   " + str(round(float(models[x-1][1]), 3)) + str(0) + str(0)
    if(len(str(round(float(models[x-1][1]), 3))) == 4):
        toPrint2 = toPrint2 + str(x) + "   " + str(round(float(models[x-1][1]), 3)) + str(0)
    if(len(str(round(float(models[x-1][1]), 3))) == 5):
        toPrint2 = toPrint2 + str(x) + "   " + str(round(float(models[x-1][1]), 3))
    if(len(str(round(float(models[x-1][1]), 3))) == 6):
        toPrint2 = toPrint2 + str(x) + "  " + str(round(float(models[x-1][1]), 3))
    if(len(str(round(float(models[x-1][1]), 3))) == 7):
        toPrint2 = toPrint2 + str(x) + " " + str(round(float(models[x-1][1]), 3))
        
    if(len(str(round(float(models[x-1][0]), 3))) == 6):
        toPrint2 = toPrint2 + " " +  str(round(float(models[x-1][0]), 3)) + str(0)
    else:
        toPrint2 = toPrint2 + " " +  str(round(float(models[x-1][0]), 3))
    #the last few columns are: colour code, user defined type of spectrum observed, volume, percent uncertainty of volume, integration method, and the rest is unused?
    toPrint2 = toPrint2 + " " +  "2 U          0.000e+00  0.00e+00  -   0 0 0 0"
    #toPrint2 = toPrint2 + " " +  "2 U          " + 0.000e+00 + "  " + "0.00e+00" " -   0 0 0 0"
    toPrint2 = toPrint2 + "\n";
    #print(toPrint2)
    
    #so far leaves linewidths and volumes blank. Can add in later if I feel like it.
    openfile.write(toPrint2)
openfile.close()