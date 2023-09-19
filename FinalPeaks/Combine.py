import sys
import os
import numpy
import math

files = os.listdir("./Models_Decon")
if(len(files) == 0):
	print("No Files")
	exit()
starting = "./Models/" + files[0]
models = []
#found = []
count = 0
dim = 2
lines = []

for x in files:
	if(x != ".DS_Store"):
		starting = "./Models_Decon/" + x
		#print(starting)
		f = open(starting, "r")
		stuff = numpy.array(f.readlines())
		f.close()
		#lines = []
		#print(stuff)
		for x in stuff:
			lines.append((x.strip('\n')).split(","))
                                
#print(lines)

openfile = open("Combined_Results.csv", "w")
for x in lines:
    toPrint = ""
    for y in x:
    	toPrint = toPrint + str(y) + ","
    toPrint = toPrint[0:len(toPrint)-1]
    toPrint = toPrint + "\n"
    openfile.write(toPrint)
openfile.close()
