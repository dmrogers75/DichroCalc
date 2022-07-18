#!/usr/bin/python
# DMRogers 2021-02-22. Replace Ham_ij elements
import sys
import numpy as np

input_fname1 = sys.argv[1]
input_fname2 = sys.argv[2]

inputH = np.loadtxt(input_fname1)
print(inputH)

print(len(inputH))

inputR = np.loadtxt(input_fname2)
print(inputR)

#i = 0
#j = 2
#print(inputH[[i],[j]])

outputH = inputH.copy()

i1 = 0
j1 = 2
outputH[[i1],[j1]] = inputR[[0],[0]]
outputH[[j1],[i1]] = inputR[[0],[0]]
i2 = 0
j2 = 3
outputH[[i2],[j2]] = inputR[[0],[1]]
outputH[[j2],[i2]] = inputR[[0],[1]]
i3 = 1
j3 = 2
outputH[[i3],[j3]] = inputR[[1],[0]]
outputH[[j3],[i3]] = inputR[[1],[0]]
i4 = 1
j4 = 3
outputH[[i4],[j4]] = inputR[[1],[1]]
outputH[[j4],[i4]] = inputR[[1],[1]]

if(len(inputH)>4):
    istop = len(inputH)
    print(istop)
    istart = j4 + 1
    print(istart)
    for i in range(istart, istop, 2):
        i1 = i1 + 2
        j1 = j1 + 2
        outputH[[i1],[j1]] = inputR[[0],[0]]
        outputH[[j1],[i1]] = inputR[[0],[0]]
        i2 = i2 + 2
        j2 = j2 + 2
        outputH[[i2],[j2]] = inputR[[0],[1]]
        outputH[[j2],[i2]] = inputR[[0],[1]]
        i3 = i3 + 2
        j3 = j3 + 2
        outputH[[i3],[j3]] = inputR[[1],[0]]
        outputH[[j3],[i3]] = inputR[[1],[0]]
        i4 = i4 + 2
        j4 = j4 + 2
        outputH[[i4],[j4]] = inputR[[1],[1]]
        outputH[[j4],[i4]] = inputR[[1],[1]]

print(outputH)

output_fname = str(input_fname1)+'.Rep.txt'
np.savetxt(output_fname,outputH,fmt='%15.6f',delimiter=' ')

