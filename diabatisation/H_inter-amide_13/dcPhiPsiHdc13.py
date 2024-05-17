# DMRogers 2024-01-09. Modify DichroCalc Hamiltonian
import sys
import numpy as np
  
input_fname1 = sys.argv[1]

inputPhiPsi = np.loadtxt(input_fname1)
outputPhiPsi = inputPhiPsi.copy()

# For N residues (AAs) and M=N-1 peptide bonds
# Interaction of peptide bonds 1 and 2 defined by (phi,psi) of AA 2; interaction of bonds (M-1) and M by (phi,psi) of AA N-1

print("No. of (phi,psi) in input file = ", len(outputPhiPsi))

istart = 0
istop = len(outputPhiPsi)

for i in range(istart, istop, 1):
    j = 0
    k = 1
    print(outputPhiPsi[[i],[j]],outputPhiPsi[[i],[k]])

    pepBond = []
    pepPhi = []
    pepPsi = []
for i in range(istart, istop, 1):
    j = 0
    k = 1
    pepBond.append(i + 1)
    pepPhi.append(outputPhiPsi[[i],[j]])
    pepPsi.append(outputPhiPsi[[i],[k]])
    print(pepBond[i],pepPhi[i],pepPsi[i])

print("No. of (phi,psi) or residues = ", len(pepBond))


# Clean above to remove 1st Nter and Cter from each chain. 1st (360,x) and all (x,360)
# As consequence of chain borders and copy back

pepBondC = []
pepPhiC = []
pepPsiC = []
jc = 0
for i in range(istart+1, istop, 1): # start at istart+1 to omit 1st Nter (360,x)
    if pepPsi[i] == 360.0: # to omit each Cter (x,360) 
        i = i + 1
    else:
        jc = jc + 1
        pepBondC.append(jc)
        pepPhiC.append(pepPhi[i]) 
        pepPsiC.append(pepPsi[i]) 
        print (jc,pepPhi[i],pepPsi[i])

for i in range(0, len(pepBondC), 1):         # to check above is OK 
    print(pepBondC[i],pepPhiC[i],pepPsiC[i]) # 

print (len(pepBondC),len(pepPhiC),len(pepPsiC))
pepBond = pepBondC
pepPhi = pepPhiC
pepPsi = pepPsiC
print("No. of cleaned peptide bonds minus final peptide bond = ", len(pepBond))


# Diamide phi, psi

diam = ('2a','2b','2c','2d','2e','2f','2g','2h','2i','2j','2k','2l','2m')
dPhi = (180.0,-120.0,-60.0,-135.0,-120.0,-120.0,-120.0,-60.0,-74.0,-48.0,-60.0,-62.0,-75.0)
dPsi = (180.0, 180.0,180.0, 135.0, 120.0,  60.0,   0.0,  0.0, -4.0,-57.0,-60.0,-41.0,145.0)

print(len(diam),len(dPhi),len(dPsi))

istartD = 0
istopD = len(diam)
for i in range(istartD, istopD, 1):
    print(diam[i],dPhi[i],dPsi[i])


# Use Euclidean distance to find closest (phi,psi). p2 is diamide, p1 is protein. p1 - p2

istart = 0
istop = len(pepBond)
pepDdist = []
pepDdLab = []
for i in range(istart, istop, 1):
    disCp = 1000.0
    disKe = 1000.0
    diamL,diamCutOff = "2x",float(sys.argv[3])  # diamCutOff is cut-off for Euclidean distance. Input argument 3
    for j in range(istartD, istopD, 1):
        distance = np.sqrt( (pepPhi[i] - float(dPhi[j]))**2 + (pepPsi[i] - float(dPsi[j]))**2 )
        if i >= 0:
            disCp = distance
        if diamCutOff == 0.0: # if diamCutOff = 0.0 then ignore
            if abs(disCp) < abs(disKe):
                disKe = disCp
                diamL = diam[j]
            print (pepBond[i], pepPhi[i], pepPsi[i], diam[j], dPhi[j], dPsi[j], distance)
        if diamCutOff != 0.0:
            if abs(disCp) < abs(disKe) and abs(disCp) < diamCutOff: # use ab initio 2x if disCp > diamCutOff
                disKe = disCp
                diamL = diam[j]
            print (pepBond[i], pepPhi[i], pepPsi[i], diam[j], dPhi[j], dPsi[j], distance)
    pepDdist.append(disKe)
    pepDdLab.append(diamL)
    print (pepBond[i], pepDdist[i], pepDdLab[i])


# Replace Hamiltonian elements for peptide bond N coupling with N+1,
# H_i,j = (2xN)-1,(2xN)+1 (2xN)-1,(2xN)+2 2xN,(2xN)+1 2xN,(2xN)+2

istart = 0
istop = len(pepBond)
for i in range(istart, istop, 1):
    npnpR = 2 * pepBond[i] - 1
    npnpC = 2 * pepBond[i] + 1
    npppR = 2 * pepBond[i] - 1
    npppC = 2 * pepBond[i] + 2
    ppnpR = 2 * pepBond[i]
    ppnpC = 2 * pepBond[i] + 1
    ppppR = 2 * pepBond[i]
    ppppC = 2 * pepBond[i] + 2 
    print (pepBond[i], pepDdLab[i], npnpR, npnpC, npppR, npppC, ppnpR, ppnpC, ppppR, ppppC)

# Diamide inter-amide couplings (cm-1) 
dCoup2a = ( -25.5787, 14.8052, 16.3827, 388.6942 )
dCoup2b = ( 341.4442, 131.2434, -920.9323, -203.1671 )
dCoup2c = ( -124.7431, 146.0569, -721.2296, 681.3303 )
dCoup2ds = ( 109.3404, -279.2747, -243.6219, -89.3696 )  # standard diabatic 2d
dCoup2dms4 = ( 109.3404, 279.2747, -243.6219, 89.3696 )  # modified diabatic 2d (monomer TDM state 4 signs)
dCoup2dms24 = ( 109.3404, 279.2747, 243.6219, -89.3696 ) # modified diabatic 2d (monomer TDM states 2,4 signs)
dCoup2dma2e = ( 16.9637, 238.7406, 242.5565, -37.2743 )  # modified diabatic 2d (adiabatic TDM signs as 2e)
dCoup2e = ( 339.3501, -27.0207, 284.8980, -906.3860 )
dCoup2f = ( 33.8472, -384.6068, 576.1425, 159.4446 )
dCoup2g = ( -11.9690, -16.8514, 541.5561, -131.2202 )
dCoup2h = ( -270.2800, -8.7818, 225.7205, -1701.0504 )
dCoup2i = ( 288.5190, -147.6903, 230.0022, 284.2714 )
dCoup2j = ( 96.5062, -808.1508, 0.1650, -191.3471 )
dCoup2k = ( -20.6389, -1707.4763, 266.7429, 172.9444 )
dCoup2l = ( -358.5142, -659.7897, -131.4717, -362.7687 )
dCoup2m = ( 80.8986, 356.3688, 145.0900, 29.2226 )

# Read in Hamiltonian to be modified and copy it
input_fname2 = sys.argv[2]
inputH = np.loadtxt(input_fname2)
print (len(inputH))
print (inputH)
outputH = inputH.copy()

# Select which 2d couplings to use (optional)
if len(sys.argv) < 5:
    couplings2d = 1  # default
if len(sys.argv) == 5:
    couplings2d = int(sys.argv[4])  # input argument 4
print (couplings2d)
if couplings2d == 0:
    dCoup2d = dCoup2ds
    suffix2d = "2d"
if couplings2d == 1:
    dCoup2d = dCoup2dms4
    suffix2d = "2ds4"
if couplings2d == 2:
    dCoup2d = dCoup2dms24
    suffix2d = "2ds24"
if couplings2d == 3:
    dCoup2d = dCoup2dma2e
    suffix2d = "2da2e"

# Loop over pepBond and diamide couplings and substitute inter-amide couplings
for i in range(istart, istop, 1):
    npnpR = 2 * pepBond[i] - 1
    npnpC = 2 * pepBond[i] + 1
    npppR = 2 * pepBond[i] - 1
    npppC = 2 * pepBond[i] + 2
    ppnpR = 2 * pepBond[i]
    ppnpC = 2 * pepBond[i] + 1
    ppppR = 2 * pepBond[i]
    ppppC = 2 * pepBond[i] + 2 
    for j in range(istartD, istopD, 1):
        if pepDdLab[i] == diam[j]:
            print (pepBond[i], pepDdLab[i], diam[j])
            if pepDdLab[i] == "2a":
                print ("Substituting 2a  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2a)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2a[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2a[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2a[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2a[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2a[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2a[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2a[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2a[3]
            if pepDdLab[i] == "2b":
                print ("Substituting 2b  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2b)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2b[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2b[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2b[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2b[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2b[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2b[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2b[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2b[3]
            if pepDdLab[i] == "2c":
                print ("Substituting 2c  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2c)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2c[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2c[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2c[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2c[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2c[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2c[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2c[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2c[3]
            if pepDdLab[i] == "2d":
                print ("Substituting 2d  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2d)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2d[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2d[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2d[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2d[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2d[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2d[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2d[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2d[3]
            if pepDdLab[i] == "2e":
                print ("Substituting 2e  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2e)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2e[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2e[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2e[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2e[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2e[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2e[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2e[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2e[3]
            if pepDdLab[i] == "2f":
                print ("Substituting 2f  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2f)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2f[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2f[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2f[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2f[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2f[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2f[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2f[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2f[3]
            if pepDdLab[i] == "2g":
                print ("Substituting 2g  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2g)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2g[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2g[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2g[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2g[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2g[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2g[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2g[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2g[3]
            if pepDdLab[i] == "2h":
                print ("Substituting 2h  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2h)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2h[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2h[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2h[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2h[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2h[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2h[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2h[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2h[3]
            if pepDdLab[i] == "2i":
                print ("Substituting 2i  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2i)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2i[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2i[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2i[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2i[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2i[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2i[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2i[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2i[3]
            if pepDdLab[i] == "2j":
                print ("Substituting 2j  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2j)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2j[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2j[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2j[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2j[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2j[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2j[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2j[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2j[3]
            if pepDdLab[i] == "2k":
                print ("Substituting 2k  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2k)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2k[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2k[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2k[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2k[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2k[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2k[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2k[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2k[3]
            if pepDdLab[i] == "2l":
                print ("Substituting 2l  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2l)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2l[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2l[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2l[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2l[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2l[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2l[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2l[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2l[3]
            if pepDdLab[i] == "2m":
                print ("Substituting 2m  H_i,j ", npnpR, npnpC," ", npppR, npppC," ", ppnpR, ppnpC," ", ppppR, ppppC," and lower triangle")
                print (dCoup2m)
                outputH[[npnpR-1],[npnpC-1]] = dCoup2m[0] # Upper triangle -1 for indexing starting at 0
                outputH[[npppR-1],[npppC-1]] = dCoup2m[1]
                outputH[[ppnpR-1],[ppnpC-1]] = dCoup2m[2]
                outputH[[ppppR-1],[ppppC-1]] = dCoup2m[3]
                outputH[[npnpC-1],[npnpR-1]] = dCoup2m[0] # Lower triangle -1 for indexing starting at 0
                outputH[[npppC-1],[npppR-1]] = dCoup2m[1]
                outputH[[ppnpC-1],[ppnpR-1]] = dCoup2m[2]
                outputH[[ppppC-1],[ppppR-1]] = dCoup2m[3]
print (outputH)

# Output modified Hamiltonian to file
#output_fname = str(input_fname2)+'.Hdc13.txt'
output_fname = str(input_fname2)+'.Hdc13_'+str(diamCutOff)+'_'+str(suffix2d)+'.txt'
np.savetxt(output_fname,outputH,fmt='%15.6f',delimiter=' ')

# Count number of diamide (phi,psi) peptide bonds in system. 13 (phi,psi)
t2a, t2b, t2c, t2d, t2e, t2f, t2g, t2h, t2i, t2j, t2k, t2l, t2m, t2x = 0,0,0,0,0,0,0,0,0,0,0,0,0,0
for i in range(istart, istop, 1):
    if pepDdLab[i] == "2a":
        t2a = t2a + 1
    if pepDdLab[i] == "2b":
        t2b = t2b + 1
    if pepDdLab[i] == "2c":
        t2c = t2c + 1
    if pepDdLab[i] == "2d":
        t2d = t2d + 1
    if pepDdLab[i] == "2e":
        t2e = t2e + 1
    if pepDdLab[i] == "2f":
        t2f = t2f + 1
    if pepDdLab[i] == "2g":
        t2g = t2g + 1
    if pepDdLab[i] == "2h":
        t2h = t2h + 1
    if pepDdLab[i] == "2i":
        t2i = t2i + 1
    if pepDdLab[i] == "2j":
        t2j = t2j + 1
    if pepDdLab[i] == "2k":
        t2k = t2k + 1
    if pepDdLab[i] == "2l":
        t2l = t2l + 1
    if pepDdLab[i] == "2m":
        t2m = t2m + 1
    if pepDdLab[i] == "2x":
        t2x = t2x + 1
print ("2a 2b 2c 2d 2e 2f 2g 2h 2i 2j 2k 2l 2m 2x")
print (t2a, t2b, t2c, t2d, t2e, t2f, t2g, t2h, t2i, t2j, t2k, t2l, t2m, t2x)
print (t2a + t2b + t2c + t2d + t2e + t2f + t2g + t2h + t2i + t2j + t2k + t2l + t2m + t2x)
