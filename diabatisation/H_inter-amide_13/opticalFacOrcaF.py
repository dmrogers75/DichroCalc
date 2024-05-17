#!/usr/bin/python
# Compute Schellman's optical factor [R_21 . mu_2 x mu_1]
import sys
import numpy as np

vector_1 = [] # mu_1
vector_2 = [] # mu_2
vector_3 = [] # R_1
vector_4 = [] # R_2

var1 = sys.argv[1] # mu_1
file1 = open(var1,"r")
for line in file1:
   line = line.rstrip('\n')
   line = line.split( ) # split using spaces
   vector_1.append(float(line[0]))
   vector_1.append(float(line[1]))
   vector_1.append(float(line[2]))
#print(vector_1)

var2 = sys.argv[2] # mu_2
file2 = open(var2,"r")
for line in file2:
   line = line.rstrip('\n')
   line = line.split( ) # split using spaces
   vector_2.append(float(line[0]))
   vector_2.append(float(line[1]))
   vector_2.append(float(line[2]))
#print(vector_2)

var3 = sys.argv[3] # R_1
file3 = open(var3,"r")
for line in file3:
   line = line.rstrip('\n')
   line = line.split( ) # split using spaces
   vector_3.append(float(line[0]))
   vector_3.append(float(line[1]))
   vector_3.append(float(line[2]))
#print(vector_3)

var4 = sys.argv[4] # R_2
file4 = open(var4,"r")
for line in file4:
   line = line.rstrip('\n')
   line = line.split( ) # split using spaces
   vector_4.append(float(line[0]))
   vector_4.append(float(line[1]))
   vector_4.append(float(line[2]))
#print(vector_4)

# Convert to atomic units if mu in debye and R in angstrom
icon = 0 # default
if len(sys.argv) > 4:
    iconInp = sys.argv[5]
icon = int(iconInp)
if icon == 0: # no conversion
    debyeAu = 1.0
    angstromBohr = 1.0
if icon == 1: # if mu in au and R in angstrom
    debyeAu = 1.0
    angstromBohr = 1.88973
if icon == 2: # if mu in debye and R in angstrom
    debyeAu = 0.393422
    angstromBohr = 1.88973

Auvector_1 = np.multiply(vector_1, debyeAu)
Auvector_2 = np.multiply(vector_2, debyeAu)
Auvector_3 = np.multiply(vector_3, angstromBohr)
Auvector_4 = np.multiply(vector_4, angstromBohr)

Aucross_product_21 = np.cross(Auvector_2, Auvector_1) # un-normalised 
Auseparation_43 = np.subtract(Auvector_4, Auvector_3) # un-normalised
Audot_product_s43cp21 = np.dot(Auseparation_43, Aucross_product_21) # un-normalised
print("\nUnits: unit conversion option icon = ",icon)
print("\nmu_2 x mu_1 (au) =",Aucross_product_21)
print("R_21 (bohr) =",Auseparation_43)
print("[R_21 . mu_2 x mu_1] (au) =",Audot_product_s43cp21,"\n")

Audot_product_12 = np.dot(Auvector_1, Auvector_2)
Audot_product_1s43 = np.dot(Auvector_1, Auseparation_43)
Audot_product_2s43 = np.dot(Auvector_2, Auseparation_43)
Aunorm_R_21 = np.linalg.norm(np.subtract(Auvector_4, Auvector_3))
if Aunorm_R_21 != 0.0:
    Vdd_au = Audot_product_12 / (Aunorm_R_21)**3
    Vdd_au = Vdd_au - ( (3.0 * Audot_product_1s43 * Audot_product_2s43) / (Aunorm_R_21)**5 )
    print("V_dd (au) =", Vdd_au,"\n")
    print("V_dd (eV) =", Vdd_au * 27.2114,"\n")
    print("V_dd (cm-1) =", Vdd_au * 219474.63,"\n")
