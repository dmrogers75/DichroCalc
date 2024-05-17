#!/usr/bin/python
import sys
import numpy as np

vector_1 = []
vector_2 = []

var1 = sys.argv[1] #("File name vector 1 (x y z): ")
file1 = open(var1,"r")
for line in file1:
   line = line.rstrip('\n')
   line = line.split( ) # split using spaces
   vector_1.append(float(line[0]))
   vector_1.append(float(line[1]))
   vector_1.append(float(line[2]))

var2 = sys.argv[2] #("File name vector 2 (x y z): ")
file2 = open(var2,"r")
for line in file2:
   line = line.rstrip('\n')
   line = line.split( ) # split using spaces
   vector_2.append(float(line[0]))
   vector_2.append(float(line[1]))
   vector_2.append(float(line[2]))

unit_vector_1 = vector_1 / np. linalg. norm(vector_1)
unit_vector_2 = vector_2 / np. linalg. norm(vector_2)
dot_product = np. dot(unit_vector_1, unit_vector_2)
angle = np. arccos(dot_product)
angled = np.degrees(angle)

print("\nDot product",dot_product) 
print("Angle in rad and deg.",angle,angled)
