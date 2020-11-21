#!/usr/bin/env python
from __future__ import division
import numpy as np


FILENAME = "ply_files/cow_mesh"
ply_lines = open(FILENAME+".ply", "r").readlines()

line_idx = 0
num_verts = 0
num_faces = 0
while(True):
    line = ply_lines[line_idx]
    if "vertex" in line:
        num_verts = int(line.split(" ")[2])
    if "face" in line:
        num_faces = int(line.split(" ")[2])
        line_idx+=3
        break
    line_idx += 1

verts_begin = line_idx
points =[]
max_norm = 0
for i in range(num_verts):
    line = ply_lines[line_idx]
    line_idx+=1
    ls = line.split(" ")
    point = [float(ls[i]) for i in range(3)]
    max_norm=max(max_norm, np.linalg.norm(np.array(point)))

print(max_norm)
line_idx = verts_begin
for i in range(num_verts):
    line = ply_lines[line_idx]
    line_idx+=1
    ls = line.split(" ")
    point = [float(ls[i]) for i in range(3)]
    point[0] *= 1.0/max_norm
    point[1] *= 1.0/max_norm
    point[2] *= 1.0/max_norm
    points.append(point)

points = np.array(points)

printout = ""
line_idx = verts_begin
for i in range(verts_begin):
    printout += ply_lines[i]
for i in range(num_verts):
    line = ply_lines[line_idx]
    line_idx+=1
    ls = line.split(" ")
    printout += str(points[i,0]) + " " + str(points[i,1]) + " " + str(points[i,2]) + " " + ls[3]+" "+ls[4]
for i in range(num_faces):
    printout += ply_lines[line_idx]
    line_idx+=1
  
outp = open(FILENAME+"_scaled.ply", "w")
outp.write(printout)
outp.close()
