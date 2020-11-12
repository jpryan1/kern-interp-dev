#!/usr/bin/env python
from __future__ import division
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.cm as cm

FILENAME="cow_mesh"
ply_lines = open("ply_files/"+FILENAME+".ply", "r").readlines()

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


points =[]
for i in range(num_verts):
    line = ply_lines[line_idx]
    line_idx+=1
    ls = line.split(" ")
    point = [float(ls[i]) for i in range(3)]
    points.append(point)

points = np.array(points)

printout = ""
for i in range(num_faces):
    line = ply_lines[line_idx]
    line_idx+=1
    ls = line.split(" ")
    tri = [int(ls[i]) for i in range(1,4)]

    corners = points[tri]
    a = np.arccos(np.dot(corners[0], corners[1])/(np.linalg.norm(corners[0])*np.linalg.norm(corners[1])))
    b = np.arccos(np.dot(corners[0], corners[2])/(np.linalg.norm(corners[0])*np.linalg.norm(corners[2])))
    c = np.arccos(np.dot(corners[1], corners[2])/(np.linalg.norm(corners[1])*np.linalg.norm(corners[2])))
    s = (a+b+c)/2
    # if(s*(s-a)*(s-b)*(s-c) < 0):
    #     print("inside ", s*(s-a)*(s-b)*(s-c), " points ", corners, " s ", s, " ", a,b,c)
    area = np.sqrt((s*(s-a)*(s-b)*(s-c)))
    mid = corners[0]+corners[1]+corners[2]
    mid = mid/3.0
    mid = mid/np.linalg.norm(mid)
    printout += str(mid[0])+","+str(mid[1])+","+str(mid[2])+","+str(area)+"\n"

outp = open("kern-interp_inputs/"+FILENAME+".txt", "w")
outp.write(printout)
outp.close()
