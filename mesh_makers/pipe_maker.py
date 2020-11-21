#!/usr/bin/env python
from __future__ import division

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.cm as cm


def sumtriangles( xyz, triangles ):
    areasum = 0
    smallest = 10
    largest = 0
    for tri in triangles:
        corners = xyz[tri]
        a = np.arccos(np.dot(corners[0], corners[1]))
        b = np.arccos(np.dot(corners[0], corners[2]))
        c = np.arccos(np.dot(corners[1], corners[2]))
        s=(a+b+c)/2.0
        lhuilier = np.sqrt(np.tan(0.5*(s))*np.tan(0.5*(s-a))*np.tan(0.5*(s-b))*np.tan(0.5*(s-c)))
        area = 4*np.arctan(lhuilier)
        if(area < smallest):
            smallest = area
        if(area>largest):
            largest = area
        areasum += area
    print("Smallest: ", smallest," largest " ,largest)
    return areasum

def bumpfun(u, N):
    #u from 0 to 1 goes from 1 to N+1 to 1
    parab = -u*(u-1)*N*4
    return int(1+parab)

#...............................................................................
if __name__ == "__main__":
    import sys
    from scipy.spatial import ConvexHull
    # fig, ax = plt.subplots(figsize=(14,14))
    # ax = fig.gca(projection='3d')
    N=40
    points = []
    cyl_height = 5
    height_disc = 100
    for j in range(int(N/2)):
        jflt = float(j)/(N-1)
        sz = bumpfun(jflt, 2*N) #1 to 41 to 1
        phi = jflt*np.pi
        for i in range(sz):
            iflt = float(i)/(sz)
            
            theta = (iflt)*2.0*np.pi
        
            # phi = np.arccos(1-2*jflt)
            points.append(  [np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), cyl_height+np.cos(phi)])


    for height in range(height_disc):
        #param 0 to 1 in range
        height_param = 1-(height/float(height_disc-1))
        bulge_param = 1 + 0.001*(np.sin((np.pi/4.)+height_param*(np.pi/2.0))-np.sqrt(2)/2.)
        jflt = int(N/2)/float(N-1)
        sz = bumpfun(jflt, 2*N)
        for i in range(sz):
            iflt = float(i)/(sz)
            theta = (iflt)*2.0*np.pi
            phi = jflt*np.pi

            # phi = np.arccos(1-2*jflt)
            points.append(  [bulge_param*np.cos(theta)*np.sin(phi), 
                bulge_param*np.sin(theta)*np.sin(phi), 
                (cyl_height*height_param)])


    for j in range(int(N/2)+1, N):
        jflt = float(j)/(N-1)
        sz = bumpfun(jflt, 2*N) #1 to 41 to 1
        phi = jflt*np.pi
        for i in range(sz):
            iflt = float(i)/(sz)
            
            theta = (iflt)*2.0*np.pi

            # phi = np.arccos(1-2*jflt)
            points.append(  [np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), np.cos(phi)])

    points = np.array(points)
    tessellation = ConvexHull( points )
    tris = tessellation.simplices  # ntri, dim+1
    num_verts = len(points)
    print(num_verts)
    num_faces = len(tris)
    printout = """ply
format ascii 1.0
element vertex """+str(num_verts)+"""
property float x
property float y
property float z
element face """+str(num_faces)+"""
property list uchar int vertex_indices
end_header\n"""
    for point in points:
        # outer prod should dot with pt to be pos
        printout+=(str(point[0])+" "+str(point[1])+" "+str(point[2]))+"\n"
    for tri in tris:
        corners = points[tri]
        v1 = corners[1]-corners[0]
        v2 = corners[2]-corners[0]
        crs = np.cross(v1, v2)
        if(np.dot(crs, corners[0])>0):
            printout+=("3 "+str(tri[0])+" "+str(tri[1])+" "+str(tri[2]))+"\n"
        else:
            printout+=("3 "+str(tri[0])+" "+str(tri[2])+" "+str(tri[1]))+"\n"

        # exit()
        # corners = points[tri]
        # a = np.arccos(np.dot(corners[0], corners[1]))
        # b = np.arccos(np.dot(corners[0], corners[2]))
        # c = np.arccos(np.dot(corners[1], corners[2]))
        # s=(a+b+c)/2.0
        # lhuilier = np.sqrt(np.tan(0.5*(s))*np.tan(0.5*(s-a))*np.tan(0.5*(s-b))*np.tan(0.5*(s-c)))
        # area = 4*np.arctan(lhuilier)
        # mid = corners[0]+corners[1]+corners[2]
        # mid = mid/3.0
        # mid = mid/np.linalg.norm(mid)
        # printout += str(mid[0])+","+str(mid[1])+","+str(mid[2])+","+str(area)+"\n"
    outp = open("ply_files/pipe_mesh.ply", "w")
    outp.write(printout)
    outp.close()
# print("N", N, " num triangles ",len(tessellation.simplices))