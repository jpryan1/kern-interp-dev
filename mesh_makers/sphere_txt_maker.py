#!/usr/bin/env python
from __future__ import division
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

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
    for N in range(10, 30, 5):
    # N=150
        points = []
        for j in range(N):
        	jflt = float(j)/(N-1)
        	sz = bumpfun(jflt, 2*N) #1 to 41 to 1
        	for i in range(sz):
        		iflt = float(i)/(sz)
        		
        		theta = (iflt)*2.0*np.pi
        		phi = jflt*np.pi
        		# phi = np.arccos(1-2*jflt)
        		points.append(  [np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), np.cos(phi)])

        points = np.array(points)
        tessellation = ConvexHull( points )
        tris = tessellation.simplices  # ntri, dim+1
        printout = ""
        for tri in tris:
            corners = points[tri]
            a = np.arccos(np.dot(corners[0], corners[1]))
            b = np.arccos(np.dot(corners[0], corners[2]))
            c = np.arccos(np.dot(corners[1], corners[2]))
            s=(a+b+c)/2.0
            lhuilier = np.sqrt(np.tan(0.5*(s))*np.tan(0.5*(s-a))*np.tan(0.5*(s-b))*np.tan(0.5*(s-c)))
            area = 4*np.arctan(lhuilier)
            mid = corners[0]+corners[1]+corners[2]
            mid = mid/3.0
            mid = mid/np.linalg.norm(mid)
            printout += str(mid[0])+","+str(mid[1])+","+str(mid[2])+","+str(area)+"\n"
        outp = open("sphere_"+str(len(tris))+"_faces.txt", "w")
        outp.write(printout)
        outp.close()
        areasum = sumtriangles( points, tris)
        print(areasum-4*np.pi)
        # print("N", N, " num triangles ",len(tessellation.simplices))
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')

        # ax.plot_trisurf(points[:,0], points[:,1], triangles=tessellation.simplices, Z=points[:,2], shade=True, cmap=matplotlib.colors.ListedColormap ( np.random.rand ( 256,3)))
        # plt.show()
