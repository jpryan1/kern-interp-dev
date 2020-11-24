import numpy as np
import sys
from copy import copy
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt


# Ex 1 Concentric Spheres
is_channel_plot = False
ARROW_LENGTH = 0.4
BORDER_WIDTH = 10
HEAD_WIDTH = 3
QUIVER_RES_X = 5
QUIVER_RES_Y = 5
BOUNDARY_RES = 5
ZOOM = 1
TICK_LABEL_SIZE = 20
TICKS = [0,0.5,1]
VMIN = 0
VMAX = 1.1
BORDER_DIST = 0.2
OUTPUT_FILE = "up_cow_cross.eps"


# Ex 2 Pipe with ellipsoid holes 
# is_channel_plot = True
# ARROW_LENGTH = 0.1
# BORDER_WIDTH = 8
# HEAD_WIDTH = 2
# QUIVER_RES_X = 1
# QUIVER_RES_Y = 1
# BOUNDARY_RES = 5
# ZOOM = 1
# TICK_LABEL_SIZE = 40
# BORDER_DIST = 1.5
# OUTPUT_FILE = "pipe_crosssection.eps"

# Ex 3 Cow
# is_channel_plot = False
# ARROW_LENGTH = 0.5
# BORDER_WIDTH = 8
# HEAD_WIDTH = 3
# QUIVER_RES_X = 2
# QUIVER_RES_Y = 2
# BOUNDARY_RES = 1
# ZOOM = 1
# TICK_LABEL_SIZE = 40
# BORDER_DIST = 0.5
# OUTPUT_FILE = "cow_crosssection.eps"

print("args: {ZOOM} ")
fig, ax = plt.subplots(figsize=(8,8))
ax.axis("off")

MASKED_VALUE = 11111.1
if(len(sys.argv) > 1):
  ZOOM = float(sys.argv[1])


solution_lines = open("output/data/up_cow_sol.txt", "r").readlines()
# boundary_lines = open("kern_interp/boundaries/meshes/cow_mesh_cross.txt", "r").readlines()
boundary_lines = open("output/data/up_cow_bound.txt", "r").readlines()

is_stokes = (len(solution_lines[0].split(",")) == 4)
if is_stokes:
	CMAP = copy(matplotlib.cm.viridis)
else:
	CMAP = copy(matplotlib.cm.hot)
CMAP.set_bad("white",1.)
solution_dim = int(np.sqrt(len(solution_lines)))
solution_grid = np.array([[MASKED_VALUE for x in range(solution_dim)] for y in range(solution_dim)])

X, Y, U, V = [], [], [], []
min_sol_x, min_sol_y, max_sol_x, max_sol_y = 10,10,-10,-10

min_val_seen = 1000
max_val_seen = -1000

for i in range(solution_dim):
	for j in range(solution_dim):
		linesplit = [float(n) for n in solution_lines[i+solution_dim*j].split(',')]
		min_sol_x = min(min_sol_x, (linesplit[0]))
		max_sol_x = max(max_sol_x, (linesplit[0]))
		min_sol_y = min(min_sol_y, (linesplit[1]))
		max_sol_y = max(max_sol_y, (linesplit[1]))

		mag = np.sqrt((linesplit[2])**2 + (linesplit[3])**2) if is_stokes \
			else linesplit[2]
		if(mag > 10): mag = 0
		solution_grid[i][j] = mag if mag!=0 else MASKED_VALUE
		if mag !=0:
			min_val_seen = min(min_val_seen, mag)
			max_val_seen = max(max_val_seen, mag)
		if(is_stokes):
			if(i % QUIVER_RES_X != 0 or j % QUIVER_RES_Y != 0 \
				or np.sqrt((linesplit[2])**2 + (linesplit[3])**2)<0.01):
				continue
			X.append((linesplit[0]))
			Y.append((linesplit[1]))
			U.append((linesplit[2]))
			V.append((linesplit[3]))

solution_grid = np.ma.masked_where(solution_grid == MASKED_VALUE, solution_grid)


imsh = ax.imshow(solution_grid,
	extent=[min_sol_x, max_sol_x, min_sol_y, max_sol_y], origin="lower", \
	cmap=CMAP, interpolation="bilinear", vmin=VMIN, vmax=VMAX)
if(is_stokes):
	quiver_scale = (10 / ARROW_LENGTH ) * ZOOM
	ax.quiver(X,Y,U,V, color="white",scale=quiver_scale, headwidth=HEAD_WIDTH)

# Boundary plot

hole_start_idx = 0
for i in range(0,len(boundary_lines)-BOUNDARY_RES,BOUNDARY_RES):
	pixel = boundary_lines[i].split(",")
	pixel = [float(pixel[0]), float(pixel[1])]
	next_pixel = boundary_lines[i+BOUNDARY_RES].split(",")
	next_pixel = [float(next_pixel[0]), float(next_pixel[1])]
	if((pixel[0]-next_pixel[0])**2+(pixel[1]-next_pixel[1])**2>BORDER_DIST**2):
		next_pixel = boundary_lines[hole_start_idx].split(",")
		next_pixel = [float(next_pixel[0]), float(next_pixel[1])]
		ax.plot([pixel[0], next_pixel[0]], [pixel[1],next_pixel[1]], \
		linewidth=BORDER_WIDTH, color="black")
		hole_start_idx = i+BOUNDARY_RES
		continue
	ax.plot([pixel[0], next_pixel[0]], [pixel[1],next_pixel[1]], \
		linewidth=BORDER_WIDTH, color="black")
	
#Now the buffer
# hole_start_idx = 0
# for i in range(0,len(boundary_lines)-BOUNDARY_RES,BOUNDARY_RES):
# 	pixel = boundary_lines[i].split(",")
# 	pixel = [float(pixel[0]), float(pixel[1])]
# 	next_pixel = boundary_lines[i+BOUNDARY_RES].split(",")
# 	next_pixel = [float(next_pixel[0]), float(next_pixel[1])]
# 	if((pixel[0]-next_pixel[0])**2+(pixel[1]-next_pixel[1])**2>BORDER_DIST**2):
# 		next_pixel = boundary_lines[hole_start_idx].split(",")
# 		next_pixel = [float(next_pixel[0]), float(next_pixel[1])]
# 		ax.plot([0.9*pixel[0], 0.9*next_pixel[0]], [0.9*pixel[1],0.9*next_pixel[1]], \
# 			linewidth=BUFFER_WIDTH, color="white")
# 		ax.plot([1.1*pixel[0], 1.1*next_pixel[0]], [1.1*pixel[1],1.1*next_pixel[1]], \
# 			linewidth=BUFFER_WIDTH, color="white")
# 		hole_start_idx = i+BOUNDARY_RES
# 		continue
# 	ax.plot([0.9*pixel[0], 0.9*next_pixel[0]], [0.9*pixel[1],0.9*next_pixel[1]], \
# 		linewidth=BUFFER_WIDTH, color="white")
# 	ax.plot([1.1*pixel[0], 1.1*next_pixel[0]], [1.1*pixel[1],1.1*next_pixel[1]], \
# 		linewidth=BUFFER_WIDTH, color="white")

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.5)
cbar= plt.colorbar(imsh, cax=cax, ticks=TICKS)
cbar.ax.tick_params(labelsize=TICK_LABEL_SIZE)


xl, xr = ax.get_xlim()
yl, yr = ax.get_ylim()
l = min(xl,yl)-0.01
r = max(xr,yr)+0.01
# ax.set_xlim((l - (r+l)/2.)/ZOOM + (r+l)/2., (r - (r+l)/2.)/ZOOM + (r+l)/2.)
# ax.set_ylim((l - (r+l)/2.)/ZOOM + (r+l)/2., (r - (r+l)/2.)/ZOOM + (r+l)/2.)
if(OUTPUT_FILE == "cow_crosssection.eps"):
	ax.set_xlim((-0.6,0.6))
	ax.set_ylim((-0.6,0.6))
if(is_channel_plot):
	if OUTPUT_FILE == "ex1_new.eps":
		xl-=0.2
		xr+=0.2
		ax.set_xlim((xl - (xr+xl)/2.)/ZOOM + (xr+xl)/2., (xr - (xr+xl)/2.)/ZOOM + (xr+xl)/2.)
	else:
		yl-=0.2
		yr+=0.2
		ax.set_ylim((yl - (yr+yl)/2.)/ZOOM + (yr+yl)/2., (yr - (yr+yl)/2.)/ZOOM + (yr+yl)/2.)

plt.savefig(OUTPUT_FILE, bbox_inches="tight", format="eps")
plt.show()
