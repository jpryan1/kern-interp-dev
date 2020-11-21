package(default_visibility = ["//visibility:public"])
load("//kern_interp:buildvars.bzl", "COPTS")

cc_library(
	name = "boundaries",
	hdrs = [
			"circle.h",
			"concentric_spheres.h",
			"pipe_mesh.h",
			"cow_mesh.h",
			"annulus.h",
			"uniformdist.h",
			"donut.h",
			"cubic_spline.h",
			"ex1boundary.h",
			"spiral_channel.h",
			"ex2boundary.h",
			"ex3boundary.h",
			"cross_section_border.h",
			"boundary.h"],
	copts= COPTS,
	srcs = [
			"circle.cpp",
			"concentric_spheres.cpp",
			"pipe_mesh.cpp",
			"cow_mesh.cpp",
			"annulus.cpp",
			"uniformdist.cpp",
			"donut.cpp",
			"cubic_spline.cpp",
			"boundary.cpp",
			"cubic_boundary.cpp",
			"spiral_channel.cpp",
			"ex1boundary.cpp",
			"ex2boundary.cpp",
			"ex3boundary.cpp",
			"cross_section_border.cpp",
			],
	deps = ["//kern_interp:ki_mat",
			"//kern_interp:pointvec",
			"//kern_interp:legendre",
			]
)
