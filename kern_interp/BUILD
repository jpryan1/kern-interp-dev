package(default_visibility = ["//visibility:public"])

load("//kern_interp:buildvars.bzl", "COPTS")

cc_library(
	name = "linear_solve",
	hdrs = ["linear_solve.h"],
  copts= COPTS,
	srcs = ["linear_solve.cpp"],
	deps = [":ki_mat",
			"//kern_interp/skel_factorization",
			"//kern_interp/quadtree",
			"//kern_interp/boundaries"],
  linkopts = ["-pthread"],

)


cc_library(
  name = "ki_mat",
  hdrs = ["ki_mat.h"],
  copts= COPTS,
  srcs = ["ki_mat.cpp",
    ],
  deps = [":pointvec"],
)



cc_library(
  name = "legendre",
  hdrs = ["legendre.h"],
    copts= COPTS,

  srcs = ["legendre.cpp",
    ],
)


cc_library(
	name = "pointvec",
	hdrs = ["pointvec.h"],
    copts= COPTS,

	srcs = ["pointvec.cpp"],
)


cc_test(
  name = "kern_interp_test",
  srcs = ["kern_interp_test.cpp"],
  copts= COPTS,
  data=["//kern_interp/boundaries/meshes/spheres:sphere_240_faces.txt","//kern_interp/boundaries/meshes/spheres:sphere_570_faces.txt",
  "//kern_interp/boundaries/meshes/spheres:sphere_1032_faces.txt",
  "//kern_interp/boundaries/meshes/spheres:sphere_2344_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_3206_faces.txt", 
  "//kern_interp/boundaries/meshes/spheres:sphere_4192_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_5322_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_6544_faces.txt", 
  "//kern_interp/boundaries/meshes/spheres:sphere_7970_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_9496_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_11146_faces.txt", 
  "//kern_interp/boundaries/meshes/spheres:sphere_12944_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_14870_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_16912_faces.txt",
   "//kern_interp/boundaries/meshes/spheres:sphere_19126_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_21424_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_23898_faces.txt",
    "//kern_interp/boundaries/meshes/spheres:sphere_26504_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_29214_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_32056_faces.txt",
     "//kern_interp/boundaries/meshes/spheres:sphere_35074_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_38184_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_41458_faces.txt",
      "//kern_interp/boundaries/meshes/spheres:sphere_44832_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_48370_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_52016_faces.txt", 
      "//kern_interp/boundaries/meshes/spheres:sphere_55802_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_59744_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_63810_faces.txt", 
      "//kern_interp/boundaries/meshes/spheres:sphere_67976_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_72322_faces.txt", "//kern_interp/boundaries/meshes/spheres:sphere_76736_faces.txt", 
      "//kern_interp/boundaries/meshes/spheres:sphere_81370_faces.txt"],
	deps = [
			":ki_mat",
			"//kern_interp/skel_factorization",
			"//kern_interp/quadtree",
			"//kern_interp/boundaries",
			":linear_solve",
			"@gtest//:gtest",
      "@gtest//:gtest_main"
    ],
  	linkopts = ["-llapacke",
  	        "-lopenblas",
			"-llapack",
			"-lblas",
			"-lgomp"]
)