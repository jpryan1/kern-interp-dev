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
  data=["//kern_interp/boundaries/sphere_triangulations:tri240.txt","//kern_interp/boundaries/sphere_triangulations:tri570.txt","//kern_interp/boundaries/sphere_triangulations:tri1032.txt","//kern_interp/boundaries/sphere_triangulations:tri1614.txt","//kern_interp/boundaries/sphere_triangulations:tri2344.txt",
  "//kern_interp/boundaries/sphere_triangulations:tri3206.txt","//kern_interp/boundaries/sphere_triangulations:tri4192.txt","//kern_interp/boundaries/sphere_triangulations:tri5322.txt","//kern_interp/boundaries/sphere_triangulations:tri6544.txt","//kern_interp/boundaries/sphere_triangulations:tri7970.txt","//kern_interp/boundaries/sphere_triangulations:tri9496.txt",
  "//kern_interp/boundaries/sphere_triangulations:tri11146.txt","//kern_interp/boundaries/sphere_triangulations:tri12944.txt","//kern_interp/boundaries/sphere_triangulations:tri14870.txt","//kern_interp/boundaries/sphere_triangulations:tri16912.txt","//kern_interp/boundaries/sphere_triangulations:tri19126.txt","//kern_interp/boundaries/sphere_triangulations:tri21424.txt",
  "//kern_interp/boundaries/sphere_triangulations:tri23898.txt","//kern_interp/boundaries/sphere_triangulations:tri26504.txt","//kern_interp/boundaries/sphere_triangulations:tri29214.txt","//kern_interp/boundaries/sphere_triangulations:tri32056.txt","//kern_interp/boundaries/sphere_triangulations:tri35074.txt","//kern_interp/boundaries/sphere_triangulations:tri38184.txt",
  "//kern_interp/boundaries/sphere_triangulations:tri41458.txt","//kern_interp/boundaries/sphere_triangulations:tri44832.txt","//kern_interp/boundaries/sphere_triangulations:tri48370.txt","//kern_interp/boundaries/sphere_triangulations:tri52016.txt","//kern_interp/boundaries/sphere_triangulations:tri55802.txt","//kern_interp/boundaries/sphere_triangulations:tri59744.txt",
  "//kern_interp/boundaries/sphere_triangulations:tri63810.txt","//kern_interp/boundaries/sphere_triangulations:tri67976.txt","//kern_interp/boundaries/sphere_triangulations:tri72322.txt","//kern_interp/boundaries/sphere_triangulations:tri76736.txt","//kern_interp/boundaries/sphere_triangulations:tri81370.txt"],
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