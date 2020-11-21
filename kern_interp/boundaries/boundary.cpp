// Copyright 2019 John Paul Ryan
#include <cmath>
#include <cassert>
#include <iostream>
#include "kern_interp/boundaries/boundary.h"


namespace kern_interp {


bool Boundary::comp_ang(PointVec a, PointVec b) {
  return atan2(a.a[1], a.a[0]) < atan2(b.a[1], b.a[0]);
}


void Boundary::set_boundary_values_size(BoundaryCondition bc) {
  int num_points = weights.size();

  switch (bc) {
    case SINGLE_ELECTRON:
    case EX3A:
    case ALL_NEG_ONES:
    case ALL_ZEROS:
    case ALL_ONES:
    case LAPLACE_CHECK_2D:
    case LAPLACE_CHECK_3D:

    case ELECTRON_3D:
    case LAPLACE_PIPE_HOLES:
    case LAPLACE_COW_MESH:
      boundary_values = ki_Mat(num_points, 1);
      break;
    // Everything above is for 1D solutions.
    // Everything below is for 2D solutions.
    case TANGENT_VEC:
    case REVERSE_TANGENT_VEC:
    case NORMAL_VEC:
    case REVERSE_NORMAL_VEC:
    case EX1:
    case LEFT_TO_RIGHT_FLOW:
    case EX3B:
    case HORIZONTAL_VEC:
    case NO_SLIP:
    case STOKES_2D_MIX:
      boundary_values = ki_Mat(2 * num_points, 1);
      break;
    case STOKES_SPHERES:
    case STOKES_PIPE_HOLES:
    case STOKES_COW_MESH:
      boundary_values = ki_Mat(3 * num_points, 1);
      break;
    case DEFAULT: {
      std::cout << "DEFAULT BoundaryCondition enum not to be given " <<
                "directly to set_boundary_value_size function." << std::endl;
      exit(-1);
      break;
    }
  }
}

// TODO(John) this shouldn't be that hard to make dimension agnostic
void Boundary::apply_boundary_condition(int start_point_idx, int end_point_idx,
                                        BoundaryCondition bc) {
  for (int point_idx = start_point_idx;
       point_idx < end_point_idx;
       point_idx++) {
    switch (bc) {
      // 1D solutions
      case SINGLE_ELECTRON: {
        boundary_values.set(point_idx, 0,
                            log(sqrt(pow(points[2 * point_idx] + 3, 2)
                                     + pow(points[2 * point_idx + 1] + 2, 2)))
                            / (2 * M_PI));
        break;
      }
      case LAPLACE_CHECK_2D: {

        double r = sqrt(pow(points[2 * point_idx] - 0.0, 2)
                        + pow(points[2 * point_idx + 1] - 0.0, 2));
        if (r < 0.6) {
          boundary_values.set(point_idx, 0, 1.);
        } else {
          boundary_values.set(point_idx, 0, 3.);
        }
        break;
      }
      case ALL_ONES: {
        boundary_values.set(point_idx, 0, 1.);
        break;
      }
      case ALL_NEG_ONES: {
        boundary_values.set(point_idx, 0, -1.);
        break;
      }
      case ALL_ZEROS: {
        boundary_values.set(point_idx, 0, 0.);
        break;
      }
      // 2D solutions.
      case TANGENT_VEC: {
        boundary_values.set(2 * point_idx, 0, -normals[2 * point_idx + 1]);
        boundary_values.set(2 * point_idx + 1, 0, normals[2 * point_idx]);
        break;
      }
      case REVERSE_TANGENT_VEC: {
        boundary_values.set(2 * point_idx, 0, normals[2 * point_idx + 1]);
        boundary_values.set(2 * point_idx + 1, 0, -normals[2 * point_idx]);
        break;
      }
      case NORMAL_VEC: {
        boundary_values.set(2 * point_idx, 0, normals[2 * point_idx ]);
        boundary_values.set(2 * point_idx + 1, 0, normals[2 * point_idx + 1]);
        break;
      }
      case REVERSE_NORMAL_VEC: {
        boundary_values.set(2 * point_idx, 0, -normals[2 * point_idx ]);
        boundary_values.set(2 * point_idx + 1, 0, -normals[2 * point_idx + 1]);
        break;
      } case HORIZONTAL_VEC: {
        boundary_values.set(2 * point_idx, 0, 1.);
        boundary_values.set(2 * point_idx + 1, 0, 0.);
        break;
      }
      case LEFT_TO_RIGHT_FLOW: {
        double x = points[2 * point_idx];
        if (x < -0.99 || x > 1.99) {
          boundary_values.set(2 * point_idx, 0, 1);
          boundary_values.set(2 * point_idx + 1, 0, 0);
        } else {
          boundary_values.set(2 * point_idx, 0, 0);
          boundary_values.set(2 * point_idx + 1, 0, 0);
        }
        break;
      }
      case EX1: {
        double x = points[2 * point_idx ];
        double y = points[2 * point_idx + 1];
        if ((y < 0.01 || y > 2.99) && x > 0.01 && x < 0.49) {
          // 0 to 0.5
          double phase = x * (4 * M_PI);
          double mag = 0.5 - 0.5 * cos(phase);
          boundary_values.set(2 * point_idx, 0, 0);
          boundary_values.set(2 * point_idx + 1, 0, -2 * mag);
        } else if (x < -0.1 && (
                     (std::abs(y - 2.75) < 0.01) ||
                     (std::abs(y - 0.5) < 0.01)
                   )) {
          // -0.5 to -0.25
          double phase = (x + 0.5) * (8 * M_PI);
          double mag = 0.5 - 0.5 * cos(phase);
          boundary_values.set(2 * point_idx, 0, 0);
          boundary_values.set(2 * point_idx + 1, 0, -2 * mag);
        } else if (x > 0.76 && (
                     (std::abs(y - 3.0) < 0.01) ||
                     (std::abs(y - 0.25) < 0.01)
                   )) {
          // 0.75 to 1
          double phase = (x - 0.75) * (8 * M_PI);
          double mag = 0.5 - 0.5 * cos(phase);
          boundary_values.set(2 * point_idx, 0, 0);
          boundary_values.set(2 * point_idx + 1, 0, -2 * mag);
        } else if ((x < -0.26 && (std::abs(y - 1.25) < 0.01))) {
          // -0.5 to -0.25
          double phase = (x + 0.5) * (8 * M_PI);
          double mag = 0.5 - 0.5 * cos(phase);
          boundary_values.set(2 * point_idx, 0, 0);
          boundary_values.set(2 * point_idx + 1, 0, -4 * mag);
        } else if ((x > 0.76 && (std::abs(y - 1.75) < 0.01))) {
          // 0.75 to 1
          double phase = (x - 0.75) * (8 * M_PI);
          double mag = 0.5 - 0.5 * cos(phase);
          boundary_values.set(2 * point_idx, 0, 0);
          boundary_values.set(2 * point_idx + 1, 0, -4 * mag);
        } else {
          boundary_values.set(2 * point_idx, 0, 0);
          boundary_values.set(2 * point_idx + 1, 0, 0);
        }
        break;
      }
      case NO_SLIP: {
        boundary_values.set(2 * point_idx, 0, 0.);
        boundary_values.set(2 * point_idx + 1, 0, 0.);
        break;
      }
      case STOKES_2D_MIX: {
        double r = sqrt(pow(points[2 * point_idx] - 0.0, 2)
                        + pow(points[2 * point_idx + 1] - 0.0, 2)
                       );
        if (r > 0.9) {
          boundary_values.set(2 * point_idx, 0, 1.0);
          boundary_values.set(2 * point_idx + 1, 0, 0);
        } else {
          boundary_values.set(2 * point_idx, 0, -STOKES_MIXER);
          boundary_values.set(2 * point_idx + 1, 0, 0);
        }
        break;
      }
      // 3D domain
      case ELECTRON_3D: {
        double r = sqrt(pow(points[3 * point_idx], 2)
                        + pow(points[3 * point_idx + 1] , 2)
                        + pow(points[3 * point_idx + 2] + 1.5, 2));
        // double r2 = sqrt(pow(points[3 * point_idx] , 2)
        //                  + pow(points[3 * point_idx + 1] , 2)
        //                  + pow(points[3 * point_idx + 2] +1.1, 2));
        boundary_values.set(point_idx, 0,
                            (-1.0 / (4 * M_PI * r)));
        // - (-1.0 / (4 * M_PI * r2)));
        break;
      } case LAPLACE_PIPE_HOLES: {
        double r = sqrt(pow(points[3 * point_idx], 2)
                        + pow(points[3 * point_idx + 1] , 2)
                        + pow(points[3 * point_idx + 2] - 2, 2));
        double r2 = sqrt(pow(points[3 * point_idx] , 2)
                         + pow(points[3 * point_idx + 1] , 2)
                         + pow(points[3 * point_idx + 2] - 5, 2));
        double r3 = sqrt(pow(points[3 * point_idx] , 2)
                         + pow(points[3 * point_idx + 1] , 2)
                         + pow(points[3 * point_idx + 2] - 8, 2));
        boundary_values.set(point_idx, 0,
                            (-1.0 / (4 * M_PI * r))
                            + (-1.0 / (4 * M_PI * r2))
                            + (-1.0 / (4 * M_PI * r3)));
        break;
      } case LAPLACE_COW_MESH: {
        double r = sqrt(pow(points[3 * point_idx], 2)
                        + pow(points[3 * point_idx + 1] , 2)
                        + pow(points[3 * point_idx + 2] , 2));
        double r2 = sqrt(pow(points[3 * point_idx] , 2)
                         + pow(points[3 * point_idx + 1] - 4, 2)
                         + pow(points[3 * point_idx + 2] , 2));
        double r3 = sqrt(pow(points[3 * point_idx] , 2)
                         + pow(points[3 * point_idx + 1] - 8, 2)
                         + pow(points[3 * point_idx + 2] , 2));
        double r4 = sqrt(pow(points[3 * point_idx] , 2)
                         + pow(points[3 * point_idx + 1] + 5, 2)
                         + pow(points[3 * point_idx + 2], 2));
        boundary_values.set(point_idx, 0,
                            (-1.0 / (4 * M_PI * r))
                            - (-1.0 / (4 * M_PI * r2))
                            + (-1.0 / (4 * M_PI * r3))
                            - 5 * (-1.0 / (4 * M_PI * r4)));
        break;
      }
      case LAPLACE_CHECK_3D: {

        double r = sqrt(pow(points[3 * point_idx] - 0.0, 2)
                        + pow(points[3 * point_idx + 1] - 0.0, 2)
                        + pow(points[3 * point_idx + 2] - 0.0, 2));
        if (r < 0.2) {
          boundary_values.set(point_idx, 0, 1.);
        } else {
          boundary_values.set(point_idx, 0, 3.);
        }
        break;
      }
      case STOKES_SPHERES: {
        double r = sqrt(pow(points[3 * point_idx] - 0.0, 2)
                        + pow(points[3 * point_idx + 1] - 0.0, 2)
                        + pow(points[3 * point_idx + 2] - 0.0, 2));
        if (r > 0.9) {
          boundary_values.set(3 * point_idx, 0, 0);
          boundary_values.set(3 * point_idx + 1, 0, 0.0);
          boundary_values.set(3 * point_idx + 2, 0, 1.0);
        } else {
          boundary_values.set(3 * point_idx, 0, 0);
          boundary_values.set(3 * point_idx + 1, 0, 0);
          boundary_values.set(3 * point_idx + 2, 0, -STOKES_MIXER);
        }
        break;
      }
      case STOKES_PIPE_HOLES: {
        if (points[3 * point_idx + 2] < 0 || points[3 * point_idx + 2] > 5) {
          boundary_values.set(3 * point_idx, 0, 0);
          boundary_values.set(3 * point_idx + 1, 0, 0);
          boundary_values.set(3 * point_idx + 2, 0, 1);
        } else {
          boundary_values.set(3 * point_idx, 0, 0);
          boundary_values.set(3 * point_idx + 1, 0, 0);
          boundary_values.set(3 * point_idx + 2, 0, 0);
        }
        break;
      }
      case STOKES_COW_MESH: {
        double r = sqrt(pow(points[3 * point_idx] - 0.0, 2)
                        + pow(points[3 * point_idx + 1] - 0.0, 2)
                        + pow(points[3 * point_idx + 2] - 0.0, 2));
        if (r > 0.9) {
          boundary_values.set(3 * point_idx, 0, 0);
          boundary_values.set(3 * point_idx + 1, 0, 1.0);
          boundary_values.set(3 * point_idx + 2, 0, 0.0);
        } 
        break;
        // bool already_set = false;
        // for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {

        //   Hole hole = holes[hole_idx];
        //   double cx = hole.center.a[0];
        //   double cy = hole.center.a[1];
        //   double cz = hole.center.a[2];
        //   double r = sqrt(pow(points[3 * point_idx] - cx, 2)
        //                   + pow(points[3 * point_idx + 1] - cy, 2)
        //                   + pow(points[3 * point_idx + 2] - cz, 2));
        //   if (r < hole.radius) {
        //     if (hole_idx == 0) {
        //       boundary_values.set(3 * point_idx, 0, normals[3 * point_idx]);
        //       boundary_values.set(3 * point_idx + 1, 0, normals[3 * point_idx + 1]);
        //       boundary_values.set(3 * point_idx + 2, 0, normals[3 * point_idx + 2]);
        //       already_set = true;
        //       break;
        //     } else if (hole_idx == 1) {
        //       already_set = true;
        //       continue;
        //     } else {
        //       boundary_values.set(3 * point_idx, 0, -normals[3 * point_idx]);
        //       boundary_values.set(3 * point_idx + 1, 0, -normals[3 * point_idx + 1]);
        //       boundary_values.set(3 * point_idx + 2, 0, -normals[3 * point_idx + 2]);
        //       already_set = true;
        //       break;
        //     }
        //   }
        // }
        // if (!already_set) {
        //   boundary_values.set(3 * point_idx, 0, 0);
        //   boundary_values.set(3 * point_idx + 1, 0, 1);
        //   boundary_values.set(3 * point_idx + 2, 0, 0);
        // }
        // break;
      }

      case EX3A:
      case EX3B:
      case DEFAULT: {
        std::cout << "BoundaryCondition enum not to be given " <<
                  "directly to apply_boundary_condition function." << std::endl;
        exit(-1);
        break;
      }
    }
  }
}

}  // namespace kern_interp
