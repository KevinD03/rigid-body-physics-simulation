#include <rigid_to_world.h>

//Input:
//  R - rotation matrix for rigid body
//  p - world space position of center-of-mass
//  X -  undeformed position at which to compute the Jacobian. 
//Output:
//  x - the world space position of the undeformed point X.

void rigid_to_world(Eigen::Vector3d &x, 
               Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
               Eigen::Ref<const Eigen::Vector3d> X) {

    // R 3*3
    // p 3*1
    // see lecture 7 slide 14

    x = R * X + p;
                       
}