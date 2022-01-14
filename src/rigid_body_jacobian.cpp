#include <rigid_body_jacobian.h>

//Input:
//  R - rotation matrix for rigid body
//  p - world space position of center-of-mass
//  X -  undeformed position at which to compute the Jacobian. 
//Output:
//  J - the rigid body jacobian acting on the undeformed space point X.

void rigid_body_jacobian(Eigen::Matrix36d &J, 
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
                         Eigen::Ref<const Eigen::Vector3d> x) {
    
    // lecture 7 slide 51

    Eigen::Matrix3d cross_product_M;
    cross_product_M << 0.0, -x(2), x(1),
                    x(2), 0.0, -x(0),
                    -x(1), x(0), 0.0;

    // X_hat_t = {[X]^T, I}
    Eigen::Matrix36d X_hat_t;

    X_hat_t << cross_product_M.transpose(), Eigen::Matrix3d::Identity();
    
    // Rt = { R^T, 0    
    //        0, I }
    Eigen::Matrix66d Rt;
    Rt << R.transpose(), Eigen::Matrix3d::Zero(),
        Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Identity();

    J = R * X_hat_t * Rt;
}

