#include <rodrigues.h>
#include <cmath>

//  R - rotation matrix 
//  omega - angular velocity vector

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega) {
	// see details in 
	// https://www.youtube.com/watch?v=1RF7j-Yc21c&list=PLTkE7n2CwG_NOFhh1We-ZKoKPoQVNf7mi&index=2&t=895s
	// 15:15

	R.setZero();

	Eigen::Matrix3d I;
	I.setIdentity();

	// theta: angle of rotation
	double theta = omega.norm();

	// a is the axis of rotation
	Eigen::Vector3d a;
	a = omega / theta;

	// lecture 7 slide 18
	Eigen::Matrix3d cross_product_M;
	cross_product_M << 0.0, -a(2), a(1),
					a(2), 0.0, -a(0),
					-a(1), a(0), 0.0;

	R = I + sin(theta) * cross_product_M + (1 - cos(theta)) * cross_product_M * cross_product_M;

}