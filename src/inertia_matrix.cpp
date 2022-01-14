#include <inertia_matrix.h>
#include <cassert>
#include <iostream>

//Input:
//  V - the nx3 matrix of vertices.
//  F - the mx3 matrix of triangle vertex indices.
//  density - the material density.
//Output:
//  I - the 3x3 angular inertia matrix
//  center - the center of mass of the object
//  mass - the total mass of the object

//compute inertia matrix and volume by integrating on surfaces

void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d & center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density) {
	
	// Lecture 7 slide 33
    I.setZero();

    // see piazza @191, @190, @131, @139

    // see A4 V_membrane_corotational for integration over triangle surface
	Eigen::Vector3d X0, X1, X2, N;

	for (int i = 0; i < F.rows(); i++) {

		X0 = V.row(F(i, 0));
		X1 = V.row(F(i, 1));
		X2 = V.row(F(i, 2));

		N = (X1 - X0).cross(X2 - X0).normalized();

        mass += 1.0 / 6.0 * N[0] * (X0[0] + X1[0] + X2[0]);
	}

    // lecture 7 slide 35 and 41
    for (int i = 0; i < F.rows(); i++) {
        X0 = V.row(F(i, 0));
        X1 = V.row(F(i, 1));
        X2 = V.row(F(i, 2));

        // norm of current triangle
        N = (X1 - X0).cross(X2 - X0).normalized();

        center[0] = 1.0 / 12.0 * N[0] * (X0[0] * X0[0] + X1[0] * X1[0] + X2[0] * X2[0] + X1[0] * X2[0] + X0[0] * X1[0] + X0[0] * X2[0]) / mass;
        center[1] = 1.0 / 12.0 * N[1] * (X0[1] * X0[1] + X1[1] * X1[1] + X2[1] * X2[1] + X1[1] * X2[1] + X0[1] * X1[1] + X0[1] * X2[1]) / mass;
        center[2] = 1.0 / 12.0 * N[2] * (X0[2] * X0[2] + X1[2] * X1[2] + X2[2] * X2[2] + X1[2] * X2[2] + X0[2] * X1[2] + X0[2] * X2[2]) / mass;

        I(0, 0) += 1.0 / 20.0 * N[1] * (X0[1] * X0[1] * X0[1] + X0[1] * X0[1] * (X1[1] + X2[1]) + (X1[1] + X2[1]) * (X1[1] * X1[1] + X2[1] * X2[1]) + X0[1] * (X1[1] * X1[1] + X1[1] * X2[1] + X2[1] * X2[1]))
                    + 1.0 / 20.0 * N[2] * (X0[2] * X0[2] * X0[2] + X0[2] * X0[2] * (X1[2] + X2[2]) + (X1[2] + X2[2]) * (X1[2] * X1[2] + X2[2] * X2[2]) + X0[2] * (X1[2] * X1[2] + X1[2] * X2[2] + X2[2] * X2[2])) - mass * (center[1] * center[1] + center[2] * center[2]);
        
        I(0, 1) += -1.0 / 60.0 * N[0] * (X0[1] * (X1[0] * X1[0] + X1[0] * X2[0] + X2[0] * X2[0]) + X1[1] * (3 * X1[0] * X1[0] + 2 * X1[0] * X2[0] + X2[0] * X2[0]) + (X1[0] * X1[0] + 2 * X1[0] * X2[0] + 3 * X2[0] * X2[0]) * X2[1] + X0[0] * X0[0] * (3 * X0[1] + X1[1] + X2[1]) + X0[0] * (2 * X0[1] * (X1[0] + X2[0]) + X1[0] * (2 * X1[1] + X2[1]) + X2[0] * (X1[1] + 2 * X2[1]))) + mass * center[0] * center[1];

        I(0, 2) += -1.0 / 60.0 * N[2] * (X0[2] * X0[2] * (X1[0] + X2[0]) + X1[2] * X1[2] * (3 * X1[0] + X2[0]) + 2 * X1[2] * (X1[0] + X2[0]) * X2[2] + (X1[0] + 3 * X2[0]) * X2[2] * X2[2] + X0[0] * (3 * X0[2] * X0[2] + X1[2] * X1[2] + X1[2] * X2[2] + X2[2] * X2[2] + 2 * X0[2] * (X1[2] + X2[2])) + X0[2] * (X1[0] * (2 * X1[2] + X2[2]) + X2[2] * (X1[2] + 2 * X2[2]))) + mass * center[2] * center[0];

        I(1, 0) += -1.0 / 60.0 * N[0] * (X0[1] * (X1[0] * X1[0] + X1[0] * X2[0] + X2[0] * X2[0]) + X1[1] * (3 * X1[0] * X1[0] + 2 * X1[0] * X2[0] + X2[0] * X2[0]) + (X1[0] * X1[0] + 2 * X1[0] * X2[0] + 3 * X2[0] * X2[0]) * X2[1] + X0[0] * X0[0] * (3 * X0[1] + X1[1] + X2[1]) + X0[0] * (2 * X0[1] * (X1[0] + X2[0]) + X1[0] * (2 * X1[1] + X2[1]) + X2[0] * (X1[1] + 2 * X2[1]))) + mass * center[0] * center[1];
        
        I(1, 1) += 1.0 / 20.0 * N[0] * (X0[0] * X0[0] * X0[0] + X0[0] * X0[0] * (X1[0] + X2[0]) + (X1[0] + X2[0]) * (X1[0] * X1[0] + X2[0] * X2[0]) + X0[0] * (X1[0] * X1[0] + X1[0] * X2[0] + X2[0] * X2[0]))
                    + 1.0 / 20.0 * N[2] * (X0[2] * X0[2] * X0[2] + X0[2] * X0[2] * (X1[2] + X2[2]) + (X1[2] + X2[2]) * (X1[2] * X1[2] + X2[2] * X2[2]) + X0[2] * (X1[2] * X1[2] + X1[2] * X2[2] + X2[2] * X2[2])) - mass * (center[2] * center[2] + center[0] * center[0]);

        I(1, 2) += -1.0 / 60.0 * N[1] * (X0[2] * (X1[1] * X1[1] + X1[1] * X2[1] + X2[1] * X2[1]) + X1[2] * (3 * X1[1] * X1[1] + 2 * X1[1] * X2[1] + X2[1] * X2[1]) + (X1[1] * X1[1] + 2 * X1[1] * X2[1] + 3 * X2[1] * X2[1]) * X2[2] + X0[1] * X0[1] * (3 * X0[2] + X1[2] + X2[2]) + X0[1] * (2 * X0[2] * (X1[1] + X2[1]) + X1[1] * (2 * X1[2] + X2[2]) + X2[1] * (X1[2] + 2 * X2[2]))) + mass * center[1] * center[2];
        
        I(2, 0) += -1.0 / 60.0 * N[2] * (X0[2] * X0[2] * (X1[0] + X2[0]) + X1[2] * X1[2] * (3 * X1[0] + X2[0]) + 2 * X1[2] * (X1[0] + X2[0]) * X2[2] + (X1[0] + 3 * X2[0]) * X2[2] * X2[2] + X0[0] * (3 * X0[2] * X0[2] + X1[2] * X1[2] + X1[2] * X2[2] + X2[2] * X2[2] + 2 * X0[2] * (X1[2] + X2[2])) + X0[2] * (X1[0] * (2 * X1[2] + X2[2]) + X2[2] * (X1[2] + 2 * X2[2]))) + mass * center[2] * center[0];

        I(2, 1) += -1.0 / 60.0 * N[1] * (X0[2] * (X1[1] * X1[1] + X1[1] * X2[1] + X2[1] * X2[1]) + X1[2] * (3 * X1[1] * X1[1] + 2 * X1[1] * X2[1] + X2[1] * X2[1]) + (X1[1] * X1[1] + 2 * X1[1] * X2[1] + 3 * X2[1] * X2[1]) * X2[2] + X0[1] * X0[1] * (3 * X0[2] + X1[2] + X2[2]) + X0[1] * (2 * X0[2] * (X1[1] + X2[1]) + X1[1] * (2 * X1[2] + X2[2]) + X2[1] * (X1[2] + 2 * X2[2]))) + mass * center[1] * center[2];

        I(2, 2) += 1.0 / 20.0 * N[0] * (X0[0] * X0[0] * X0[0] + X0[0] * X0[0] * (X1[0] + X2[0]) + (X1[0] + X2[0]) * (X1[0] * X1[0] + X2[0] * X2[0]) + X0[0] * (X1[0] * X1[0] + X1[0] * X2[0] + X2[0] * X2[0]))
                    + 1.0 / 20.0 * N[1] * (X0[1] * X0[1] * X0[1] + X0[1] * X0[1] * (X1[1] + X2[1]) + (X1[1] + X2[1]) * (X1[1] * X1[1] + X2[1] * X2[1]) + X0[1] * (X1[1] * X1[1] + X1[1] * X2[1] + X2[1] * X2[1])) - mass * (center[0] * center[0] + center[1] * center[1]);
                
    }

}