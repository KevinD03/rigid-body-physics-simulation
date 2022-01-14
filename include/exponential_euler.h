#include <Eigen/Dense>
#include <EigenTypes.h>
#include <rodrigues.h>
#include <iostream>

//Input:
//  q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles. 
//      The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//  qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and 
//         the second 3 are the world space linear velocity.
//  dt - the integration time step
//  masses - a vector to mass matrices for each rigid body
//  forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//           while the second 3 doubles are the linear forces.
//Output:
//  q - updated generalized coordinates 
//  qdot - updated generalized velocities 

inline void exponential_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces) {

    // formula on lecture 7 slide 64
    Eigen::Matrix3d R, RTR;
    Eigen::Vector3d pt, pdot, w, external_torque, external_force;

    for (int i = 0; i < q.rows() / 12; i++) {
        // use eigenn map as handout instructed.
        R = Eigen::Map<const Eigen::Matrix3d>(q.segment(12 * i, 9).data());
        double mass = masses[i].block(3, 3, 3, 3)(0, 0);
        RTR = R * masses[i].block(0 ,0, 3, 3) * R.transpose();

        external_torque = forces.segment(6 * i, 3);
        external_force = forces.segment(6 * i + 3, 3);
        pt = q.segment(12 * i + 9, 3);
        pdot = qdot.segment(6 * i + 3, 3);
        w = qdot.segment(6 * i, 3);

        // velocity update
        // see velocity update equation
        Eigen::Vector3d omega_next = w - (RTR.inverse() * (dt * w.cross(RTR * w)) + dt * external_torque);
        Eigen::Vector3d pdot_next = (mass * pdot + dt * external_force) / mass;

        qdot.segment(6 * i, 3) = omega_next;
        qdot.segment(6 * i + 3, 3) = pdot_next;
 
        // position update
        // see position update equation
        Eigen::Matrix3d rotation, Rt_next;
        rodrigues(rotation, w * dt);
        Rt_next = rotation * R;

        for (int j = 0; j < 3; j++) {
            q.segment(j * 3, 3) = Rt_next.block(0, j, 3, 1);
        }
        q.segment(9, 3) = pt + dt * pdot;

    }

}