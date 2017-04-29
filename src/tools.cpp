#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RMSE here.
   */

  // The RMSE is a 4 elements vector which
  // contains the RMSE of: x, y, vx and vy.
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // Check the validity of the inputs:
  // 1. The estimation vector size should not be zero.
  if (estimations.size() == 0) {
    std::cout << "The estimation vector size should not be zero.\n";
    return rmse;
  }

  // 2. The estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()) {
    std::cout << "The estimation vector size should equal ground truth vector size.\n";
    return rmse;
  }

  // Accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];

    // Coefficient-wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculate a Jacobian here.
   */

  MatrixXd Hj(3,4);

  // Recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // Check division by zero
  float px2_py2 = px*px + py*py;
  float sqrt_px2_py2 = sqrt(px2_py2);
  float px2_py2_3_2 = pow(px2_py2, 1.5);

  if (px2_py2 == 0) {
    std::cout << __func__ << ": Cannot divide by zero.\n";
    return Hj;
  }

  float h_00 = px / sqrt_px2_py2;
  float h_01 = py / sqrt_px2_py2;
  float h_10 = - py / px2_py2;
  float h_11 = px / px2_py2;
  float h_20 = py * (vx*py - vy*px) / px2_py2_3_2;
  float h_21 = px * (vy*px - vx*py) / px2_py2_3_2;
  float h_22 = h_00;
  float h_23 = h_01;

  // Compute the Jacobian matrix
  Hj << h_00, h_01, 0.0, 0.0,
    h_10, h_11, 0.0, 0.0,
    h_20, h_21, h_22, h_23;

  return Hj;
}
