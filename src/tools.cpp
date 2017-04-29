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
  TODO:
    * Calculate a Jacobian here.
  */
}
