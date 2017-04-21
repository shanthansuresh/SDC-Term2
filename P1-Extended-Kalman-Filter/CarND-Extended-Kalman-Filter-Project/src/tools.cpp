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
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  int N;

  rmse << 0,0,0,0;
  N = estimations.size();
  if (N == 0) {
    std::cout << "Estimation vector size is 0" << std::endl;
    return rmse;
  }
    
  if (estimations.size() != ground_truth.size()) {
    std::cout << "Estimation and ground truth size should be same" << std::endl;
    return rmse;
  }

  //accumulate squared residuals
  VectorXd a(4);
  VectorXd b(4);
  for(int i=0; i < estimations.size(); ++i){
    a = estimations[i] - ground_truth[i];
    b = a.array()*a.array();
    rmse = rmse + b;
  }

  //calculate the mean
  rmse = (1.0/N)*rmse;
    
  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  //check division by zero
  if(fabs(c1) < 0.0001){
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }

  //compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}
