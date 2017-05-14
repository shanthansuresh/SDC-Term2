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
