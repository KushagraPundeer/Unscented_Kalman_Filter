#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size() == 0){
  	cout << "Input Empty" << endl;
  	return rmse;
  }

  // check the validity of the following inputs:
  if (estimations.size() != ground_truth.size() ){
    cout << "Invalid estimation/ground_truth data" << endl;
    return rmse;
  }

  // Accumulate squared residuals
  for(unsigned int i = 0; i < estimations.size(); ++i){
  	VectorXd residual = estimations[i] - ground_truth[i];

  	// coefficient-wise multiplication
  	residual = residual.array() * residual.array();
  	rmse = rmse + residual;
  }

  // Mean Calculation
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt(); // Squared Root
  
  return rmse;


}