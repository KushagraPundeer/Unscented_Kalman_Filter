#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define EPS 0.001

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Set is_initialization as False
  is_initialized_ = false; 

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.75;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // Set state dim.
  n_x_ = 5;     
  // Set augmented dim.        
  n_aug_ = n_x_ + 2;    

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, 5);

  // Start time
  time_us_ = 0;          

  // Design Parameter Lambda
  lambda_ = 3 - n_aug_;


  // Predicted sigma Points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1); 

  // Weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < (2 * n_aug_ + 1); i++) {
      double weight = 0.5/ (n_aug_ + lambda_);
      weights_(i) = weight;
  }

  // Intialize Noise Matrix
  R_radar_ = MatrixXd(3,3); 
  R_lidar_ = MatrixXd(2,2); 

  R_radar_ << std_radr_*std_radr_, 0, 0,
             0, std_radphi_*std_radphi_, 0,
             0, 0,std_radrd_*std_radrd_;

  R_lidar_ << std_laspx_*std_laspx_,0,
             0,std_laspy_*std_laspy_;


  NIS_radar_ = 0;        // NIS for Radar
  NIS_lidar_ = 0;        // NIS for Laser

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_){

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){

      double rho = meas_package.raw_measurements_[0];      // Range rho
      double phi = meas_package.raw_measurements_[1];      // Bearing phi
      double rho_dot = meas_package.raw_measurements_[2];  // Rho dot (velocity)

      // Polar to cartesian
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);

      double v = sqrt(vx*vx+vy*vy);
      x_ << px, py, v, 0, 0;

    }

    else if (meas_package.sensor_type_ == MeasurementPackage::LASER){

      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

      if (fabs(x_(0)) < EPS and fabs(x_(1)) < EPS){
        x_(0) = EPS;
        x_(1) = EPS;
      }

    }

    // Set P_
    P_ << 0.1, 0, 0, 0, 0,
          0, 0.1, 0, 0, 0,
          0, 0, 20, 0, 0,
          0, 0, 0, 20, 0,
          0, 0, 0, 0, 20;
    
    time_us_ = meas_package.timestamp_;
    is_initialized_= true;
    
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; // delta_t in seconds
  time_us_ = meas_package.timestamp_;

  if (delta_t > 0.001){
    Prediction(delta_t);
  }
  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
    UpdateRadar(meas_package);
  } 
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Lesson 7 Section 18: Assign. - 2
  // Create Augmented Sigma Points

  // Augmented mean vector
  VectorXd x_aug_ = VectorXd(7);
  //Augmented State covariance matrix
  MatrixXd P_aug_ = MatrixXd(7,7);

  // Sigma Point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  double delta_t_2 = delta_t * delta_t;

  // Create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  // Create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;

  // Create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  // Create augmented sigma points
  Xsig_aug_.col(0)= x_aug_;
  for (int i = 0; i < n_aug_; i++){
    Xsig_aug_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  // Predict Sigma Points
  for (int i = 0; i < (2 * n_aug_+ 1); i++)
  {
    // Extract values for better readablity
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    // Predicted State values
    double px_p, py_p;

    // Avoid division by zero 
    if (fabs(yawd) > EPS)
    {
      px_p = p_x + (v/yawd) * (sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + (v/yawd) * (cos(yaw) - cos(yaw+yawd*delta_t));

    }
    else
    {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);

    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // Add noise
    px_p = px_p + (0.5 * nu_a* delta_t_2 * cos(yaw));
    py_p = py_p + (0.5 * nu_a* delta_t_2* sin(yaw));
    v_p = v_p + nu_a * delta_t;
    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t_2;
    yawd_p = yawd_p + nu_yawdd * delta_t;


    // Write predicted sigma points into right col.
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;

  }
  // Lesson 7, section 24: Prediction Mean and Covariance Assign - 2
  // Set weights 
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i < (2 * n_aug_ + 1); i++) {
      double weight = 0.5 / (n_aug_ + lambda_);
      weights_(i) = weight;
  }

  // Predicted State Mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    // Iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // Predicated state covariance matrix
  P_.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++){
    // Iterate over sigma points
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle Noramalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();

  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // Initialize measurement dimension
  int n_z = 2;

  // Matrix for sigma points
  MatrixXd z_sig = MatrixXd(n_z, 2 * n_aug_ + 1);
  z_sig = Xsig_pred_.topRows(2);

  // Covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // Predicted mean measurement
  VectorXd z_pred = VectorXd(n_z);
  VectorXd z = meas_package.raw_measurements_;
  z_sig.fill(0.0);
  z_pred.fill(0.0);

  // Set weights 
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i< (2 * n_aug_ + 1); i++) {
      double weight = 0.5/ (n_aug_ + lambda_);
      weights_(i) = weight;
  }

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    z_sig(0, i) = px;
    z_sig(1, i) = py;

    // Predicted mean measurement
    z_pred = z_pred + weights_(i) * z_sig.col(i);
  } 

  // Cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  S.fill(0.0);
  Tc.fill(0.0);

  //Covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = z_sig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  S = S + R_lidar_;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - z_pred;

  // Update state mean and Covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  //NIS Lidar Update
  NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // Lesson 7, section 27 - Predict Radar Measurement Assign - 2

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd z_sig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  { //2n+1 sigma points

    // extract values for better readability
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    float rho = sqrt(p_x * p_x + p_y * p_y);                             //rho
    float theta = atan2(p_y, p_x);                                       //theta
    float rho_dot = (p_x * v1 + p_y * v2) / rho;                         //rho_dot

    //handling if px and py are near zero
    if (p_x < 0.001 && p_y < 0.001)
    {
      theta = 0;
      rho_dot = 0;
    }

    z_sig(0, i) = rho;
    z_sig(1, i) = theta;
    z_sig(2, i) = rho_dot;
  }

  
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    z_pred = z_pred + weights_(i) * z_sig.col(i);
  }

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  { //2n+1 sigma points
    //residual
    VectorXd z_diff = z_sig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1) >  M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_radar_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  { //2n+1 sigma points

    //residual
    VectorXd z_diff = z_sig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1) >  M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3) >  M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;


  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1) >  M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

}
