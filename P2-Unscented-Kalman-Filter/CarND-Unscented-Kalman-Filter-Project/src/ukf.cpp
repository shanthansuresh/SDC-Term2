#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = 3.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = 1;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
 
  // Dimension of the state vector
  n_x_ = 5;

  //Dimension of the augmentation vector
  n_aug_ = (n_x_ + 2);

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Predicted sigma point matrix
  Xsig_pred_ = MatrixXd(n_x_, (2*n_aug_) + 1);

  // Weights of sigma points
  weights_ = VectorXd((2*n_aug_)+1);

  // set weights
  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight_0 = 0.5/(n_aug_+lambda_);
    weights_(i) = weight_0;
  }

  // Initialize the state vector
  x_ = VectorXd(5, 1);
  x_ << 0, 0, 0, 0, 0;

  // Initialize the state covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
}

UKF::~UKF() {}

double normalize_angle(double phi) {
  double phi_norm;

  if (phi > M_PI) {
      phi_norm = (int(phi-M_PI)%int(2*M_PI)) - M_PI;
  } else if (phi < (-1*M_PI)) {
      phi_norm = (int(phi+M_PI)%int(2*M_PI)) + M_PI;
  }

  return phi_norm;
}

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

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    x_ = VectorXd(5);
    x_ << 1, 1, 1, 1, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float px_ = meas_package.raw_measurements_[0] *  cos(meas_package.raw_measurements_[1]);
      float py_ = meas_package.raw_measurements_[0] *  sin(meas_package.raw_measurements_[1]);
      x_ << px_, py_, 0, 0, 0;  
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;  
    }

    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */
  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  //Call the predict function.
  Prediction(dt); 

  /*
   * Update
   */
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package); 
  } else {
      UpdateLidar(meas_package); 
  }

  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  
  /*
   * (a) Initialize the matrices for sigma points calculation.
   */
  
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //square root matrix
  MatrixXd A = MatrixXd(n_aug_, n_aug_);

  //create row vector of ones
  MatrixXd ones_nAug = MatrixXd(1, n_aug_);
  ones_nAug.setOnes();

  /*
   * (b) Calculate the augmented points.
   */

  //create augmented mean state
  x_aug.setZero();
  x_aug.head(n_x_) << x_;

  //create augmented covariance matrix
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //calculate square root of P
  A = P_aug.llt().matrixL();

  //create augmented sigma points matrix
  //The matrix version is taken from https://github.com/vxy10
  Xsig_aug << x_aug,
              (x_aug * ones_nAug) + (sqrt(lambda_+n_aug_) * A),
              (x_aug * ones_nAug) - (sqrt(lambda_+n_aug_) * A);

  /*
   * (c) Predict sigma points
   */
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;
    //double yaw_p;
    double dt2 = delta_t * delta_t;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    } 

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*dt2 * cos(yaw);
    py_p = py_p + 0.5*nu_a*dt2 * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*dt2; //NEW
    
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  /*
   * (d) Predict state mean and covariance
   */

  //state mean
  x_ = Xsig_pred_ * weights_;

  //state covariance
  MatrixXd Wts_diag = MatrixXd(2*n_aug_+1,2*n_aug_+1);
  Wts_diag  = MatrixXd(weights_.asDiagonal());
  MatrixXd Ones_nA = MatrixXd(1, 2*n_aug_+1);
  Ones_nA.setOnes();
 
  P_ =(Xsig_pred_-x_*Ones_nA)*Wts_diag*(Xsig_pred_-x_*Ones_nA).transpose();
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

  //set measurement dimension, lidar can measure px, py
  int n_z = 2;

  //create matrix for sigma points in measurement space  
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  MatrixXd S = MatrixXd(n_z,n_z);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  MatrixXd R = MatrixXd(n_z,n_z);
  VectorXd z_pred = VectorXd(n_z);
  VectorXd z = VectorXd(n_z);
  VectorXd z_diff = VectorXd(n_z);
  MatrixXd Ones_A = MatrixXd(1,2*n_aug_+1);
  Ones_A.setOnes();
  MatrixXd Z_diff = MatrixXd(n_z, 2*n_aug_+1);
  MatrixXd K = MatrixXd(n_x_, n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < ((2*n_aug_) + 1); i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }
 
  //mean predicted measurement
  z_pred.fill(0.0);
  z_pred = Zsig * weights_;


  S =(Zsig-z_pred*Ones_A)*MatrixXd(weights_.asDiagonal())*(Zsig-z_pred*Ones_A).transpose();
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;

  Z_diff = (Zsig-z_pred*Ones_A);
  Tc = (Xsig_pred_-x_*Ones_A)*MatrixXd(weights_.asDiagonal())*Z_diff.transpose();
  K = Tc * S.inverse();

  z =  meas_package.raw_measurements_;
  z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  NIS_laser_ = (z - z_pred).transpose()*S*(z - z_pred);

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

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space  
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    double c1 = p_x*p_x+p_y*p_y;
    
    // measurement model
    Zsig(0,i) = sqrt(c1);

    //Check if c2 is really small.
    if (Zsig(0,i) < 0.00001) {
      Zsig(0,i) = 0.00001;
      Zsig(1,i) = atan2(0.001, 0.001);         //phi
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,1);   //r_dot
    } else {
      Zsig(1,i) = atan2(p_y,p_x);                   //phi
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);   //r_dot
    }
  }
 
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred = Zsig * weights_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  MatrixXd S = MatrixXd(n_z,n_z);
  VectorXd z_diff = VectorXd(n_z);
  VectorXd x_diff = VectorXd(n_x_);
  VectorXd z = VectorXd(n_z, 1);
  MatrixXd R = MatrixXd(n_z,n_z);
  MatrixXd Ones_A = MatrixXd(1,2*n_aug_+1);
  Ones_A.setOnes();
  MatrixXd Z_diff = MatrixXd(n_z, 2*n_aug_+1);

  /*
   * calculate measurement covariance matrix S.
   * calculate cross correlation matrix.
   */
  S =(Zsig-z_pred*Ones_A)*MatrixXd(weights_.asDiagonal())*(Zsig-z_pred*Ones_A).transpose();
  Z_diff = (Zsig-z_pred*Ones_A);

  for (int i=0;i<n_z;i++){
      Z_diff(i,1) = atan2(sin(Z_diff(i,2)),cos(Z_diff(i,2)));

  }
  
  Tc = (Xsig_pred_-x_*Ones_A)*MatrixXd(weights_.asDiagonal())*Z_diff.transpose();
  
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = normalize_angle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
}
