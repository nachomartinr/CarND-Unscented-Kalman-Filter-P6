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
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.7; /* TODO: Tune */

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6; /* TODO: Tune */

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

  // initially set to false
  is_initialized_ = false;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2*n_aug_+1);

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Sigma points weights 
  weights_ = VectorXd(2*n_aug_+1);  
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  weights_.tail(2*n_aug_).fill(0.5 / (lambda_ + n_aug_));


  // the current NIS for rada
  NIS_radar_ = 0;

  // the current NIS for laser
  NIS_laser_ = 0;

  // time when the state is true, in us
  time_us_ = 0;

  //laser measurement dimensions
  n_z_las_ = 2;

  //radar measurement dimensions
  n_z_rad_ = 3;

  // lidar measurement noise covariance matrix
  R_laser_ = MatrixXd(n_z_las_, n_z_las_);
  R_laser_ << std_laspx_*std_laspx_,                    0,
                                 0, std_laspy_*std_laspy_;

  R_radar_ = MatrixXd(n_z_rad_, n_z_rad_);
  R_radar_ << std_radr_*std_radr_,                       0,                       0,
                                0, std_radphi_*std_radphi_,                       0,
                                0,                       0, std_radrd_ * std_radrd_;
}

UKF::~UKF() {}


void UKF::InitializeState(const MeasurementPackage meas_package) {

    P_ <<   1,  0,  0,   0,  0,
            0,  1,  0,   0,  0,
            0,  0,  15,  0,  0,
            0,  0,  0,  15,  0,
            0,  0,  0,   0,  0.5;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);

      double px = rho * cos(phi);
      double py = rho * sin(phi);

      // we don't have enough information to initialize v (tangential velocity) from rho_dot.
      x_ << px, py, 0, 0, 0;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);

      x_ << px, py, 0, 0, 0;

      /* With a lidar measurement, the initial variance or uncertainty in p​x and py
​​         would be equal to the variance of lidar sensor*/        
      P_(0,0) = std_laspx_*std_laspx_;
      P_(1,1) = std_laspy_*std_laspy_;
    }

    time_us_ = meas_package.timestamp_;

    is_initialized_ = true; 
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage meas_package) {

  if (!is_initialized_) {		
      InitializeState(meas_package);
      return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;

  Prediction(delta_t);  

  time_us_ = meas_package.timestamp_;

  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_) {

    UpdateRadar(meas_package);


  }
  else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_) {

    UpdateLidar(meas_package);

  }

}

/**
 * Generate augmented sigma points from current, a posteriori state
 * @param {MatrixXd} Xsig_aug Generated augmented sigma points matrix,
 * with size (n_aug_, 2*n_aug_+1).
 */
void UKF::GenerateSigmaPoints(MatrixXd& Xsig_aug) const {  

  //create augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  //calculate the scaled (by the spread factor) sigma points vectors
  MatrixXd P_sigmas =  (sqrt(lambda_+n_aug_) * L);

  //create augmented sigma points (vectorized)
  Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.middleCols(1, n_aug_) = P_sigmas.colwise() + x_aug;
  Xsig_aug.middleCols(n_aug_+1, n_aug_) = (-P_sigmas).colwise() + x_aug;

}

/**
 * Apply the process transformation function to the sigma points, including
 * the process noise.
 * @param {MatrixXd} Xsig_aug Augmented sigma points matrix.
 * @param {double} delta_t The change in time (in seconds) since the last
 * a posteriori state.
 */
void UKF::PredictSigmaPoints(const MatrixXd& Xsig_aug, const double delta_t) {

  //predict sigma points
  for (int i = 0; i < (2 * n_aug_ + 1); ++i) {
      //extract values for better readability
      double px       = Xsig_aug(0,i);
      double py       = Xsig_aug(1,i);
      double v        = Xsig_aug(2,i);
      double yaw      = Xsig_aug(3,i);
      double yawd     = Xsig_aug(4,i);
      double nu_a     = Xsig_aug(5,i);
      double nu_yawdd =  Xsig_aug(6,i);      
      
      double px_p;      
      double py_p;
      
      //avoid division by zero
      if (fabs(yawd) > 0.001) {
        px_p = px + (v/yawd) * (sin(yaw + yawd*delta_t) - sin(yaw));
        py_p = py + (v/yawd) * (-cos(yaw + yawd*delta_t) + cos(yaw));
      }
      else {
        px_p = px + v*cos(yaw)*delta_t;
        py_p = py + v*sin(yaw)*delta_t;
      }
      
      double v_p = v;
      double yaw_p = yaw + yawd*delta_t;
      double yawd_p = yawd;
      
      //add process noise 
      px_p += 0.5*delta_t*delta_t*cos(yaw)*nu_a;
      py_p += 0.5*delta_t*delta_t*sin(yaw)*nu_a;
      v_p += delta_t*nu_a;
      yaw_p += 0.5*delta_t*delta_t*nu_yawdd;
      yawd_p += delta_t*nu_yawdd;
      
      Xsig_pred_(0,i) = px_p;
      Xsig_pred_(1,i) = py_p;
      Xsig_pred_(2,i) = v_p;
      Xsig_pred_(3,i) = yaw_p;
      Xsig_pred_(4,i) = yawd_p;      
  }
}

/* Constraint angle to [-PI, PI]
 * @param {double} a Angle to wrap */
static inline double WrapAngle(const double & a) {

  double wrapped;
  if (( a >= -M_PI ) && (a <= M_PI )) {
    wrapped = a;
  } 
  else if ( a < 0) {
    wrapped = fmod(a - M_PI, 2*M_PI) + M_PI;
  }
  else {
    wrapped = fmod(a + M_PI, 2*M_PI) - M_PI;
  }

  return wrapped;
}

/**
 * Calculate the predicted mean and covariance of the a priori state from the predicted
 * sigma points.
 */
void UKF::PredictedMeanAndCovariance() {

  //predict state mean (vectorized)
  x_ = Xsig_pred_*weights_; // result is a vector of means
  
  //predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < (2*n_aug_ + 1); ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    /* normalize yaw angle to [-pi, pi] */
    x_diff(3) = WrapAngle(x_diff(3));
    
    P_ += weights_(i) * x_diff * x_diff.transpose();
    
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(const double delta_t) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // generate augmented sigma points from the last 'a posteriori' state
  MatrixXd Xsig_aug;
  GenerateSigmaPoints(Xsig_aug);

  /* When dt is small, skip the prediction step.*/
  if (delta_t >= 0.001) {
    // apply f(x,nu) to augmented sigma points to get the predicted sigma points.  
    PredictSigmaPoints(Xsig_aug, delta_t); 

    // Get mean and covariance from predicted sigma points
    PredictedMeanAndCovariance();
  }
  else {

    /* Sigma points must be updated
    after each new 'a posteriori' state, even if
    no  since they are used in the measurement step, */
    Xsig_pred_ = Xsig_aug.topRows(n_x_);
  }
}



void UKF::PredictLidarMeasurement(VectorXd& z_pred,
                                  MatrixXd& S,
                                  MatrixXd& S_inv,
                                  MatrixXd& K) {

  //transform sigma points into measurement space  
  MatrixXd Zsig = MatrixXd(n_z_las_, 2*n_aug_+1);
  for (int i=0; i<(2*n_aug_+1); ++i) {
    // sigma points in measurement space 
    Zsig(0,i)  = Xsig_pred_(0,i); /* px */
    Zsig(1,i)  = Xsig_pred_(1,i); /* py */
  }
  
  // calculate mean predicted measurement (vectorized)
  z_pred = Zsig*weights_;
  
  // measurement covariance matrix S
  S = MatrixXd::Zero(n_z_las_, n_z_las_);

  // Cross-correlation matrix of sigma points in
  // state space and measurmement space, Tc.
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_las_);  

  for (int i = 0; i<(2*n_aug_+1); ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;      
            
    S += weights_(i)*z_diff*z_diff.transpose();

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc += weights_(i)*x_diff*z_diff.transpose();
  }  

  S += R_laser_;

  S_inv = S.inverse();

  //Kalman gain K
  K = Tc * S_inv;
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  Calculate the lidar NIS.
  */

  VectorXd z_pred;
  MatrixXd S, S_inv;
  MatrixXd K;
  PredictLidarMeasurement(z_pred, S, S_inv, K);

  VectorXd z = VectorXd(n_z_las_);
  z << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1);

  //residual (innovation)
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // calculate the lidar Normalized Innovation Squared (NIS)
  NIS_laser_ = z_diff.transpose() * S_inv * z_diff;
}


void UKF::PredictRadarMeasurement(VectorXd& z_pred,
                                  MatrixXd& S,
                                  MatrixXd& S_inv,
                                  MatrixXd& K) {

  //transform sigma points into measurement space  
  MatrixXd Zsig = MatrixXd(n_z_rad_, 2 * n_aug_ + 1);
  for (int i=0; i<(2*n_aug_+1); ++i) {
    // sigma points in measurement space 
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;    

    /* avoid division by 0 */
    if (fabs(px) < 0.0001)
          px = (px < 0)? -0.0001 : 0.0001;
    if (fabs(py) < 0.0001)
          py = (py < 0)? -0.0001 : 0.0001;

    // measurement model
    Zsig(0,i) = sqrt(px*px + py*py);                     //rho
    Zsig(1,i) = atan2(py, px);                           //phi
    Zsig(2,i) = (px*v1 + py*v2) / sqrt(px*px + py*py);   //rho_dot
  }
  
  // calculate mean predicted measurement (vectorized)
  z_pred = Zsig*weights_;
  
  // measurement covariance matrix S
  S = MatrixXd::Zero(n_z_rad_, n_z_rad_);

  // Cross-correlation matrix of sigma points in
  // state space and measurmement space, Tc.
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_rad_);  

  for (int i = 0; i<(2*n_aug_+1); ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;      
            
    S += weights_(i)*z_diff*z_diff.transpose();

    //angle normalization
    z_diff(1) = WrapAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    x_diff(3) = WrapAngle(x_diff(3));

    Tc += weights_(i)*x_diff*z_diff.transpose();
  }  

  S += R_radar_;

  S_inv = S.inverse();

  //Kalman gain K
  K = Tc * S_inv;
}



/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  Calculate the radar NIS.
  */
  VectorXd z_pred;
  MatrixXd S, S_inv;
  MatrixXd K;
  PredictRadarMeasurement(z_pred, S, S_inv, K);

  VectorXd z = VectorXd(n_z_rad_);
  z << meas_package.raw_measurements_(0),
       meas_package.raw_measurements_(1),
       meas_package.raw_measurements_(2);

  //residual (innovation)
  VectorXd z_diff = z - z_pred;

  /* normalize phi angle to [-pi, pi] */
  z_diff(1) = WrapAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // calculate the radar Normalized Innovation Squared (NIS)
  NIS_radar_ = z_diff.transpose() * S_inv * z_diff;
}
