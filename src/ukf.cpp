#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

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

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;
  n_sig_ = 2 * n_aug_ + 1;

  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  weights_ = VectorXd(n_sig_);
  weights_(0) = lambda_ / (lambda_+n_aug_);
  for (int i = 1; i < n_sig_; i++) {
    weights_(i) = 0.5 / (lambda_+n_aug_);
  }

  R_laser_ = MatrixXd(2,2);
  R_laser_ << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    x_ = VectorXd(5);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      float px = rho * cos(theta);
      float py = rho * sin(theta);
      float vx = rho_dot * cos(theta);
      float vy = rho_dot * sin(theta);
      float v  = sqrt(vx * vx + vy * vy);
      x_ << px, py, v, 0, 0; //vx, vy;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ <<  meas_package.raw_measurements_[0],
                  meas_package.raw_measurements_[1],
                  0,0,0;
      if (x_(0) < 0.001 && x_(1) < 0.001) {
        x_(0) = x_(1) = 0.001;
      }
    }
    P_ = MatrixXd::Identity(5,5);
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);
  if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER && use_laser_) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR && use_radar_) {
    UpdateRadar(meas_package);
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
  //cout << "Prediction 1" << endl;
  // sigma points generation
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  MatrixXd A = P_aug.llt().matrixL();          // sqrt
  //cout << "Prediction 1b" << endl;

  Xsig_aug.fill(0.0);  
  Xsig_aug.col(0) = x_aug;
  //cout << "Prediction 1c " << endl;

  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1)         = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i+1+n_aug_)  = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  double delta_t2 = delta_t * delta_t;
  //cout << "Prediction 2" << endl;

  for (int i = 0; i < n_sig_; i++) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double cy = cos(yaw);
    double sy = sin(yaw);
    double arg = yaw+yawd*delta_t;

    double px2, py2;
    if (fabs(yawd) > 0.001) {
      px2 = px + v/yawd*(sin(arg)-sy) + 0.5*delta_t2*cy*nu_a;
      py2 = py + v/yawd*(cy-cos(arg)) + 0.5*delta_t2*sy*nu_a;
    } else {
      px2 = px + v*cy*delta_t + 0.5*delta_t2*cy*nu_a;
      py2 = py + v*sy*delta_t + 0.5*delta_t2*sy*nu_a;
    }
    double v2 = v + delta_t*nu_a;
    double yaw2 = arg + 0.5*delta_t2*nu_yawdd;
    double yawd2 = yawd + delta_t*nu_yawdd;

    Xsig_pred_(0,i) = px2;
    Xsig_pred_(1,i) = py2;
    Xsig_pred_(2,i) = v2;
    Xsig_pred_(3,i) = yaw2;
    Xsig_pred_(4,i) = yawd2;
  }
  //cout << "Prediction 3" << endl;

  // state and state covariance
  x_.fill(0.0);
  //for (int i = 0; i < n_sig_; i++) {
  //  x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  //}
  x_ = Xsig_pred_ * weights_;

  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    VectorXd xDiff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(xDiff(3));
    P_ = P_ + weights_(i) * xDiff * xDiff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //cout << "Update Lidar" << endl;
  int n_z = 2;
  // create matrix for sigma ponts in measurement space
  MatrixXd Zsig  = Xsig_pred_.block(0,0,n_z,n_sig_);
  UpdateUKF(meas_package, Zsig, n_z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // transform sigma points into measurement space
  //cout << "Update Radar" << endl;

  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);
  for (int i = 0; i < n_sig_; i++) {
    float px = Xsig_pred_(0, i);
    float py = Xsig_pred_(1, i);
    float v = Xsig_pred_(2, i);
    float yaw = Xsig_pred_(3, i);
    //float yawd = Xsig_pred_(4,i);

    float phi = sqrt(px*px+py*py);
    float rho = atan2(py,px);
    float phid = (px * cos(yaw) * v + py * sin(yaw) * v) / phi;
    Zsig(0,i) = phi;
    Zsig(1,i) = rho;
    Zsig(2,i) = phid;
  }
  UpdateUKF(meas_package, Zsig, n_z);
}

void UKF::UpdateUKF(MeasurementPackage& meas_package, MatrixXd& Zsig, int n_z) {
  // z_pred, S, R
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  MatrixXd *R;

  //z_pred.fill(0.0);
  //for (int i = 0; i < n_sig_; i++) {
  //  z_pred = z_pred + weights_(i) * Zsig.col(i);
  //}
  z_pred = Zsig * weights_;

  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    VectorXd zd = Zsig.col(i) - z_pred;
      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        NormalizeAngle(zd(1));
      }
    S = S + weights_(i) * zd * zd.transpose();
  }

  R = (meas_package.sensor_type_ == MeasurementPackage::LASER) ? &R_laser_ : &R_radar_;

  S = S + *R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    
    VectorXd vX = Xsig_pred_.col(i) - x_;
    NormalizeAngle(vX(3));

    VectorXd vZ = Zsig.col(i) - z_pred;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      NormalizeAngle(vZ(1));
    }
    Tc = Tc + weights_(i) * vX * vZ.transpose();
  }

  VectorXd z = meas_package.raw_measurements_;
  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - z_pred;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    NormalizeAngle(z_diff(1));
  }
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  //cout << "Update UKF " << x_ << endl << P_ <<endl;
}

void UKF::NormalizeAngle(double& theta) {
  if (theta < -M_PI || theta > M_PI) {
    theta = atan2(sin(theta), cos(theta));
  }
}