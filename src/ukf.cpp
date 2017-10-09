#include "ukf.h"
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

  // set intialized status to 0;
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_ << 0, 0, 0, 0, 0;

  // initial covariance matrix
  P_ = MatrixXd(5, 5);


  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, .1, 0,
        0, 0, 0, 0, .1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5; //1.55; // 5; //1.8;   // 5 // .5 // .2 // 1.8 // .55

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .5; //0.7;  // .2 //.05 // .2 // .7 // .55

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

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_(0) = (lambda_/(n_aug_ + lambda_));
  for (int i = 1; i < n_aug_*2 + 1; i++){
    weights_(i) = (.5/(n_aug_ + lambda_));
  }

  // initialize Xsig_pred to right size
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
  Xsig_pred_.fill(0.0);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (is_initialized_ == false){
    // If Radar, initialize on radar data
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      double ro = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double ro_d = meas_package.raw_measurements_[2];
      x_ << ro*cos(phi), ro*sin(phi), 0, 0, 0;
    }

    // If Lidar, initialize on lidar data
    if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      double x = meas_package.raw_measurements_[0];
      double y = meas_package.raw_measurements_[1];
      x_ << x, y, 0, 0, 0;
    }

    // initialize timer
    time_us_ = meas_package.timestamp_;

    // Set Initialized to true
    is_initialized_ = true;

  } else {
		// Calculate dt
    double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
    time_us_ = meas_package.timestamp_;
    // Call prediction
    Prediction(delta_t);

    // Call Lidar or Radar Update
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      UpdateRadar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      UpdateLidar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // --------- Generate Sigma Points --------------
  // Set up our augmented matrix/vector -for noise-
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);
  Xsig_aug.fill(0.0);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) =  0;
  x_aug(6) =  0;
  // MatrixXd Q = MatrixXd(2,2);
  // Q <<    pow(std_a_, 2), 0,
  //        0 , pow(std_yawdd_, 2);

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  //P_aug.bottomRightCorner(2,2) = Q;		// Code seems to hang for a moment using this
  P_aug(5,5) = std_a_*std_a_;					// These might be faster
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  MatrixXd A = P_aug.llt().matrixL();

  // Calculate Sigma Points
  Xsig_aug.col(0) = (x_aug);
  for (int i = 0; i < n_aug_ ; i++){
    Xsig_aug.col(i+1)           = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i+1 + n_aug_)  = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  // -------- Predict Sigma Points ---------------
  for (int i = 0 ; i < 2*n_aug_ + 1; i ++){
    // predicted mean
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double p_v = Xsig_aug(2, i);
    double p_yaw = Xsig_aug(3, i);
    double p_yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double px_p, py_p;

    // calculate x integral
    if (fabs(p_yawd) > 0.00001){
      px_p = p_x + p_v/p_yawd*( sin(p_yaw + p_yawd*delta_t) - sin(p_yaw));
      py_p = p_y + p_v/p_yawd*( - cos(p_yaw + p_yawd*delta_t) + cos(p_yaw));
    } else {
      px_p = p_x + p_v*cos(p_yaw)*delta_t;
      py_p = p_y + p_v*sin(p_yaw)*delta_t;
    }

    double v_p = p_v;
    double yaw_p = p_yaw + p_yawd*delta_t;
    double yawd_p = p_yawd;

    // noise
    px_p += 0.5*nu_a*cos(p_yaw)*pow(delta_t,2);
    py_p += 0.5*nu_a*sin(p_yaw)*pow(delta_t,2);
    v_p += nu_a*delta_t;

		yaw_p += 0.5*nu_yawdd*pow(delta_t,2);
		yawd_p += nu_yawdd*delta_t;

    // Set our predicted Sigma Points
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // -------------- Predict Mean and Covariance Matrix -------------------
  x_.fill(0.0);
  // Calculate Predicted Mean & store it in x_
  for (int i = 0; i < 2 * n_aug_ + 1; i ++){
    x_ +=  weights_(i)* Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  // Calculate Predicted Covariance
  for (int i = 0; i < 2 * n_aug_ + 1; i ++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
		x_diff(3) = normalizeAngle(x_diff(3));

		// Store our values in P_
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  MatrixXd Zsig = MatrixXd(2, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  // Create Lidar Measurement Model Sigma Points
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }

	VectorXd z_pred = VectorXd(2);
	// Calculate Mean for Measurement Model
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    z_pred += weights_(i) * Zsig.col(i);
  }


	MatrixXd Tc = MatrixXd(n_x_, 2);
	Tc.fill(0.0);
	// Calculate Cross Correlation matrix
  MatrixXd S = MatrixXd(2,2);
	S.fill(0.0);
	// Calculate Covariance Matrix for Measurement Model
  for (int i = 0; i < 2 * n_aug_ + 1; i++){

		// residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
		z_diff(1) = normalizeAngle(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();

		// state difference
		VectorXd x_dif = Xsig_pred_.col(i) - x_;
		x_dif(3) = normalizeAngle(x_dif(3));

		Tc += weights_(i) * x_dif * z_diff.transpose();
  }

  MatrixXd R_las_ = MatrixXd(2, 2); // does not compile if moved to global.. weird
  // Add noise to our Covariance Matrix
  R_las_ <<    pow(std_laspx_, 2), 0,
          0, pow(std_laspy_, 2);
  S += R_las_;


  // Calculate Kalman Gain
  MatrixXd K = Tc * S.inverse();

  // Use measurement
  VectorXd z = VectorXd(2);
	z = meas_package.raw_measurements_;
  // residual
  VectorXd z_diff = z - z_pred;
	z_diff(1) = normalizeAngle(z_diff(1));

  x_ += K * z_diff;
	x_(3) = normalizeAngle(x_(3));

  P_ -= K*S*K.transpose();

  //std::cout << " P_ \n" << P_ << std::endl;
  //std::cout << " x_ \n" << x_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  MatrixXd Zsig = MatrixXd(3, 2* n_aug_ + 1);
  Zsig.fill(0.0);
  // Create Radar Measurement Model Sigma Points
  for (int i =0; i < 2 * n_aug_ + 1; i++){
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double p_v = Xsig_pred_(2,i);
    double p_yaw = Xsig_pred_(3,i);
    double p_yawd = Xsig_pred_(4,i);

    double v_x = cos(p_yaw)*p_v;
    double v_y = sin(p_yaw)*p_v;

    double r =  sqrt(pow(p_x, 2) + pow(p_y, 2));
    // check for /0
		if (fabs(r) < 0.000001) r = 0.000001;
    Zsig(0,i) = r;
    Zsig(1,i) = atan2(p_y, p_x);
    Zsig(2,i) = (p_x*v_x + p_y*v_y)/r;
  }

  // Calculate Mean for Measurement Model
  VectorXd z_pred = VectorXd(3);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    z_pred += weights_(i) * Zsig.col(i);
  }

	MatrixXd Tc = MatrixXd(n_x_, 3);
	Tc.fill(0.0);
  MatrixXd S = MatrixXd(3,3);
	S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
		z_diff(1) = normalizeAngle(z_diff(1));

		// Calculate Covariance Matrix for Measurement Model
    S += weights_(i) * z_diff * z_diff.transpose();

		VectorXd x_dif = Xsig_pred_.col(i) - x_;
		x_dif(3) = normalizeAngle(x_dif(3));

		// Calculate Cross Correlation Matrix
		Tc += weights_(i) * x_dif * z_diff.transpose();
  }

	// Add Sensor Noise
  MatrixXd R = MatrixXd(3, 3);
  // Add noise to our Covariance Matrix
  R <<  pow(std_radr_, 2), 0, 0,
        0, pow(std_radphi_, 2), 0,
        0, 0, pow(std_radrd_, 2);
  S += R;

  // Calculate Kalman Gain
  MatrixXd K = Tc * S.inverse();

  // Use Sensor Measurements
  VectorXd z = VectorXd(3);
  double ro = meas_package.raw_measurements_[0];
  double phi = meas_package.raw_measurements_[1];
  double ro_d = meas_package.raw_measurements_[2];
  //z << ro, phi, ro_d;
  z << meas_package.raw_measurements_;
	// residual
  VectorXd z_diff = z - z_pred;
  z_diff(1) = normalizeAngle(z_diff(1));

  x_ += K * z_diff;
  x_(3) = normalizeAngle(x_(3));

	// Update Covariance Matrix
  P_ -= K*S*K.transpose();

  //std::cout << " P_ \n" << P_ << std::endl;
  //std::cout << " x_ \n" << x_ << std::endl;

}

double UKF::normalizeAngle(double phi){
  return  atan2(sin(phi), cos(phi));
}