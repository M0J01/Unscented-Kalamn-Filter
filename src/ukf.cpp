#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>
#include <time.h>


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
  // Will be further initialized during first measurement
  x_ = VectorXd(5);
  // initial covariance matrix
  // Initialized to values of 1
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;



  // --------- Measurement Noises. Set by sensor MFG ----------------
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

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
	weights_ = VectorXd(n_aug_*2+1);

	// set weights
	double weight_0 = lambda_/(lambda_*1.0 + n_aug_);
	weights_(0) = weight_0;
	for (int i =1; i < 2*n_aug_+1; i++){
		double weight = 0.5/(n_aug_+lambda_);
		weights_(i) = weight;
	}

	// px, py, v, yaw, yawd
	// 1, 1, 1, 1, .1 = ~1 rmse for most catagories
	// .1, .1, .1, .1, .1 = ~ decent RMSE until it threads the middle loop. Could be due to the yawdd

	P_ <<   .15, 0, 0, 0, 0,
					0, .15, 0, 0, 0,
					0, 0, 5, 0, 0,
					0, 0, 0, .1, 0,
					0, 0, 0, 0, .1;

	/*
	P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
					-0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
					0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
					-0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
					-0.0020,    0.0060,    0.0008,    0.0100,    0.0123;
	*/
	R_las_ = MatrixXd(2,2);
	R_las_ << 	std_laspx_ * std_laspx_, 0,
							0, std_laspy_ * std_laspy_;

	R_radar_ = MatrixXd(3,3);
	R_radar_ << std_radr_*std_radr_, 0, 0,
							0, std_radphi_*std_radphi_, 0,
							0, 0, std_radrd_*std_radrd_;

	NIS_las_;
	NIS_radar_;

	ofstream myfile_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


	//time_t time1;
	//time(&time1);

	// If Not Initialized
  if (!is_initialized_){
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
			double ro = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);
			//double ro_dot = meas_package.raw_measurements_(2);
			x_ << cos(phi)*ro, sin(phi)*ro, 2.2049, 0, 0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
			double px = meas_package.raw_measurements_(0);
			double py = meas_package.raw_measurements_(1);
			x_ << px, py, 2.2049, 0.5015, 0.3528;
		}

		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
		std::cout << "Initialized\n";

		//time_t time2;
		//time(&time2);
		//std::cout <<"Initialization took  " << (time2 - time1) << " Seconds \n";

	}

	// If Is Initialized
	else if (is_initialized_){

		//Predict
		double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
		time_us_ = meas_package.timestamp_;
		Prediction(delta_t);


		//std::cout<< "delta_t = "<< delta_t << std::endl;

		//time_t time3;
		//time(&time3);
		//std::cout <<"Prediction took  " << (time3 - time1) << " Seconds \n";


		// Update
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
			// Update Radar
			UpdateRadar(meas_package);
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
			// Update Lidar
			UpdateLidar(meas_package);
		}
	}

}

double NormalizeAngle(double  phi)
{
	return  atan2(sin(phi), cos(phi));
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {


	// ---------------------------------------------------------- Generate Sigma Points

	//time_t time1;
	//time(&time1);

	// Set up our augmentation containers
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	MatrixXd P_aug = MatrixXd(7,7);
	VectorXd x_aug = VectorXd(7);

	// create augmented mean state
	x_aug(0) = x_(0);
	x_aug(1) = x_(1);
	x_aug(2) = x_(2);
	x_aug(3) = x_(3);
	x_aug(4) = x_(4);
	//x_aug.head(5) = x_();

	x_aug(5) = 0;
	x_aug(6) = 0;

	// Create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5,5) = P_;
	P_aug(5,5) = std_a_*std_a_;
	P_aug(6,6) = std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	// create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i< n_aug_; i++){
		Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_)*L.col(i);
		Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*L.col(i);
	}

	//time_t time2;
	//time(&time2);
	//std::cout <<"Generating Sigma Points took  " << (time2 - time1) << " Seconds \n";




	// ----------------------------------------------------------  Predict Sigma Points

	for (int i = 0; i < 2*n_aug_ + 1; i++) {

		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		double px_p, py_p;

			// Calculate Sigma Points

			// Avoid divide by 0
			if(fabs(yawd) > 0.0001){
				px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
				py_p = p_y + v/yawd * (-cos(yaw + yawd*delta_t) + cos(yaw));
			}
			else{
				px_p = p_x + v*delta_t*cos(yaw);
				py_p = p_y + v*delta_t*sin(yaw);
			}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;
		// Add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p += 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p += nu_yawdd*delta_t;

		// Store Sigma Points
		Xsig_pred_(0,i) = px_p;
		Xsig_pred_(1,i) = py_p;
		Xsig_pred_(2,i) = v_p;
		Xsig_pred_(3,i) = yaw_p;
		Xsig_pred_(4,i) = yawd_p;

	}

	//time_t time3;
	//time(&time3);
	//std::cout <<"Predicting Sigma Points took  " << (time3 - time2) << " Seconds \n";
	
	// ----------------------------------------------------------  Predict Mean and Covariance


	VectorXd x_mean = VectorXd(n_x_);
	//MatrixXd P_predicted = MatrixXd(n_x_, n_x_);
	// predicted state mean

	x_mean.fill(0.0);
	for (int i = 0; i < 2*n_aug_ + 1; i++){
		// double weight = weights_(i);
		x_mean += weights_(i)*Xsig_pred_.col(i);
	}

	//time_t timem;
	//time(&timem);
	//std::cout <<"Filling the Mean took  " << (timem - time3) << " Seconds \n";

	P_.fill(0.0);
	for (int i = 1; i < 2*n_aug_ + 1; i++) {
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - Xsig_pred_.col(0);
		//normalize angles
		//std::cout << "x_mean during p_ calc = "<< x_mean(3) << std::endl;
		//std::cout << "x_diff before normalization = "<< x_diff(3) << std::endl;

		x_diff(3) = NormalizeAngle(x_diff(3));
		//while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;		// This takes a lot of time for some reason.
		//while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;		// 4 billion... okay, yeah can take some time.
		// store P_predicted
		P_ += weights_(i) * x_diff * x_diff.transpose();

		//std::cout << "x_diff after normalization = "<< x_diff(3) << std::endl;
	}

/*
	for (int i = 0; i < 2*n_aug_ + 1; i++) {
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_mean;
		//normalize angles
		std::cout << "x_mean during p_ calc = "<< x_mean(3) << std::endl;
		std::cout << "x_diff before normalization = "<< x_diff(3) << std::endl;
		while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;		// This takes a lot of time for some reason.
		while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;		// 4 billion... okay, yeah can take some time.
		// store P_predicted
		P_ += weights_(i) * x_diff * x_diff.transpose();

		std::cout << "x_diff after normalization = "<< x_diff(3) << std::endl;
	}
*/
	//time_t timepmcv;
	//time(&timepmcv);
	//std::cout <<"subtracting the Mean took  " << (timepmcv - timem) << " Seconds \n";
	//std::cout <<"Predicting Mean and Covariance took  " << (timepmcv - time3) << " Seconds \n";

}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{

	//time_t timeUL1;
	//time(&timeUL1);


	// Radar Space State
	int n_z = 2;
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	// Radar Measurement Values
	VectorXd z = VectorXd(n_z);
	double px = meas_package.raw_measurements_(0);
	double py = meas_package.raw_measurements_(1);
	z << px, py;

	// Calculate Sigma Points in Measurement Space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		// grab our predicted sigma values
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw) * v;  // v_x
		double v2 = sin(yaw) * v; // v_y

		// measurement model
		Zsig(0, i) = p_x;
		Zsig(1, i) = p_y;
	}

	//time_t timeUL2;
	//time(&timeUL2);
	//std::cout <<"Calculating the LU Sigma points took  " << (timeUL2 - timeUL1) << " Seconds \n";

	// Calculate Mean Predicted Measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}


	//time_t timeUL3;
	//time(&timeUL3);
	//std::cout <<"Calculating the LU Mean predicted Measurement took  " << (timeUL3 - timeUL2) << " Seconds \n";

	// Create Covariance Matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	// Create Cross Correlation Matrix Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);


	// Calculate S and Tc
	for (int i = 1; i < 2 * n_aug_ + 1; i++) {

		// Calculate Measurement Covariance Matrix S
		// residual
		VectorXd z_diff = Zsig.col(i) - Zsig.col(0);
		// angle norm
		//z_diff(1) = NormalizeAngle(z_diff(1));
		//while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
		//while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
		S += weights_(i) * z_diff * z_diff.transpose();

		// Calculate Cross Correlation Matrix Tc
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - Xsig_pred_.col(0);
		// angle norm
		x_diff(3) = NormalizeAngle(x_diff(3));
		//while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
		//while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

	}

	//time_t timeUL4;
	//time(&timeUL4);
	//std::cout <<"Calculating the LU S and TC took  " << (timeUL4 - timeUL3) << " Seconds \n";

	// add measurement noise covariance matrix

	S = S + R_las_;


	// Calculate Kalman Gain
	MatrixXd K = Tc * S.inverse();

	// Calculate residual
	VectorXd z_diff = z - z_pred;
	// Normalize residual
	//z_diff(1) = NormalizeAngle(z_diff(1));
	//while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
	//while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;


	// Update state and covariance matrix
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();
	NIS_las_ = (z_pred - z).transpose() * S.inverse() * (z_pred - z);

	std::cout << "NIS Laser : " << NIS_las_ << std::endl;




	//time_t timeUL5;
	//time(&timeUL5);
	//std::cout <<"Calculating the LU x and P took  " << (timeUL5 - timeUL4) << " Seconds \n";
}



/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

	//time_t timeUR1;
	//time(&timeUR1);

	// Radar Space State
	int n_z = 3;
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	// Radar Measurement Values
	VectorXd z = VectorXd(n_z);
	double ro = meas_package.raw_measurements_(0);
	double phi = meas_package.raw_measurements_(1);
	double ro_dot = meas_package.raw_measurements_(2);
	z << ro, phi, ro_dot;

	// Calculate Sigma Points in Measurement Space
	for (int i =0; i < 2*n_aug_ + 1; i++){
		// grab our predicted sigma values
		double p_x = Xsig_pred_(0,i);
		double p_y = Xsig_pred_(1,i);
		double v = Xsig_pred_(2,i);
		double yaw = Xsig_pred_(3,i);

		//double v1 = sin(yaw)*v;	// v_x
		//double v2 = cos(yaw)*v; // v_y
		double v1 = cos(yaw)*v;	// v_x
		double v2 = sin(yaw)*v; // v_y

		// measurement model
		Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);

		if (abs(p_x)>1e-5){ Zsig(1,i) = atan2(p_y,p_x);}
		else { Zsig(1,i) = M_PI/2;}

		if (fabs(Zsig(0,i)) < 1e-5){
			Zsig(2,i) = 0.0;
		} else {
			Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
			//Zsig(2,i) = v*(p_x*cos(yaw)+p_y*sin(yaw))/Zsig(0,i);
		}
	}

  //time_t timeUR2;
	//time(&timeUR2);
	//std::cout <<"Calculating the UR Sigma points took  " << (timeUR2 - timeUR1) << " Seconds \n";

	// Calculate Mean Predicted Measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i =0; i < 2*n_aug_+1; i++){
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//time_t timeUR3;
	//time(&timeUR3);
	//std::cout <<"Calculating the UR Mean took  " << (timeUR3 - timeUR2) << " Seconds \n";

	// Create Covariance Matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	// Create Cross Correlation Matrix Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);


	// Calculate S and Tc
	for (int i = 1; i < 2*n_aug_ + 1; i++){

		// Calculate Measurement Covariance Matrix S
		// residual
		VectorXd z_diff = Zsig.col(i) - Zsig.col(0);
		// angle norm
		z_diff(1) = NormalizeAngle(z_diff(1));
		//while ( z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
		//while ( z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
		S += weights_(i)*z_diff*z_diff.transpose();

		// Calculate Cross Correlation Matrix Tc
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - Xsig_pred_.col(0);
		// angle norm

		x_diff(3) = NormalizeAngle(x_diff(3));
		//while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
		//while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

	}


	//time_t timeUR4;
	//time(&timeUR4);
	//std::cout <<"Calculating the UR S and Tc took  " << (timeUR4 - timeUR3) << " Seconds \n";

	// add measurement noise covariance matrix

	S = S + R_radar_;


	// Calculate Kalman Gain
	MatrixXd K = Tc * S.inverse();

	// Calculate residual
	VectorXd z_diff = z - z_pred;
	// Normalize residual
	z_diff(1) = NormalizeAngle(z_diff(1));
	//while ( z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
	//while ( z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;


	// Update state and covariance matrix
	x_ +=  K*z_diff;
	P_ -= K*S*K.transpose();

	NIS_radar_ = (z_pred - z).transpose() * S.inverse() * (z_pred - z);

	std::cout << "NIS Radar : " << NIS_radar_ << std::endl;

	//time_t timeUR5;
	//time(&timeUR5);
	//std::cout <<"Calculating the UR x and P took  " << (timeUR5 - timeUR4) << " Seconds \n";

}


void UKF::write_file_data(ofstream &myfile_){

	double v1 = cos(x_(3))*x_(2);
	double v2 = sin(x_(3))*x_(2);

	myfile_ << v1 << ", ";
	myfile_ << v2 << ", ";
	myfile_ << NIS_radar_ << ", ";
	myfile_ << NIS_las_ << ", ";
	myfile_ << x_(0) << ", ";
	myfile_ << x_(1) << ", ";
	myfile_ << x_(2) << ", ";
	myfile_ << x_(3) << ", ";
	myfile_ << x_(4); //<< ", ";

	myfile_ << " \r\n";
	//myfile_ << std::endl;


}