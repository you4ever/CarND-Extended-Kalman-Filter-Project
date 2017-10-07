#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#define EPS 0.0001 // A very small number

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  Hj_ << 0, 0, 0, 0,
	  0, 0, 0, 0,
	  0, 0, 0, 0;
  
  //efk_.Intl();

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


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
    //cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	  double rho = measurement_pack.raw_measurements_[0];
	  double phi = measurement_pack.raw_measurements_[1];
	  double rho_dot = measurement_pack.raw_measurements_[2];
	  double x = rho * cos(phi);
	  double y = rho * sin(phi);
	  double x_dot = rho_dot * sin(phi);
	  double y_dot = rho_dot* cos(phi);
	  ekf_.x_ << x, y, x_dot, y_dot;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	  double x = measurement_pack.raw_measurements_[0];
	  double y = measurement_pack.raw_measurements_[1];
	  ekf_.x_ << x, y, 0, 0;
    }

	if (fabs(ekf_.x_(0)) < EPS && fabs(ekf_.x_(1)) < EPS) {
		ekf_.x_(0) = EPS;
		ekf_.x_(1) = EPS;
	}

	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1000, 0,
		0, 0, 0, 1000;

	previous_timestamp_ = measurement_pack.timestamp_;
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
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  double dt_2 = dt * dt;
  double dt_3 = dt_2*dt;
  double dt_4 = dt_2*dt_2;
  double noise_ax = 9.0;
  double noise_ay = 9.0;

  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
	  0, 1, 0, dt,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
	  0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
	  dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
	  0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	  //cout << "   Radar update " << endl;

    // Radar updates
	  Tools tools;
	  Hj_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.H_ = Hj_;
	  ekf_.R_ = R_radar_;
	  VectorXd z = measurement_pack.raw_measurements_;	
	  ekf_.UpdateEKF(z);
  } else {
	  //cout << "   Laser update " << endl;
    // Laser updates
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << "xNew_ = " << ekf_.x_ << endl;
  //cout << "PNew_ = " << ekf_.P_ << endl;
  //cout << "end update " << endl;

}
