#include "kalman_filter.h"
#include <iostream>

#define PI 3.141592653589793238462643383279502884197169399375105820974
#define TWO_PI 6.283185307179586476925286766559005768394338798750211641949

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;


KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
	MatrixXd Ft = F_.transpose();
	x_ = F_ * x_;
	P_ = F_ * P_ * Ft + Q_;
	//cout << "   KalmanFilter::Predict::x_ = " << x_.transpose() << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
	//cout << "   KalmanFilter::Update::x_ = " << x_.transpose() << endl;
	//cout << endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	MatrixXd I3 = MatrixXd::Identity(3, 3);

	float rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
	float phi = atan2(x_(1), x_(0));
	float rho_dot;
	if (fabs(rho) < 0.0001) {
		rho_dot = 0;
	}
	else {
		rho_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;
	}

	VectorXd z_pred(3);
	z_pred << rho, phi, rho_dot;
	VectorXd y = z - z_pred;

	// Difference of angle shall be in the range (-pi, pi)
	if (y(1) > PI) {
		y(1) -= TWO_PI;
	}
	else if (y(1) < -PI) {
		y(1) += TWO_PI;
	}
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_ + I3*0.000001;
	MatrixXd St = S.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd PHtt = PHt.transpose();
	MatrixXd Kt = St.colPivHouseholderQr().solve(PHtt);
	MatrixXd K = Kt.transpose();

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
