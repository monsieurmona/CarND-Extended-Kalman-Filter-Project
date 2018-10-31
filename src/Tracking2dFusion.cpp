#include "Tracking2dFusion.hpp"

#include <iostream>

Tracking2dFusion::Tracking2dFusion()
   : Kalman(4,4, // F size - Motion Model
            2,4, // H size
            4,4, // Q
            2,2, // R
            4,4, // I - Identity Matrix
            4,   // x
            4, 4 // P
            )
   , m_noise_ax(9.0)
   , m_noise_ay(9.0)
   , m_R_lidar(2,2)
   , m_R_radar(3,3)
{
   // motion model
   m_F << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

   // initialize measurement space transformation matrix
   m_H << 1, 0, 0, 0,
          0, 1, 0, 0;

   //measurement covariance matrix - laser
   m_R_lidar << 0.0225, 0,
                0, 0.0225;

   //measurement covariance matrix - radar
   m_R_radar << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;

   //state covariance matrix P
   m_P << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;

   // converts predicted state to predicted radar measurement
   // (cartesian coordinates to polar coordinates)
   // and substructs it from radar measurement
   m_f_h_radar = [](const Eigen::VectorXd & z, const Eigen::VectorXd & x_pred)
   {
      Eigen::VectorXd z_pred(3);
      z_pred.setZero();
      const double dist = sqrt(pow(x_pred(0), 2) + pow(x_pred(1), 2));

      if (fabs(dist) > 0.0001)
      {
         // see also
         // https://stackoverflow.com/questions/639283/how-should-i-tackle-arctan-and-prevent-deviding-by-zero-without-the-convenienc
         const double phi = atan2(x_pred(1), x_pred(0));

         z_pred <<
               dist,
               phi,
               (x_pred(0) * x_pred(2) + x_pred(1) * x_pred(3)) / dist;
      }

      Eigen::VectorXd y = z - z_pred;

      // normalize phi to -/+ pi
      y(1) = y(1) - TWO_PI * floor((y(1) + M_PI) / TWO_PI);

      return y;
   };

   // converts predicted state to predicted lidar measurement
   // and substructs it from lidar measurement
   m_f_h_lidar = [this](const Eigen::VectorXd & z, const Eigen::VectorXd & x_pred)
   {
      const Eigen::VectorXd z_pred = this->m_H * x_pred;
      Eigen::VectorXd y = z - z_pred;
      return y;
   };
}
void Tracking2dFusion::init(const Eigen::VectorXd &_x, const int64_t timestamp)
{
   m_x << _x(0), _x(1), 0, 0;
   setNewTime(timestamp);
}

inline Eigen::MatrixXd Tracking2dFusion::calculateJacobianRadar(const Eigen::VectorXd& x_state) const {

    Eigen::MatrixXd Hj(3,4);
    Hj.setZero();

    //recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    //pre-compute a set of terms to avoid repeated calculation
    double c1 = px*px+py*py;
    double c2 = sqrt(c1);
    double c3 = (c1*c2);

    //check division by zero
    if(fabs(c1) < 0.0001){
        std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
        return Hj;
    }

    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
            -(py/c1), (px/c1), 0, 0,
            py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj;
}

void Tracking2dFusion::track(const Eigen::VectorXd &meas, const MeasurementPackage::SensorType sensorType, const int64_t timestamp)
{
   const int64_t dt_us = setNewTime(timestamp);

   if (dt_us > 0)
   {
      // predic state, if some time past since last prediction
      const double dt_s = dt_us / 1000000.0;	//dt - expressed in seconds

      m_F(0,2) = dt_s;
      m_F(1,3) = dt_s;

      const double dt4by4 = pow(dt_s, 4) / 4.0;
      const double dt3by2 = pow(dt_s, 3) / 2.0;
      const double dt2    = pow(dt_s, 2);

      // Process Noise - Covariance
      m_Q << dt4by4 * m_noise_ax, 0,                   dt3by2 * m_noise_ax, 0,
             0,                   dt4by4 * m_noise_ay, 0,                   dt3by2 * m_noise_ay,
             dt3by2 * m_noise_ax, 0,                   dt2 * m_noise_ax,    0,
             0,                   dt3by2 * m_noise_ay, 0,                   dt2 * m_noise_ay;

      predict();
   }

   switch (sensorType) {
   case MeasurementPackage::LASER:
      updateWithLaserMeasurement(meas);
      break;
   case MeasurementPackage::RADAR:
      updateWithRadarMeasurement(meas);
      break;
   }
}

inline void Tracking2dFusion::updateWithLaserMeasurement(const Eigen::VectorXd &meas)
{
   update(meas, m_f_h_lidar, m_H, m_R_lidar);
}

inline void Tracking2dFusion::updateWithRadarMeasurement(const Eigen::VectorXd &meas)
{
   update(meas, m_f_h_radar, calculateJacobianRadar(m_x_pred), m_R_radar);
}

