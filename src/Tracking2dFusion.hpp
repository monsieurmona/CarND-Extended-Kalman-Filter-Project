#ifndef TRACKING2DFUSION_HPP
#define TRACKING2DFUSION_HPP

#include "Kalman.hpp"
#include "measurement_package.h"

class Tracking2dFusion : public Kalman
{
public:
   Tracking2dFusion();
   void init(const Eigen::VectorXd &_x, const int64_t timestamp);

   void track(const Eigen::VectorXd &meas, const MeasurementPackage::SensorType sensorType, const int64_t timestamp);

private:
   std::function<const Eigen::VectorXd(const Eigen::VectorXd & z, const Eigen::VectorXd &)> m_f_h_radar;
   std::function<const Eigen::VectorXd(const Eigen::VectorXd & z, const Eigen::VectorXd &)> m_f_h_lidar;

   inline Eigen::MatrixXd calculateJacobianRadar(const Eigen::VectorXd& x_state) const;
   inline void updateWithLaserMeasurement(const Eigen::VectorXd &meas);
   inline void updateWithRadarMeasurement(const Eigen::VectorXd &meas);

   //acceleration noise components
   double m_noise_ax;
   double m_noise_ay;

   // measurement noise - covariance matrix
   Eigen::MatrixXd m_R_lidar;
   Eigen::MatrixXd m_R_radar;
};

#endif // TRACKING2DFUSION_HPP
