#ifndef KALMAN_HPP
#define KALMAN_HPP

#include "Eigen/Dense"

class Kalman
{
public:
   void filter(const Eigen::VectorXd &meas);
   int64_t setNewTime(const int64_t timestamp);

protected:
   Kalman(const int32_t fHeight, const int32_t fWidth,
          const int32_t hHeight, const int32_t hWidth,
          const int32_t qHeight, const int32_t qWidth,
          const int32_t rHeight, const int32_t rWidth,
          const int32_t iHeight, const int32_t iWidth,
          const int32_t xLength,
          const int32_t pHeight, const int32_t pWidth);

   // predict step
   void predict();

   // update step
   void update(const Eigen::VectorXd &meas, std::function<const Eigen::VectorXd (const Eigen::VectorXd &, const Eigen::VectorXd &)>f_h, const Eigen::MatrixXd &H, const Eigen::MatrixXd &R);

   // Motion Model
   Eigen::MatrixXd m_F;

   // Measurement Space Transformation Matrix
   Eigen::MatrixXd m_H;

   // Process Noise - Covariance
   Eigen::MatrixXd m_Q;

   // Measurement Noise - Covariance
   Eigen::MatrixXd m_R;

   // Identity
   Eigen::MatrixXd m_I;

   // object state
   Eigen::VectorXd m_x;

   // object covariance matrix
   Eigen::MatrixXd m_P;

   int64_t m_previousTimestamp;

   // predicted state
   Eigen::VectorXd m_x_pred;

public:
   // Read only interface
   const Eigen::VectorXd & x;
   const Eigen::MatrixXd & P;
   static constexpr double TWO_PI = 2 * M_PI;

private:
   // predicted process noise
   Eigen::MatrixXd m_P_pred;
};

#endif // KALMAN_HPP
