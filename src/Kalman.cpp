#include "Kalman.hpp"

Kalman::Kalman(
      const int32_t fHeight, const int32_t fWidth,
      const int32_t hHeight, const int32_t hWidth,
      const int32_t qHeight, const int32_t qWidth,
      const int32_t rHeight, const int32_t rWidth,
      const int32_t iHeight, const int32_t iWidth,
      const int32_t xLength,
      const int32_t pHeight, const int32_t pWidth)
   : m_F(fHeight,fWidth)
   , m_H(hHeight,hWidth)
   , m_Q(qHeight,qWidth)
   , m_R(rHeight,rWidth)
   , m_I(Eigen::MatrixXd::Identity(iHeight,iWidth))
   , m_x(xLength)
   , m_P(pHeight, pWidth)
   , m_previousTimestamp(0)
   , x(m_x)
   , P(m_P)

{
}

int64_t Kalman::setNewTime(const int64_t timestamp)
{
   //compute the time elapsed between the current and previous measurements
   int64_t dt = timestamp - m_previousTimestamp;
   m_previousTimestamp = timestamp;
   return dt;
}

void Kalman::filter(const Eigen::VectorXd &meas)
{
   predict();
   std::function<const Eigen::VectorXd &(const Eigen::VectorXd &, const Eigen::VectorXd &)>f_h = [this](const Eigen::VectorXd & z, const Eigen::VectorXd & x_pred)
   {
      const Eigen::VectorXd z_pred = this->m_H * x_pred;
      Eigen::VectorXd y = z - z_pred;
      return y;
   };
   update(meas, f_h, m_H, m_R);
}

void Kalman::predict()
{
   // prediction
   m_x_pred = m_F * m_x;
   m_P_pred = m_F * m_P * m_F.transpose() + m_Q;
}

void Kalman::update(const Eigen::VectorXd &meas, std::function<const Eigen::VectorXd (const Eigen::VectorXd &, const Eigen::VectorXd &)>f_h, const Eigen::MatrixXd &H, const Eigen::MatrixXd &R)
{
   const Eigen::VectorXd & z = meas; // m_H * meas;
   //const Eigen::VectorXd y = z - m_H * m_x_pred;
   const Eigen::VectorXd & y = f_h(z, m_x_pred);
   const Eigen::MatrixXd S = H * m_P_pred * H.transpose() + R;
   const Eigen::MatrixXd K = m_P_pred * H.transpose() * S.inverse();
   m_x = m_x_pred + K * y;
   m_P = (m_I - K * H) * m_P_pred;
}


