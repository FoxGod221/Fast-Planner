/**
* This file is part of Fast-Planner.
*
* Copyright 2019 Boyu Zhou, Aerial Robotics Group, Hong Kong University of Science and Technology, <uav.ust.hk>
* Developed by Boyu Zhou <bzhouai at connect dot ust dot hk>, <uv dot boyuzhou at gmail dot com>
* for more information see <https://github.com/HKUST-Aerial-Robotics/Fast-Planner>.
* If you use this code, please cite the respective publications as
* listed on the above website.
*
* Fast-Planner is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Fast-Planner is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Fast-Planner. If not, see <http://www.gnu.org/licenses/>.
*/



#ifndef _NON_UNIFORM_BSPLINE_H_
#define _NON_UNIFORM_BSPLINE_H_

#include <Eigen/Eigen>
#include <algorithm>
#include <iostream>

using namespace std;

namespace fast_planner {
// An implementation of non-uniform B-spline with different dimensions
// It also represents uniform B-spline which is a special case of non-uniform
class NonUniformBspline {     //非均匀B样条类
private:
  // control points for B-spline with different dimensions.
  // Each row represents one single control point
  // The dimension is determined by column number
  // e.g. B-spline with N points in 3D space -> Nx3 matrix
  Eigen::MatrixXd control_points_;      //控制点的行数表示控制点的数量，列数表示空间的维度（例如，对于三维空间，列数为3）
  /*
  首尾各有 p+1 个重复的节点，这样的设置保证了曲线始终经过第一个和最后一个控制点。 也因此完整B样条的轨迹对应的时间区间是(tpb,t(m-pb)) 而四个控制点之间的时间区间划分为m段
 
  前p_+1个节点被称为开头节点或者起始节点，他们决定了B样条的起始点，这样定义出来的起始节点为小于0的区间
  中间的p_ + 1 到 m_ - p_ 之间节点为中间节点，由累加而成。确保连续
  后m_ - p_ 个节点称为尾节点（trailing nodes）或结束节点，它们决定了 B 样条曲线的结束点，初始化时，这些节点的值与前一个节点相等，确保在这一段区间内的节点是连续的。这个区间的长度是 (m_ - p_) 
  自己写出来很好理解
  */
  int             p_, n_, m_;  // p degree, n+1 control points, m = n+p+1  p表示b样条的次数 n表示b样条的阶数 m节点数为阶数+次数
  Eigen::VectorXd u_;          // knots vector  这是B样条的节点向量，定义了每个控制点对曲线的影响范围和顺序。
  double          interval_;   // knot span \delta t  这是节点向量的跨度

  Eigen::MatrixXd getDerivativeControlPoints();   //该方法用于计算并返回B样条的导数控制点 如速度、加速度

  double limit_vel_, limit_acc_, limit_ratio_;  // physical limits and time adjustment ratio 变量表示物理限制 以及时间调整比例

public:
  NonUniformBspline() {}
  /// @brief 非均匀b样条构造函数
  /// @param points 控制点的矩阵，每行表示一个控制点，每列表示一个维度
  /// @param order  B样条的阶数
  /// @param interval   节点间隔
  NonUniformBspline(const Eigen::MatrixXd& points, const int& order, const double& interval);   
  ~NonUniformBspline();

  // initialize as an uniform B-spline
  void setUniformBspline(const Eigen::MatrixXd& points, const int& order, const double& interval);

  // get / set basic bspline info

  void                                   setKnot(const Eigen::VectorXd& knot);
  Eigen::VectorXd                        getKnot();
  Eigen::MatrixXd                        getControlPoint();
  double                                 getInterval();
  void                                   getTimeSpan(double& um, double& um_p);
  pair<Eigen::VectorXd, Eigen::VectorXd> getHeadTailPts();

  // compute position / derivative

  Eigen::VectorXd   evaluateDeBoor(const double& u);   // use u \in [up, u_mp]
  Eigen::VectorXd   evaluateDeBoorT(const double& t);  // use t \in [0, duration]
  NonUniformBspline getDerivative();

  // 3D B-spline interpolation of points in point_set, with boundary vel&acc
  // constraints
  // input : (K+2) points with boundary vel/acc; ts
  // output: (K+6) control_pts
  static void parameterizeToBspline(const double& ts, const vector<Eigen::Vector3d>& point_set,
                                    const vector<Eigen::Vector3d>& start_end_derivative,
                                    Eigen::MatrixXd&               ctrl_pts);

  /* check feasibility, adjust time */

  void   setPhysicalLimits(const double& vel, const double& acc);
  bool   checkFeasibility(bool show = false);
  double checkRatio();
  void   lengthenTime(const double& ratio);
  bool   reallocateTime(bool show = false);

  /* for performance evaluation */

  double getTimeSum();
  double getLength(const double& res = 0.01);
  double getJerk();
  void   getMeanAndMaxVel(double& mean_v, double& max_v);
  void   getMeanAndMaxAcc(double& mean_a, double& max_a);

  void recomputeInit();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // namespace fast_planner
#endif