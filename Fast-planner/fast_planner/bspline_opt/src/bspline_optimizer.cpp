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



#include "bspline_opt/bspline_optimizer.h"
#include <nlopt.hpp>
// using namespace std;

namespace fast_planner {
//用来表示当前阶段
const int BsplineOptimizer::SMOOTHNESS  = (1 << 0); //1 平滑度
const int BsplineOptimizer::DISTANCE    = (1 << 1); //2 距离
const int BsplineOptimizer::FEASIBILITY = (1 << 2); //4 可行性
const int BsplineOptimizer::ENDPOINT    = (1 << 3); //8 端点
const int BsplineOptimizer::GUIDE       = (1 << 4); //16  引导优化
const int BsplineOptimizer::WAYPOINTS   = (1 << 6); //64 关键点
//表示在引导阶段需要同时优化平滑度和引导。用位运算符 |（按位或）将两个常量合并，结果是 1 | 16 = 17，即 GUIDE_PHASE 的值为 17。
const int BsplineOptimizer::GUIDE_PHASE = BsplineOptimizer::SMOOTHNESS | BsplineOptimizer::GUIDE;
//表示在正常阶段需要同时优化平滑度、距离和可行性。用 | 操作符将三个常量合并 距离指碰撞优化 可行性指动力学约束
const int BsplineOptimizer::NORMAL_PHASE =
    BsplineOptimizer::SMOOTHNESS | BsplineOptimizer::DISTANCE | BsplineOptimizer::FEASIBILITY;

//以下各个参数意义在.h中有说明 主要是为计算代价函数的时候使用的权重系数
//从 ROS 参数服务器获取名为 optimization/lambda1 的参数，并将其存储到变量中。如果该参数没有被设置或者无法获取，那么就会使用默认值 -1.0。
void BsplineOptimizer::setParam(ros::NodeHandle& nh) {
  nh.param("optimization/lambda1", lambda1_, -1.0);
  nh.param("optimization/lambda2", lambda2_, -1.0);
  nh.param("optimization/lambda3", lambda3_, -1.0);
  nh.param("optimization/lambda4", lambda4_, -1.0);
  nh.param("optimization/lambda5", lambda5_, -1.0);
  nh.param("optimization/lambda6", lambda6_, -1.0);
  nh.param("optimization/lambda7", lambda7_, -1.0);
  nh.param("optimization/lambda8", lambda8_, -1.0);

  nh.param("optimization/dist0", dist0_, -1.0);
  nh.param("optimization/max_vel", max_vel_, -1.0);
  nh.param("optimization/max_acc", max_acc_, -1.0);
  nh.param("optimization/visib_min", visib_min_, -1.0);
  nh.param("optimization/dlmin", dlmin_, -1.0);
  nh.param("optimization/wnl", wnl_, -1.0);

  nh.param("optimization/max_iteration_num1", max_iteration_num_[0], -1);
  nh.param("optimization/max_iteration_num2", max_iteration_num_[1], -1);
  nh.param("optimization/max_iteration_num3", max_iteration_num_[2], -1);
  nh.param("optimization/max_iteration_num4", max_iteration_num_[3], -1);
  nh.param("optimization/max_iteration_time1", max_iteration_time_[0], -1.0);
  nh.param("optimization/max_iteration_time2", max_iteration_time_[1], -1.0);
  nh.param("optimization/max_iteration_time3", max_iteration_time_[2], -1.0);
  nh.param("optimization/max_iteration_time4", max_iteration_time_[3], -1.0);

  nh.param("optimization/algorithm1", algorithm1_, -1);
  nh.param("optimization/algorithm2", algorithm2_, -1);
  nh.param("optimization/order", order_, -1);
}
//下面是一些参数设置
void BsplineOptimizer::setEnvironment(const EDTEnvironment::Ptr& env) {
  this->edt_environment_ = env;
}

void BsplineOptimizer::setControlPoints(const Eigen::MatrixXd& points) {
  control_points_ = points;
  dim_            = control_points_.cols();
}

void BsplineOptimizer::setBsplineInterval(const double& ts) { bspline_interval_ = ts; }

void BsplineOptimizer::setTerminateCond(const int& max_num_id, const int& max_time_id) {
  max_num_id_  = max_num_id;      //这个变量保存最大迭代次数的限制
  max_time_id_ = max_time_id;
}
//设置目标函数并打印出当前所选的优化目标。择优化过程中使用的目标函数类型。
void BsplineOptimizer::setCostFunction(const int& cost_code) {
  cost_function_ = cost_code;

  // print optimized cost function 通过操作（&），可以检查 cost_function_ 中是否包含某个特定的目标函数类型 用来确定都需要对哪些项进行优化
  string cost_str;
  if (cost_function_ & SMOOTHNESS) cost_str += "smooth |";
  if (cost_function_ & DISTANCE) cost_str += " dist  |";
  if (cost_function_ & FEASIBILITY) cost_str += " feasi |";
  if (cost_function_ & ENDPOINT) cost_str += " endpt |";
  if (cost_function_ & GUIDE) cost_str += " guide |";
  if (cost_function_ & WAYPOINTS) cost_str += " waypt |";

  ROS_INFO_STREAM("cost func: " << cost_str);
}
//引导路径  应该指的是前端A*搜索出来的的路径
void BsplineOptimizer::setGuidePath(const vector<Eigen::Vector3d>& guide_pt) { guide_pts_ = guide_pt; }

void BsplineOptimizer::setWaypoints(const vector<Eigen::Vector3d>& waypts,
                                    const vector<int>&             waypt_idx) {
  waypoints_ = waypts;
  waypt_idx_ = waypt_idx;
}
//！！！！！！！！！！！！！！！！！！！！！！！！整个优化流程
Eigen::MatrixXd BsplineOptimizer::BsplineOptimizeTraj(const Eigen::MatrixXd& points, const double& ts,
                                                      const int& cost_function, int max_num_id,
                                                      int max_time_id) {
  setControlPoints(points);   //设置控制点
  setBsplineInterval(ts);     //时间间隔  
  setCostFunction(cost_function); //优化目标函数
  setTerminateCond(max_num_id, max_time_id);  //终止条件

  optimize();       //优化
  return this->control_points_;   //返回优化后的控制点
}
//！！！！！！！！！！！！！！！！！！！！！！！！！优化实现部分
void BsplineOptimizer::optimize() {
  /* initialize solver */
  iter_num_        = 0;     //迭代次数
  min_cost_        = std::numeric_limits<double>::max();  //应该是优化到这个值就可以接受
  const int pt_num = control_points_.rows();  //控制点个数
  //修改梯度值值存放空间大小
  g_q_.resize(pt_num);
  g_smoothness_.resize(pt_num);
  g_distance_.resize(pt_num);
  g_feasibility_.resize(pt_num);
  g_endpoint_.resize(pt_num);
  g_waypoints_.resize(pt_num);
  g_guide_.resize(pt_num);

  if (cost_function_ & ENDPOINT) {
    //表示优化过程中控制变量的数量。在进行路径优化时，控制点的数量会影响变量的数量。 考虑端点约束 使用因此只减去一个阶数
    variable_num_ = dim_ * (pt_num - order_);
    // end position used for hard constraint    
    //端点为硬约束 因为必须要经过 并且要这段计算使用了一个加权平均公式，这个公式通常用于计算 B-spline 或贝塞尔曲线的端点位置。 通常用于平滑曲线的端点。
    end_pt_ = (1 / 6.0) *
        (control_points_.row(pt_num - 3) + 4 * control_points_.row(pt_num - 2) +
         control_points_.row(pt_num - 1));
  } else {
    variable_num_ = max(0, dim_ * (pt_num - 2 * order_)) ;
  }

  /* do optimization using NLopt slover */
  //定义一个优化器对象 如果 isQuadratic() 返回 true，使用 algorithm1_（通常表示一种针对二次优化问题的算法）。
  nlopt::opt opt(nlopt::algorithm(isQuadratic() ? algorithm1_ : algorithm2_), variable_num_);
  opt.set_min_objective(BsplineOptimizer::costFunction, this);    //优化的目标函数是 this 作为 func_data 传递给目标函数。
  opt.set_maxeval(max_iteration_num_[max_num_id_]);   //设置优化时的最大次数与时间
  opt.set_maxtime(max_iteration_time_[max_time_id_]);
  opt.set_xtol_rel(1e-5); //设置优化算法停止的精度条件。1e-5 表示相对精度容忍度，即如果在连续迭代之间的变量变化小于 1e-5，优化器将认为已经收敛，并停止优化。

  vector<double> q(variable_num_);    //用于存储优化过程中的控制点
  //遍历控制点的索引，从 order_ 到 pt_num。 order_ 是 B 样条的阶数，通常是样条的控制点之间的多项式次数，因此从 order_ 开始是因为 B 样条的第一个有效控制点从 order_ 索引开始
  for (int i = order_; i < pt_num; ++i) {
    if (!(cost_function_ & ENDPOINT) && i >= pt_num - order_) continue;   //如果没有使用 ENDPOINT 约束，代码将跳过控制点的最后部分，因为这些控制点可能是约束的边界点。
    for (int j = 0; j < dim_; j++) {
      //将这个坐标值赋值给优化变量 q 中相应的位置。
      q[dim_ * (i - order_) + j] = control_points_(i, j); //计算当前控制点在 q 中的对应位置。i - order_ 是因为 B 样条的前 order_ 个点不参与优化，实际优化的控制点从 order_ 开始。
    }
  }
  //这段代码的目的是设置优化变量 q 的上下边界（即变量的约束范围），并将这些边界应用到 nlopt 优化器 opt 中。
  if (dim_ != 1) {
    vector<double> lb(variable_num_), ub(variable_num_); //用来存储q中变量的上下界
    const double   bound = 10.0;
    //设置上下界
    for (int i = 0; i < variable_num_; ++i) {
      lb[i] = q[i] - bound;
      ub[i] = q[i] + bound;
    }
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
  }
  //使用 try-catch 块来捕获并处理在优化过程中可能抛出的异常。try 块中是主要的优化过程，如果发生异常（如优化过程中的错误），会跳转到 catch 块中。

  try {
    // cout << fixed << setprecision(7);
    // vec_time_.clear();
    // vec_cost_.clear();
    // time_start_ = ros::Time::now();

    double        final_cost;
    //！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    //优化器的使用 这句代码之前都是创建优化器 做优化准备 同时在下面分多个函数计算了各个项的代价值 以及梯度 然后将代价值全加起来 进行整体的优化
    nlopt::result result = opt.optimize(q, final_cost);   //是优化结果的类型，包含了优化是否成功、达到的最优值等信息。  优化对象是q 调用优化是在这一步

    /* retrieve the optimization result */
    // cout << "Min cost:" << min_cost_ << endl;
  } catch (std::exception& e) {   
    ROS_WARN("[Optimization]: nlopt exception");
    cout << e.what() << endl;
  }

  for (int i = order_; i < control_points_.rows(); ++i) {
    if (!(cost_function_ & ENDPOINT) && i >= pt_num - order_) continue; //外部循环从 order_ 开始，遍历所有控制点。order_ 表示 B-spline 的阶数，它用于计算控制点的有效范围。
    for (int j = 0; j < dim_; j++) {
      control_points_(i, j) = best_variable_[dim_ * (i - order_) + j];  //优化过程通过 nlopt 找到的最优解被存储在 best_variable_ 中，并用这些最优解来更新控制点。
    }
  }

  if (!(cost_function_ & GUIDE)) ROS_INFO_STREAM("iter num: " << iter_num_);  //如果 cost_function_ 中不包含 GUIDE，则输出当前的迭代次数 iter_num_。
}
//计算平滑项代价   采用jerk最小化,代码中采用的是差分的形式 
void BsplineOptimizer::calcSmoothnessCost(const vector<Eigen::Vector3d>& q, double& cost,
                                          vector<Eigen::Vector3d>& gradient) {
  
  cost = 0.0;
  Eigen::Vector3d zero(0, 0, 0);
  std::fill(gradient.begin(), gradient.end(), zero);
  Eigen::Vector3d jerk, temp_j;

  for (int i = 0; i < q.size() - order_; i++) {
    /* evaluate jerk */
    jerk = q[i + 3] - 3 * q[i + 2] + 3 * q[i + 1] - q[i];   //计算jerk 对这4个控制点做一个差分操作，来估算加速度的变化率
    cost += jerk.squaredNorm();     //使用jerk的平方来估计平滑项
    temp_j = 2.0 * jerk;  //是计算 jerk 的梯度的系数。
    /* jerk gradient */
    gradient[i + 0] += -temp_j; //对于不同控制点 梯度的增值不同
    gradient[i + 1] += 3.0 * temp_j;
    gradient[i + 2] += -3.0 * temp_j;
    gradient[i + 3] += temp_j;
  }
}
//计算距离代价 也就是碰撞代价
void BsplineOptimizer::calcDistanceCost(const vector<Eigen::Vector3d>& q, double& cost,
                                        vector<Eigen::Vector3d>& gradient) {
  cost = 0.0;
  Eigen::Vector3d zero(0, 0, 0);
  std::fill(gradient.begin(), gradient.end(), zero);

  double          dist;
  Eigen::Vector3d dist_grad, g_zero(0, 0, 0); //用于存储距离的梯度（即该点的梯度向量）。

  int end_idx = (cost_function_ & ENDPOINT) ? q.size() : q.size() - order_;//是否包含端点

  for (int i = order_; i < end_idx; i++) {
    edt_environment_->evaluateEDTWithGrad(q[i], -1.0, dist, dist_grad);   // 算当前控制点 q[i] 到障碍物的距离 dist 及其梯度 dist_grad
    if (dist_grad.norm() > 1e-4) dist_grad.normalize(); //判断梯度的范数是否大于一个小阈值（1e-4），如果是，则将梯度向量归一化。这样可以避免梯度过大或不稳定。

    if (dist < dist0_) {      //如果小于指定的安全距离 进行误差累计
      cost += pow(dist - dist0_, 2);
      gradient[i] += 2.0 * (dist - dist0_) * dist_grad; //根据距离的变化（dist - dist0_）和距离的梯度（dist_grad）更新控制点的梯度。
    }
  }
}
//可行性误差计算 主要是动力学约束 是否满足
void BsplineOptimizer::calcFeasibilityCost(const vector<Eigen::Vector3d>& q, double& cost,
                                           vector<Eigen::Vector3d>& gradient) {
  cost = 0.0;
  Eigen::Vector3d zero(0, 0, 0);
  std::fill(gradient.begin(), gradient.end(), zero);

  /* abbreviation */
  double ts, vm2, am2, ts_inv2, ts_inv4;
  vm2 = max_vel_ * max_vel_;
  am2 = max_acc_ * max_acc_;

  ts      = bspline_interval_;
  ts_inv2 = 1 / ts / ts;    //是与时间间隔相关的常数，分别是 ts 的倒数平方和四次方，这些常数在计算速度和加速度的约束时用到。
  ts_inv4 = ts_inv2 * ts_inv2;

  /* velocity feasibility */
  for (int i = 0; i < q.size() - 1; i++) {
    Eigen::Vector3d vi = q[i + 1] - q[i];   //用于计算后面的速度

    for (int j = 0; j < 3; j++) {
      double vd = vi(j) * vi(j) * ts_inv2 - vm2;    //vi*vi*ts_in2表示速度的平方 vd用来判断是否超过指定限制
      if (vd > 0.0) {
        cost += pow(vd, 2);   //误差累计

        double temp_v = 4.0 * vd * ts_inv2;
        gradient[i + 0](j) += -temp_v * vi(j);    //梯度计算
        gradient[i + 1](j) += temp_v * vi(j);
      }
    }
  }

  /* acceleration feasibility  加速度同理*/
  for (int i = 0; i < q.size() - 2; i++) {
    Eigen::Vector3d ai = q[i + 2] - 2 * q[i + 1] + q[i];

    for (int j = 0; j < 3; j++) {
      double ad = ai(j) * ai(j) * ts_inv4 - am2;
      if (ad > 0.0) {
        cost += pow(ad, 2);

        double temp_a = 4.0 * ad * ts_inv4;
        gradient[i + 0](j) += temp_a * ai(j);
        gradient[i + 1](j) += -2 * temp_a * ai(j);
        gradient[i + 2](j) += temp_a * ai(j);
      }
    }
  }
}
//计算端点值的代价
void BsplineOptimizer::calcEndpointCost(const vector<Eigen::Vector3d>& q, double& cost,
                                        vector<Eigen::Vector3d>& gradient) {
  cost = 0.0;
  Eigen::Vector3d zero(0, 0, 0);
  std::fill(gradient.begin(), gradient.end(), zero);

  // zero cost and gradient in hard constraints 通常，B-spline 曲线的末端点并不仅仅是最后一个控制点，而是由最后几个控制点的加权平均值计算出来的，通常采用 B-spline 样条的端点公式。
  Eigen::Vector3d q_3, q_2, q_1, dq;  //q_3, q_2, q_1 分别表示控制点向量 q 中倒数第三、第二和第一的点。端点约束通常是与这些控制点相关的。
  q_3 = q[q.size() - 3];
  q_2 = q[q.size() - 2];
  q_1 = q[q.size() - 1];
  //end_pt_ 是预设的目标端点，它表示优化时希望 B-spline 曲线到达的期望终点位置
  dq = 1 / 6.0 * (q_3 + 4 * q_2 + q_1) - end_pt_;
  cost += dq.squaredNorm();

  gradient[q.size() - 3] += 2 * dq * (1 / 6.0);   
  gradient[q.size() - 2] += 2 * dq * (4 / 6.0);
  gradient[q.size() - 1] += 2 * dq * (1 / 6.0);
}
//关键点的代价计算 与上面计算端点相同
void BsplineOptimizer::calcWaypointsCost(const vector<Eigen::Vector3d>& q, double& cost,
                                         vector<Eigen::Vector3d>& gradient) {
  cost = 0.0;
  Eigen::Vector3d zero(0, 0, 0);
  std::fill(gradient.begin(), gradient.end(), zero);

  Eigen::Vector3d q1, q2, q3, dq;

  // for (auto wp : waypoints_) {
  for (int i = 0; i < waypoints_.size(); ++i) {
    Eigen::Vector3d waypt = waypoints_[i];
    int             idx   = waypt_idx_[i];

    q1 = q[idx];
    q2 = q[idx + 1];
    q3 = q[idx + 2];

    dq = 1 / 6.0 * (q1 + 4 * q2 + q3) - waypt;
    cost += dq.squaredNorm();

    gradient[idx] += dq * (2.0 / 6.0);      // 2*dq*(1/6)
    gradient[idx + 1] += dq * (8.0 / 6.0);  // 2*dq*(4/6)
    gradient[idx + 2] += dq * (2.0 / 6.0);
  }
}

/* use the uniformly sampled points on a geomertic path to guide the
 * trajectory. For each control points to be optimized, it is assigned a
 * guiding point on the path and the distance between them is penalized */
//计算当前B样条的整体路径与想要贴近的路径的误差
void BsplineOptimizer::calcGuideCost(const vector<Eigen::Vector3d>& q, double& cost,
                                     vector<Eigen::Vector3d>& gradient) {
  cost = 0.0;
  Eigen::Vector3d zero(0, 0, 0);
  std::fill(gradient.begin(), gradient.end(), zero);

  int end_idx = q.size() - order_;

  for (int i = order_; i < end_idx; i++) {
    Eigen::Vector3d gpt = guide_pts_[i - order_];
    cost += (q[i] - gpt).squaredNorm();
    gradient[i] += 2 * (q[i] - gpt);
  }
}
//计算所有的误差结合
void BsplineOptimizer::combineCost(const std::vector<double>& x, std::vector<double>& grad,
                                   double& f_combine) {
  /* convert the NLopt format vector to control points. */

  // This solver can support 1D-3D B-spline optimization, but we use Vector3d to store each control point
  // For 1D case, the second and third elements are zero, and similar for the 2D case.
  //通过 control_points_（通常是 B-spline 曲线的控制点）初始化 g_q_。g_q_ 是一个存储控制点坐标的容器，大小为 order_。
  //这里将 B-spline 的前几个控制点从 control_points_ 中拷贝到 g_q_。 这个过程会将维度大于 1 的控制点（如3D）按 1D、2D 或 3D 填充，如果维度小于 3，则剩余维度设置为0。
  for (int i = 0; i < order_; i++) {
    for (int j = 0; j < dim_; ++j) {
      g_q_[i][j] = control_points_(i, j);
    }
    for (int j = dim_; j < 3; ++j) {
      g_q_[i][j] = 0.0;
    }
  }
  /*
    x 是由优化算法（如 NLopt）传递的变量，它包含了需要优化的控制点的坐标。
    该段代码从 x 中提取相应的控制点坐标，并将其存入 g_q_ 中的相应位置。g_q_ 将保存从 order_ 开始的控制点。
    同样，维度小于 3 的情况下，将剩余的维度填充为0。
  */
  for (int i = 0; i < variable_num_ / dim_; i++) {
    for (int j = 0; j < dim_; ++j) {
      g_q_[i + order_][j] = x[dim_ * i + j];
    }
    for (int j = dim_; j < 3; ++j) {
      g_q_[i + order_][j] = 0.0;
    }
  }
  //这段代码的作用是处理 B-spline 曲线的端点约束，特别是在没有启用端点约束时，将 B-spline 曲线的尾部控制点从 control_points_ 中复制到优化过程中使用的控制点数组 g_q_ 中。
  if (!(cost_function_ & ENDPOINT)) {
    for (int i = 0; i < order_; i++) {

      for (int j = 0; j < dim_; ++j) {
        g_q_[order_ + variable_num_ / dim_ + i][j] =
            control_points_(control_points_.rows() - order_ + i, j);
      }
      for (int j = dim_; j < 3; ++j) {
        g_q_[order_ + variable_num_ / dim_ + i][j] = 0.0;
      }
    }
  }

  f_combine = 0.0;
  grad.resize(variable_num_);
  fill(grad.begin(), grad.end(), 0.0);

  /*  evaluate costs and their gradient  */
  double f_smoothness, f_distance, f_feasibility, f_endpoint, f_guide, f_waypoints;
  f_smoothness = f_distance = f_feasibility = f_endpoint = f_guide = f_waypoints = 0.0;
  //计算各种误差以及与权重系数相乘 同时计算梯度 将各项代价值累加以及梯度累加 用于优化
  if (cost_function_ & SMOOTHNESS) {
    calcSmoothnessCost(g_q_, f_smoothness, g_smoothness_);
    f_combine += lambda1_ * f_smoothness;
    for (int i = 0; i < variable_num_ / dim_; i++)
      for (int j = 0; j < dim_; j++) grad[dim_ * i + j] += lambda1_ * g_smoothness_[i + order_](j);   //更新梯度
  }
  if (cost_function_ & DISTANCE) {
    calcDistanceCost(g_q_, f_distance, g_distance_);
    f_combine += lambda2_ * f_distance;
    for (int i = 0; i < variable_num_ / dim_; i++)
      for (int j = 0; j < dim_; j++) grad[dim_ * i + j] += lambda2_ * g_distance_[i + order_](j);
  }
  if (cost_function_ & FEASIBILITY) {
    calcFeasibilityCost(g_q_, f_feasibility, g_feasibility_);
    f_combine += lambda3_ * f_feasibility;
    for (int i = 0; i < variable_num_ / dim_; i++)
      for (int j = 0; j < dim_; j++) grad[dim_ * i + j] += lambda3_ * g_feasibility_[i + order_](j);
  }
  if (cost_function_ & ENDPOINT) {
    calcEndpointCost(g_q_, f_endpoint, g_endpoint_);
    f_combine += lambda4_ * f_endpoint;
    for (int i = 0; i < variable_num_ / dim_; i++)
      for (int j = 0; j < dim_; j++) grad[dim_ * i + j] += lambda4_ * g_endpoint_[i + order_](j);
  }
  if (cost_function_ & GUIDE) {
    calcGuideCost(g_q_, f_guide, g_guide_);
    f_combine += lambda5_ * f_guide;
    for (int i = 0; i < variable_num_ / dim_; i++)
      for (int j = 0; j < dim_; j++) grad[dim_ * i + j] += lambda5_ * g_guide_[i + order_](j);
  }
  if (cost_function_ & WAYPOINTS) {
    calcWaypointsCost(g_q_, f_waypoints, g_waypoints_);
    f_combine += lambda7_ * f_waypoints;
    for (int i = 0; i < variable_num_ / dim_; i++)
      for (int j = 0; j < dim_; j++) grad[dim_ * i + j] += lambda7_ * g_waypoints_[i + order_](j);
  }
  /*  print cost  */
  // if ((cost_function_ & WAYPOINTS) && iter_num_ % 10 == 0) {
  //   cout << iter_num_ << ", total: " << f_combine << ", acc: " << lambda8_ * f_view
  //        << ", waypt: " << lambda7_ * f_waypoints << endl;
  // }

  // if (optimization_phase_ == SECOND_PHASE) {
  //  << ", smooth: " << lambda1_ * f_smoothness
  //  << " , dist:" << lambda2_ * f_distance
  //  << ", fea: " << lambda3_ * f_feasibility << endl;
  // << ", end: " << lambda4_ * f_endpoint
  // << ", guide: " << lambda5_ * f_guide
  // }
}
//在调用 BsplineOptimizer 类的优化函数时，func_data 被设置为指向该优化器实例的指针。  个人理解 在使用优化库进行优化时 目标函数应该按照指定格式编写
//这样，通过 reinterpret_cast，我们可以将 func_data 转换回指向 BsplineOptimizer 对象的指针，从而访问该对象的成员函数和数据。
//我们就需要将 func_data 转换回原来的类型（BsplineOptimizer*），以便访问优化器的成员。通过 reinterpret_cast
double BsplineOptimizer::costFunction(const std::vector<double>& x, std::vector<double>& grad,
                                      void* func_data) {
  BsplineOptimizer* opt = reinterpret_cast<BsplineOptimizer*>(func_data); //使用 reinterpret_cast 将 func_data 转换为 BsplineOptimizer* 类型。
  double            cost;
  opt->combineCost(x, grad, cost);
  opt->iter_num_++;

  /* save the min cost result */
  if (cost < opt->min_cost_) {
    opt->min_cost_      = cost;
    opt->best_variable_ = x;
  }
  return cost;

  // /* evaluation */
  // ros::Time te1 = ros::Time::now();
  // double time_now = (te1 - opt->time_start_).toSec();
  // opt->vec_time_.push_back(time_now);
  // if (opt->vec_cost_.size() == 0)
  // {
  //   opt->vec_cost_.push_back(f_combine);
  // }
  // else if (opt->vec_cost_.back() > f_combine)
  // {
  //   opt->vec_cost_.push_back(f_combine);
  // }
  // else
  // {
  //   opt->vec_cost_.push_back(opt->vec_cost_.back());
  // }
}
//矩阵转换为向量
vector<Eigen::Vector3d> BsplineOptimizer::matrixToVectors(const Eigen::MatrixXd& ctrl_pts) {
  vector<Eigen::Vector3d> ctrl_q;
  for (int i = 0; i < ctrl_pts.rows(); ++i) {
    ctrl_q.push_back(ctrl_pts.row(i));
  }
  return ctrl_q;
}
//获取控制点
Eigen::MatrixXd BsplineOptimizer::getControlPoints() { return this->control_points_; }
//根据当前的 cost_function_ 值判断是否是二次优化问题。具体的代码逻辑如下：
bool BsplineOptimizer::isQuadratic() {
  if (cost_function_ == GUIDE_PHASE) {
    return true;
  } else if (cost_function_ == (SMOOTHNESS | WAYPOINTS)) {
    return true;
  }
  return false;
}

}  // namespace fast_planner