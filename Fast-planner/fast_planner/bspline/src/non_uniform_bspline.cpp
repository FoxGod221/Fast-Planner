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
/*
时间节点向量将整个b样条曲线进行了分段 共


*/


#include "bspline/non_uniform_bspline.h"
#include <ros/ros.h>

namespace fast_planner {

NonUniformBspline::NonUniformBspline(const Eigen::MatrixXd& points, const int& order,
                                     const double& interval) {
  setUniformBspline(points, order, interval);     //直接将样条初始化为均匀B样条。
}

NonUniformBspline::~NonUniformBspline() {}
//均匀B样条设置
void NonUniformBspline::setUniformBspline(const Eigen::MatrixXd& points, const int& order,
                                          const double& interval) {
                                            
  //三阶bspline p=3 假设n=3(控制点数-1)  m=3+3+1=7   n是控制点个数-1  并非一定是3 这里举个例子 是整段b样条曲线的控制点个数减1
  control_points_ = points;
  p_              = order;
  interval_       = interval;

  n_ = points.rows() - 1;     //n为控制点的个数-1  其中控制点向量矩阵的行数表示控制点的个数
  m_ = n_ + p_ + 1;
  /*初始化存储节点的向量
  前p_+1个节点被称为开头节点或者起始节点，他们决定了B样条的起始点，这样定义出来的起始节点为小于0的区间
  中间的p_ + 1 到 m_ - p_ 之间节点为中间节点，由累加而成。确保连续
  后m_ - p_ 个节点称为尾节点（trailing nodes）或结束节点，它们决定了 B 样条曲线的结束点，初始化时，这些节点的值与前一个节点相等，确保在这一段区间内的节点是连续的。这个区间的长度是 (m_ - p_) 
  节点向量用于参数化样条曲线的范围，确保每个控制点都在某个参数区间内产生影响。  前后部分重复节点(前后p+1个节点)用于起止控制。
  */
  u_ = Eigen::VectorXd::Zero(m_ + 1);   //时间节点为从t0到tm    一维向量
  for (int i = 0; i <= m_; ++i) {       //(-3,-2,-1,0,1,2,3,4) 假设只有4个控制点时

    if (i <= p_) {
      u_(i) = double(-p_ + i) * interval_;    //i<=3 -3 -2 -1 0
    } else if (i > p_ && i <= m_ - p_) {
      u_(i) = u_(i - 1) + interval_;          //i==4 1
    } else if (i > m_ - p_) {
      u_(i) = u_(i - 1) + interval_;          //i>4 2 3 4 
    }
  }
}

void NonUniformBspline::setKnot(const Eigen::VectorXd& knot) { this->u_ = knot; }   //设置时间节点

Eigen::VectorXd NonUniformBspline::getKnot() { return this->u_; }     //获取当前时间节点

void NonUniformBspline::getTimeSpan(double& um, double& um_p) {       //获取时间节点的横跨范围
  um   = u_(p_);  //第 p+1 个节点，对应样条的头部点。
  um_p = u_(m_ - p_); //倒数第p+1 个节点，对应样条的尾部点。
}

Eigen::MatrixXd NonUniformBspline::getControlPoint() { return control_points_; }    //返回控制点

pair<Eigen::VectorXd, Eigen::VectorXd> NonUniformBspline::getHeadTailPts() {      //
  //De Boor递归算法，计算传入参数u处的B样条曲线值，返回B-样条曲线在参数 u 处的B样条曲线函
  Eigen::VectorXd head = evaluateDeBoor(u_(p_));    //表示点的坐标，可以适应任意维度的 B 样条
  Eigen::VectorXd tail = evaluateDeBoor(u_(m_ - p_));
  return make_pair(head, tail);   //打包返回
}
/// @brief 这段代码是 NonUniformBspline 类中 evaluateDeBoor 方法的实现利用 De Boor 算法 在给定参数u下计算非均匀 B 样条曲线的值
/// @param u  样条曲线的参数值（介于节点向量范围内）。
/// @return 返回一个 Eigen::VectorXd，表示在参数u 处的曲线值（位置坐标）。
Eigen::VectorXd NonUniformBspline::evaluateDeBoor(const double& u) {

  double ub = min(max(u_(p_), u), u_(m_ - p_));   //将u限制在节点向量的合法范围内。 B 样条的定义只在上述范围内有效，超出范围可能导致非法计算。

  // determine which [ui,ui+1] lay in
  int k = p_;
  while (true) {                
    if (u_(k + 1) >= ub) break;  //找到u所在的区间 通过节点大小来判断 例如u(4)大于等于ub 也就是在在第四个区间内 也就是在u(k)之后 即第k+1的区间 用来确定时哪4个控制点影响他
    ++k;
  }

  /* deBoor's alg */
  vector<Eigen::VectorXd> d;
  //初始化控制点，用于 De Boor 算法的递推。
  for (int i = 0; i <= p_; ++i) {     
    d.push_back(control_points_.row(k - p_ + i));   //控制点索引范围为 [k−p,k]，共 p+1 个控制点 每相邻的四个控制点控制一个时间区间内的轨迹。 通过k找到与当前节点相关的控制点k-p_到k
    // cout << d[i].transpose() << endl;
  }
  //递推次数为p 次（样条的阶数） 根据的是b样条的数学递推关系
  for (int r = 1; r <= p_; ++r) {
    for (int i = p_; i >= r; --i) {
      double alpha = (ub - u_[i + k - p_]) / (u_[i + 1 + k - r] - u_[i + k - p_]);
      // cout << "alpha: " << alpha << endl;
      d[i] = (1 - alpha) * d[i - 1] + alpha * d[i];
    }
  }

  return d[p_];   //返回计算得到的 B-样条曲线在参数 u 处的点
}
//它调用了 evaluateDeBoor 方法，并对输入参数进行了偏移操作 原理是b样条的特性
Eigen::VectorXd NonUniformBspline::evaluateDeBoorT(const double& t) {
  return evaluateDeBoor(t + u_(p_));
}

//它计算并返回 B 样条曲线的导数控制点。 计算方法为b样条的导数计算公式 数学推导 利用原控制点进行计算
Eigen::MatrixXd NonUniformBspline::getDerivativeControlPoints() {
  // The derivative of a b-spline is also a b-spline, its order become p_-1
  // control point Qi = p_*(Pi+1-Pi)/(ui+p_+1-ui+1)
  Eigen::MatrixXd ctp = Eigen::MatrixXd::Zero(control_points_.rows() - 1, control_points_.cols());  //导数控制点个数比原曲线控制点个数少1 是根据它的导数计算公式得来的
  for (int i = 0; i < ctp.rows(); ++i) {
    ctp.row(i) =
        p_ * (control_points_.row(i + 1) - control_points_.row(i)) / (u_(i + p_ + 1) - u_(i + 1));
  }
  return ctp;
}

//计算并返回当前 B 样条曲线的导数（即 B 样条的导数曲线）。
NonUniformBspline NonUniformBspline::getDerivative() {
  Eigen::MatrixXd   ctp = getDerivativeControlPoints();   //首先得到导数的控制点
  NonUniformBspline derivative(ctp, p_ - 1, interval_);   //构造新的b样条 阶数比原曲线低1阶 时间间隔保持相同

  /* cut the first and last knot */
  Eigen::VectorXd knot(u_.rows() - 2);    //因为阶数减1 控制点也减1 所以节点向量比原来少2
  knot = u_.segment(1, u_.rows() - 2);    //导数曲线的节点向量需要删除原始节点向量的第一个和最后一个节点（因为导数曲线的节点范围是去掉原始曲线的首尾节点）
  derivative.setKnot(knot);   //设置节点向量

  return derivative;
}
//获取时间间隔  
double NonUniformBspline::getInterval() { return interval_; }
//设置动力学限制
void NonUniformBspline::setPhysicalLimits(const double& vel, const double& acc) {
  limit_vel_   = vel;
  limit_acc_   = acc;
  limit_ratio_ = 1.1;
}
//检查可行性
bool NonUniformBspline::checkFeasibility(bool show) {
  bool fea = true;
  // SETY << "[Bspline]: total points size: " << control_points_.rows() << endl;

  Eigen::MatrixXd P         = control_points_;
  int             dimension = control_points_.cols();

  /* check vel feasibility and insert points */
  double max_vel = -1.0;  //计算两个相邻点啊之间的速度
  for (int i = 0; i < P.rows() - 1; ++i) {
    Eigen::VectorXd vel = p_ * (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1));

    if (fabs(vel(0)) > limit_vel_ + 1e-4 || fabs(vel(1)) > limit_vel_ + 1e-4 ||
        fabs(vel(2)) > limit_vel_ + 1e-4) {

      if (show) cout << "[Check]: Infeasible vel " << i << " :" << vel.transpose() << endl;   //如果超过了限制 选择是否打印是哪个点超出了
      fea = false;

      for (int j = 0; j < dimension; ++j) {
        max_vel = max(max_vel, fabs(vel(j)));   //记录了超过限制的最大速度
      }
    }
  }

  /* acc feasibility  计算方法以及判断方法与上面相同*/
  double max_acc = -1.0;
  for (int i = 0; i < P.rows() - 2; ++i) {

    Eigen::VectorXd acc = p_ * (p_ - 1) *
        ((P.row(i + 2) - P.row(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
         (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
        (u_(i + p_ + 1) - u_(i + 2));

    if (fabs(acc(0)) > limit_acc_ + 1e-4 || fabs(acc(1)) > limit_acc_ + 1e-4 ||
        fabs(acc(2)) > limit_acc_ + 1e-4) {

      if (show) cout << "[Check]: Infeasible acc " << i << " :" << acc.transpose() << endl;
      fea = false;

      for (int j = 0; j < dimension; ++j) {
        max_acc = max(max_acc, fabs(acc(j)));
      }
    }
  }

  double ratio = max(max_vel / limit_vel_, sqrt(fabs(max_acc) / limit_acc_));   //计算速度和加速度的比率  在这个函数内没用任何作用

  return fea;
}

//检查比率 这个函数与上面的函数功能有些重叠
double NonUniformBspline::checkRatio() {
  Eigen::MatrixXd P         = control_points_;
  int             dimension = control_points_.cols();

  // find max vel
  double max_vel = -1.0;
  for (int i = 0; i < P.rows() - 1; ++i) {
    Eigen::VectorXd vel = p_ * (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1));
    for (int j = 0; j < dimension; ++j) {
      max_vel = max(max_vel, fabs(vel(j)));         //记录最大速度
    }
  }
  // find max acc
  double max_acc = -1.0;
  for (int i = 0; i < P.rows() - 2; ++i) {
    Eigen::VectorXd acc = p_ * (p_ - 1) *
        ((P.row(i + 2) - P.row(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
         (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
        (u_(i + p_ + 1) - u_(i + 2));
    for (int j = 0; j < dimension; ++j) {
      max_acc = max(max_acc, fabs(acc(j)));       //记录最大加速度
    }
  }
  double ratio = max(max_vel / limit_vel_, sqrt(fabs(max_acc) / limit_acc_));       //记录最大速度与限制速度的比率
  ROS_ERROR_COND(ratio > 2.0, "max vel: %lf, max acc: %lf.", max_vel, max_acc);     //超速报错

  return ratio;
}
//时间重分配 将B样条时间重分配为非均匀b样条
bool NonUniformBspline::reallocateTime(bool show) {
  // SETY << "[Bspline]: total points size: " << control_points_.rows() << endl;
  // cout << "origin knots:\n" << u_.transpose() << endl;
  bool fea = true;

  Eigen::MatrixXd P         = control_points_;
  int             dimension = control_points_.cols();

  double max_vel, max_acc;

  /* check vel feasibility and insert points */
  for (int i = 0; i < P.rows() - 1; ++i) {
    Eigen::VectorXd vel = p_ * (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1));

    if (fabs(vel(0)) > limit_vel_ + 1e-4 || fabs(vel(1)) > limit_vel_ + 1e-4 ||
        fabs(vel(2)) > limit_vel_ + 1e-4) {

      fea = false;                                                  //可行性检查
      if (show) cout << "[Realloc]: Infeasible vel " << i << " :" << vel.transpose() << endl;

      max_vel = -1.0;
      for (int j = 0; j < dimension; ++j) {
        max_vel = max(max_vel, fabs(vel(j)));
      }

      double ratio = max_vel / limit_vel_ + 1e-4;     //计算最大速度与限制速度的比率 只需要对最大速度进行限制 其他速度一定会满足条件
      if (ratio > limit_ratio_) ratio = limit_ratio_;       //设置最大速度比率限制

      double time_ori = u_(i + p_ + 1) - u_(i + 1);     //计算原始时间间隔  
      double time_new = ratio * time_ori;           //计算重新分配的新的时间间隔
      double delta_t  = time_new - time_ori;      //计算重分配的插值
      double t_inc    = delta_t / double(p_);     //计算均匀地将调整的时间差 delta_t 分配到相邻的 p_ 个时间节点之间，从而使得时间节点间隔更加平滑。

      for (int j = i + 2; j <= i + p_ + 1; ++j) {     //只修改有超出限制的部分的时间分配 因此从i+2开始进行修改
        u_(j) += double(j - i - 1) * t_inc;
        if (j <= 5 && j >= 1) {
          // cout << "vel j: " << j << endl;
        }
      }

      for (int j = i + p_ + 2; j < u_.rows(); ++j) {    //将后面不受影响的节点向后整体平移
        u_(j) += delta_t;
      }
    }
  }

  /* acc feasibility */
  for (int i = 0; i < P.rows() - 2; ++i) {

    Eigen::VectorXd acc = p_ * (p_ - 1) *
        ((P.row(i + 2) - P.row(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
         (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
        (u_(i + p_ + 1) - u_(i + 2));

    if (fabs(acc(0)) > limit_acc_ + 1e-4 || fabs(acc(1)) > limit_acc_ + 1e-4 ||
        fabs(acc(2)) > limit_acc_ + 1e-4) {

      fea = false;
      if (show) cout << "[Realloc]: Infeasible acc " << i << " :" << acc.transpose() << endl;

      max_acc = -1.0;
      for (int j = 0; j < dimension; ++j) {
        max_acc = max(max_acc, fabs(acc(j)));
      }

      double ratio = sqrt(max_acc / limit_acc_) + 1e-4;   //记录加速度最超出的比率 与速度最大超出比率的计算略有不同 原因是因为速度的影响区间比加速度的影响区间要小
      if (ratio > limit_ratio_) ratio = limit_ratio_;
      // cout << "ratio: " << ratio << endl;

      double time_ori = u_(i + p_ + 1) - u_(i + 2);
      double time_new = ratio * time_ori;
      double delta_t  = time_new - time_ori;
      double t_inc    = delta_t / double(p_ - 1); //与速度分配也有不同

      if (i == 1 || i == 2) {     //加速度影响两个区间 因此做特殊处理
        // cout << "acc i: " << i << endl;
        for (int j = 2; j <= 5; ++j) {
          u_(j) += double(j - 1) * t_inc;
        }

        for (int j = 6; j < u_.rows(); ++j) {
          u_(j) += 4.0 * t_inc;   //后续节点需要按一个固定的比例来调整，以保持平滑过渡。 这里4.0也是相当于对i==1 i==2时的特殊处理
        }

      } else {

        for (int j = i + 3; j <= i + p_ + 1; ++j) {
          u_(j) += double(j - i - 2) * t_inc;
          if (j <= 5 && j >= 1) {
            // cout << "acc j: " << j << endl;
          }
        }

        for (int j = i + p_ + 2; j < u_.rows(); ++j) {      //这里与速度调整策略保持一致
          u_(j) += delta_t;
        }
      }
    }
  }

  return fea;
}

//通过将时间节点重新分布，使得某个区间内的时间间隔变得更长，从而“拉伸”或“延长”该区间的时间。 这个函数调整了后面所有的节点 程序中也没有用到
void NonUniformBspline::lengthenTime(const double& ratio) {
  int num1 = 5;
  int num2 = getKnot().rows() - 1 - 5;

  double delta_t = (ratio - 1.0) * (u_(num2) - u_(num1));
  double t_inc   = delta_t / double(num2 - num1);
  for (int i = num1 + 1; i <= num2; ++i) u_(i) += double(i - num1) * t_inc;
  for (int i = num2 + 1; i < u_.rows(); ++i) u_(i) += delta_t;
}

void NonUniformBspline::recomputeInit() {}

//作用是根据给定的一组点（point_set）和起始/结束的导数（start_end_derivative），以及时间步长（ts），来计算对应的 B样条曲线控制点（ctrl_pts）。
//在这段代码的开头部分，它主要做了一些输入检查。
//这个函数很重要，起着承上启下的作用，在第一部分动态搜索得到了离散的路经点，此函数作用为求取控制点矩阵，该控制点矩阵生成的B样条曲线拟合路径点的路径。

/*
函数参数：ts：轨迹执行时间、point_set：原轨迹点集合、start_end_derivative：起点与终点的高阶约束
输出：更新ctrl_pts控制点矩阵的值
函数功能：将给定点集和起始/终止导数转换为 B-样条曲线的控制点矩阵，通过对原始轨迹的拟合得到B样条轨迹的控制点
*/

void NonUniformBspline::parameterizeToBspline(const double& ts, const vector<Eigen::Vector3d>& point_set,
                                              const vector<Eigen::Vector3d>& start_end_derivative,
                                              Eigen::MatrixXd&               ctrl_pts) {
  if (ts <= 0) {
    cout << "[B-spline]:time step error." << endl;
    return;
  }

  if (point_set.size() < 2) {
    cout << "[B-spline]:point set have only " << point_set.size() << " points." << endl;
    return;
  }
  //起始和结束的导数）是否包含 4 个元素。在 B样条曲线的构造过程中，通常需要在起始和结束点上提供一阶和二阶导数的值（例如，起始点的导数和结束点的导数）。
  //如果导数的数量不等于 4（即 2 个点的导数，每个点包含 2 个导数值），则输出错误信息。
  if (start_end_derivative.size() != 4) {
    cout << "[B-spline]:derivatives error." << endl;
  }
  //
//记录原曲线路经点个数
  int K = point_set.size();

  // write A 矩阵 A 的构造，矩阵 A 将用于计算 B样条的控制点。
  //三个 3D 向量，它们分别表示用于计算矩阵 A 中的权重系数：
  //该矩阵用来构建线性方程组， B-样条曲线参数化过程中求解控制点
  //一阶 B-样条基函数的系数向量、一阶 B-样条基函数的导数系数向量、二阶 B-样条基函数的系数向量
  //取该数值的原因是B样条曲线矩阵表示论文中提出的，以满足 B-样条曲线的一些良好性质。 为了将时间区间都映射到0到4上 便于计算 可以参考B样条原论文

  /*
  
  第一行表示了位置约束关系，第二行是速度、第三行是加速度。当传入的离散轨迹点为K时，有K-1段轨迹，由于是三次B样条曲线，
  所以算上头部3个节点，此外为了保证曲线的连续性后面添2个节点（因为要保证速度与加速度约束），共有K+5个时间节点，可以得到所求控制点个数应为K+5-次数3为K+2个。
  然后通过上式构建AX=B方程即可求解出控制点向量X

  相当于给定了前端A*的路径轨迹点 以及需要的动力学约束数据 以及确定的时间间隔来构建方程 AX=B  求解得到控制点的位置
  
  */
 //这是M矩阵  对应B样条计算p(s(t))=s(t)MQm 中的M矩阵 其中A矩阵就是s(t)M  Qm也就是X  B就是前端传入的点集 所以想要得到控制点的方法 就是解AX=B的方程
  Eigen::Vector3d prow(3), vrow(3), arow(3);
  prow << 1, 4, 1;
  vrow << -1, 0, 1;
  arow << 1, -2, 1;
  //A 是一个 K + 4 行、K + 2 列的矩阵，初始化为零。矩阵的大小比点集的大小大，因为在矩阵中还包含了速度和加速度的信息。  因为要包含起始点和目标点 
  //以及起始点与目标点的速度加速度约束  这里列数表示k+2个点 行数表示k+4个约束
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(K + 4, K + 2);

  for (int i = 0; i < K; ++i) A.block(i, i, 1, 3) = (1 / 6.0) * prow.transpose();

  A.block(K, 0, 1, 3)         = (1 / 2.0 / ts) * vrow.transpose();
  A.block(K + 1, K - 1, 1, 3) = (1 / 2.0 / ts) * vrow.transpose();

  A.block(K + 2, 0, 1, 3)     = (1 / ts / ts) * arow.transpose();
  A.block(K + 3, K - 1, 1, 3) = (1 / ts / ts) * arow.transpose();
  // cout << "A:\n" << A << endl;

  // A.block(0, 0, K, K + 2) = (1 / 6.0) * A.block(0, 0, K, K + 2);
  // A.block(K, 0, 2, K + 2) = (1 / 2.0 / ts) * A.block(K, 0, 2, K + 2);
  // A.row(K + 4) = (1 / ts / ts) * A.row(K + 4);
  // A.row(K + 5) = (1 / ts / ts) * A.row(K + 5);

  // write b
  Eigen::VectorXd bx(K + 4), by(K + 4), bz(K + 4);
  for (int i = 0; i < K; ++i) {
    bx(i) = point_set[i](0);
    by(i) = point_set[i](1);
    bz(i) = point_set[i](2);
  }

  for (int i = 0; i < 4; ++i) {
    bx(K + i) = start_end_derivative[i](0);
    by(K + i) = start_end_derivative[i](1);
    bz(K + i) = start_end_derivative[i](2);
  }

  // solve Ax = b
  Eigen::VectorXd px = A.colPivHouseholderQr().solve(bx);
  Eigen::VectorXd py = A.colPivHouseholderQr().solve(by);
  Eigen::VectorXd pz = A.colPivHouseholderQr().solve(bz);

  // convert to control pts
  ctrl_pts.resize(K + 2, 3);
  ctrl_pts.col(0) = px;
  ctrl_pts.col(1) = py;
  ctrl_pts.col(2) = pz;

  // cout << "[B-spline]: parameterization ok." << endl;
}
//获取B样条时间跨度
double NonUniformBspline::getTimeSum() {
  double tm, tmp;
  getTimeSpan(tm, tmp);
  return tmp - tm;
}
//获取B样条曲线的总长度 方法是通过对曲线上的每个小时间段（根据给定的分辨率 res）进行离散采样，并计算相邻采样点之间的距离，然后将这些距离累加得到总长度。
double NonUniformBspline::getLength(const double& res) {
  double          length = 0.0;
  double          dur    = getTimeSum();
  Eigen::VectorXd p_l    = evaluateDeBoorT(0.0), p_n;   //p_l存储当前的曲线点，初始化为曲线在时间点 0.0 处的值。  p_n用于存储下一个时间点的曲线点。
  for (double t = res; t <= dur + 1e-4; t += res) {
    p_n = evaluateDeBoorT(t);
    length += (p_n - p_l).norm();       //累加求和
    p_l = p_n;
  }
  return length;
}
//jerk的计算方式同样是利用b样条的倒数递推性质来进行计算
double NonUniformBspline::getJerk() {
  NonUniformBspline jerk_traj = getDerivative().getDerivative().getDerivative();  //方法三次，依次计算了 B样条曲线的导数（速度）、二阶导数（加速度）和三阶导数（jerk）。

  Eigen::VectorXd times     = jerk_traj.getKnot();    //times 获取的是 B样条曲线的时间节点（即节点向量）。
  Eigen::MatrixXd ctrl_pts  = jerk_traj.getControlPoint();    //获取的是jerk B样条曲线的控制点矩阵
  int             dimension = ctrl_pts.cols();

  double jerk = 0.0;
  for (int i = 0; i < ctrl_pts.rows(); ++i) {
    for (int j = 0; j < dimension; ++j) {
      jerk += (times(i + 1) - times(i)) * ctrl_pts(i, j) * ctrl_pts(i, j);    //该表达式的含义是时间差（即相邻时间节点之间的间隔）乘以控制点在该维度的值的平方。
    }
  }

  return jerk;    //返回积分的jerk 
}
//计算非均匀 B 样条曲线的平均速度（mean_v）和最大速度（max_v）
void NonUniformBspline::getMeanAndMaxVel(double& mean_v, double& max_v) {
  NonUniformBspline vel = getDerivative();
  double            tm, tmp;
  vel.getTimeSpan(tm, tmp);

  double max_vel = -1.0, mean_vel = 0.0;
  int    num = 0;
  for (double t = tm; t <= tmp; t += 0.01) {
    Eigen::VectorXd vxd = vel.evaluateDeBoor(t);
    double          vn  = vxd.norm();     //当前速度值

    mean_vel += vn;       //累加为了求平均
    ++num;
    if (vn > max_vel) {               //选择排序找最大值
      max_vel = vn;
    }
  }

  mean_vel = mean_vel / double(num);    //  求平均值
  mean_v   = mean_vel;
  max_v    = max_vel;
}
//计算加速度的平均值以及最大值 与上面类似
void NonUniformBspline::getMeanAndMaxAcc(double& mean_a, double& max_a) {
  NonUniformBspline acc = getDerivative().getDerivative();
  double            tm, tmp;
  acc.getTimeSpan(tm, tmp);

  double max_acc = -1.0, mean_acc = 0.0;
  int    num = 0;
  for (double t = tm; t <= tmp; t += 0.01) {
    Eigen::VectorXd axd = acc.evaluateDeBoor(t);
    double          an  = axd.norm();

    mean_acc += an;
    ++num;
    if (an > max_acc) {
      max_acc = an;
    }
  }

  mean_acc = mean_acc / double(num);
  mean_a   = mean_acc;
  max_a    = max_acc;
}
}  // namespace fast_planner
