#ifndef _KINODYNAMIC_ASTAR_H
#define _KINODYNAMIC_ASTAR_H

// #include <path_searching/matrix_hash.h>
#include <ros/console.h>
#include <ros/ros.h>
#include <Eigen/Eigen>
#include <boost/functional/hash.hpp>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <unordered_map>
#include <utility>
#include "plan_env/edt_environment.h"

namespace fast_planner {
// #define REACH_HORIZON 1
// #define REACH_END 2
// #define NO_PATH 3
#define IN_CLOSE_SET 'a'
#define IN_OPEN_SET 'b'
#define NOT_EXPAND 'c'
#define inf 1 >> 30

class PathNode {
 public:
  /* -------------------- */
  Eigen::Vector3i index;
  Eigen::Matrix<double, 6, 1> state;      //fp中的状态量包括x y z vx vy vz
  double g_score, f_score;                //定义g值以及 f值
  Eigen::Vector3d input;                  //fp 中控制输入为加速度
  double duration;            
  double time;  // dyn
  int time_idx;
  PathNode* parent;     //父节点
  char node_state;      //当前节点状态 用来判断是再close_list中还是再open_list中

// 构造函数 每个节点构造时父节点先清空 状态表示为未扩展

  PathNode() {
    parent = NULL;
    node_state = NOT_EXPAND;
  }
  ~PathNode(){};
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
typedef PathNode* PathNodePtr; //新名字

//定义节点大小比较运算 定义的是函数对象  当创建了一个NodeComparator对象A 直接调用A（a,b） 会直接返回两个f的大小结果 相当于重载了()
class NodeComparator {
 public:
  bool operator()(PathNodePtr node1, PathNodePtr node2) {
    return node1->f_score > node2->f_score;
  }
};

//定义模板T
/*
matrix_hash 是一个自定义的哈希函数， 获取的是整个矩阵对应的哈希值
用于生成矩阵类型对象的哈希值。它通过遍历矩阵中的每个元素并利用元素的哈希值、一些常数和位运算来计算最终的哈希值。
这使得矩阵可以作为容器（例如 std::unordered_map）的键。
std::size_t operator()(T const& matrix) const {             //operator() 是重载了的函数调用运算符，使得 matrix_hash 可以像普通函数一样被调用。
  //它接受一个常量引用类型的 matrix 参数（即矩阵类型 T 的对象），并返回一个 std::size_t 类型的哈希值。
*/
template <typename T>
struct matrix_hash : std::unary_function<T, size_t> {
  std::size_t operator()(T const& matrix) const {             //operator() 是重载了的函数调用运算符，使得 matrix_hash 可以像普通函数一样被调用。
  //它接受一个常量引用类型的 matrix 参数（即矩阵类型 T 的对象），并返回一个 std::size_t 类型的哈希值。
    size_t seed = 0;
    for (size_t i = 0; i < matrix.size(); ++i) {
      auto elem = *(matrix.data() + i);
      seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) +
              (seed >> 2);
    }
    return seed;
  }
};

/***
 * 
 * 这段代码定义了一个名为 NodeHashTable 的类，它实现了一个专用的哈希表，支持对三维（Eigen::Vector3i）和四维（Eigen::Vector4i）索引的管理，
 * 存储与这些索引相关联的指针类型值（PathNodePtr）。以下是对代码的详细解析：
 * 
 * std::unordered_map 是 C++ 标准模板库（STL）中的一个关联容器，
 * 提供 键值对（key-value pairs） 的存储和快速查找功能。它是基于哈希表（Hash Table）实现的，具有 常数平均时间复杂度 的插入、查找和删除操作。
 * 其中键的类型是Eigen::Vector3i 值的类型是PathNodePtr 而键值对之间的哈希映射关系 未后面的matrix_hash函数
 */
class NodeHashTable {
 private:
  /* data */
  std::unordered_map<Eigen::Vector3i, PathNodePtr, matrix_hash<Eigen::Vector3i>>
      data_3d_;
  std::unordered_map<Eigen::Vector4i, PathNodePtr, matrix_hash<Eigen::Vector4i>>
      data_4d_;

 public:
  NodeHashTable(/* args */) {}
  ~NodeHashTable() {}                   //由于成员变量都是由 标准库STL构建出来的  因此不需要手动管理内存  

  //插入操作 也是直接调用stl的一些库函数
  void insert(Eigen::Vector3i idx, PathNodePtr node) { 
    data_3d_.insert(std::make_pair(idx, node));
  }
  //将 idx 与 time_idx 组合为一个四维索引 Eigen::Vector4i，再插入到 data_4d_ 中。  只是给三维索引多加了一个时间戳的索引
  void insert(Eigen::Vector3i idx, int time_idx, PathNodePtr node) {
    data_4d_.insert(std::make_pair(
        Eigen::Vector4i(idx(0), idx(1), idx(2), time_idx), node));
  }
  //查找函数 未找到返回null 找到了返回Node指针
  PathNodePtr find(Eigen::Vector3i idx) {
    auto iter = data_3d_.find(idx);
    return iter == data_3d_.end() ? NULL : iter->second;
  }
  PathNodePtr find(Eigen::Vector3i idx, int time_idx) {
    auto iter =
        data_4d_.find(Eigen::Vector4i(idx(0), idx(1), idx(2), time_idx));
    return iter == data_4d_.end() ? NULL : iter->second;
  }
  //清空数据
  void clear() {
    data_3d_.clear();
    data_4d_.clear();
  }
};
//A*所用到的数据结构
class KinodynamicAstar {
 private:
  /* ---------- main data structure ---------- */
  vector<PathNodePtr> path_node_pool_;    //路径节点内存池
  int use_node_num_, iter_num_;           //记录已经被使用作为路径节点的个数 以及算法的迭代次数
  NodeHashTable expanded_nodes_;          //使用一个哈希表存储已经扩展过的节点
  std::priority_queue<PathNodePtr, std::vector<PathNodePtr>, NodeComparator>      //使用优先队列作为open_list 参数分别是存储的元素数据类型  数据的序列存储格式为vector 以及存储规则按何种大小
      open_set_;
  std::vector<PathNodePtr> path_nodes_;     //存储规划路径的节点序列。 这个应该是最后需要的路径

  /* ---------- record data ---------- */
  Eigen::Vector3d start_vel_, end_vel_, start_acc_;       //初末状态速度 以及加速度
  Eigen::Matrix<double, 6, 6> phi_;  // state transit matrix    //状态变量矩阵 状态转移矩阵，表示从一个状态（位置、速度等）到下一个状态的动力学关系
  // shared_ptr<SDFMap> sdf_map;
  EDTEnvironment::Ptr edt_environment_;     //环境模型（基于 EDTEnvironment），提供对距离场或地图的访问，用于检测障碍物或估算代价。
  bool is_shot_succ_ = false;                 // 是否直接检测到目标成功   标志变量，指示是否成功计算了到目标的平滑连接轨迹（"shot trajectory"）。 投机one shot尝试
  Eigen::MatrixXd coef_shot_;       //记录平滑连接轨迹的系数和所需时间。
  double t_shot_;
  bool has_path_ = false;       //标志变量，指示是否成功找到一条有效路径。  非平滑 指图搜索的结果

  /* ---------- parameter ---------- */
  /* search */
  double max_tau_, init_max_tau_;     //控制状态传播的时间步长。
  double max_vel_, max_acc_;          //最大动力学参数
  double w_time_, horizon_, lambda_heu_;      //时间权重 horizon_搜索的时间范围 lambda_heu_启发式函数的权重，用来控制贪心程度
  int allocate_num_, check_num_;
  double tie_breaker_;      //A* 算法中打破优先队列中代价相同的节点顺序 obvp问题 
  bool optimistic_;         //启发式估计是否使用乐观估计。  

  /* map */  
  double resolution_, inv_resolution_, time_resolution_, inv_time_resolution_;    //空间分辨率及其倒数，用于将位置坐标映射为离散索引。  时间分辨率及其倒数。
  Eigen::Vector3d origin_, map_size_3d_;      //地图的原点和三维大小。
  double time_origin_;        //时间的起点，用于时间索引计算。

  /* helper */
  Eigen::Vector3i posToIndex(Eigen::Vector3d pt);     //将连续的三维空间坐标（Eigen::Vector3d）转换为离散网格索引
  int timeToIndex(double time);       //将时间点转换为离散的时间索引。
  void retrievePath(PathNodePtr end_node);      //从终点节点回溯，提取出整条路径并存储到 path_nodes_ 中

  /* shot trajectory */
  vector<double> cubic(double a, double b, double c, double d);     //计算三次多项式的系数，用于生成平滑轨迹。
  vector<double> quartic(double a, double b, double c, double d, double e);   //四次
  bool computeShotTraj(Eigen::VectorXd state1, Eigen::VectorXd state2,          //计算从当前状态到目标状态的平滑连接轨迹。
                       double time_to_goal);
  double estimateHeuristic(Eigen::VectorXd x1, Eigen::VectorXd x2,      //估计从一个状态到目标状态的启发式代价。
                           double& optimal_time);

  /* state propagation 
  实现状态的传播计算，给定当前状态 state0 和控制输入 um，计算经过时间步长 tau 后的新状态 state1。
  */
  void stateTransit(Eigen::Matrix<double, 6, 1>& state0,
                    Eigen::Matrix<double, 6, 1>& state1, Eigen::Vector3d um,
                    double tau);

 public:
  KinodynamicAstar(){};
  ~KinodynamicAstar();

  enum { REACH_HORIZON = 1, REACH_END = 2, NO_PATH = 3, NEAR_END = 4 };

  /* main API */
  void setParam(ros::NodeHandle& nh);   //从 ROS 参数服务器中获取算法参数并初始化相关变量。
  void init();
  void reset();
  //搜索函数，输入起点和终点的状态（位置、速度、加速度等），并搜索符合约束的路径。        返回值是一个整数，表示搜索状态：
  int search(Eigen::Vector3d start_pt, Eigen::Vector3d start_vel,
             Eigen::Vector3d start_acc, Eigen::Vector3d end_pt,
             Eigen::Vector3d end_vel, bool init, bool dynamic = false,
             double time_start = -1.0);

  void setEnvironment(const EDTEnvironment::Ptr& env);      //设置环境地图（EDTEnvironment），提供障碍物信息和距离场。

  std::vector<Eigen::Vector3d> getKinoTraj(double delta_t);   //返回规划的轨迹，每隔 delta_t 采样一个点。

  void getSamples(double& ts, vector<Eigen::Vector3d>& point_set,         //返回规划的轨迹，每隔 delta_t 采样一个点。 
                  vector<Eigen::Vector3d>& start_end_derivatives);          //getSamples：返回轨迹的采样点及其导数信息。

  std::vector<PathNodePtr> getVisitedNodes();     //getVisitedNodes：返回搜索过程中访问过的节点。

  typedef shared_ptr<KinodynamicAstar> Ptr;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace fast_planner

#endif