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

#include <path_searching/kinodynamic_astar.h>
#include <sstream>
#include <plan_env/sdf_map.h>

using namespace std;
using namespace Eigen;

namespace fast_planner
{
//析构函数 用来释放内存
KinodynamicAstar::~KinodynamicAstar()
{
  for (int i = 0; i < allocate_num_; i++)
  {
    delete path_node_pool_[i];
  }
}
//混合Astar 启发函数采用的是庞特里亚金计算的最优路径计算出来的代价 同时还有最优时间也作为代价的一项 f=g+h
int KinodynamicAstar::search(Eigen::Vector3d start_pt, Eigen::Vector3d start_v, Eigen::Vector3d start_a,
                             Eigen::Vector3d end_pt, Eigen::Vector3d end_v, bool init, bool dynamic, double time_start)
{
  start_vel_ = start_v;
  start_acc_ = start_a;
  //从 path_node_pool_（一个预分配的节点池）中取出第一个节点，作为起点节点 cur_node。 内存池 相当于已经提前申请好了一块内存节点
  //这种方式避免了频繁的动态内存分配，提高效率。
  PathNodePtr cur_node = path_node_pool_[0];
  cur_node->parent = NULL;        //无父节点
  cur_node->state.head(3) = start_pt; //给该节点的初始状态赋值
  cur_node->state.tail(3) = start_v;
  cur_node->index = posToIndex(start_pt); //获取当前第一个节点在地图中的索引
  cur_node->g_score = 0.0;  //g值记为0
  //末状态也就是目标节点初始化
  Eigen::VectorXd end_state(6);
  Eigen::Vector3i end_index;
  double time_to_goal;    //记录从开始到目标节点花费的时间

  end_state.head(3) = end_pt;
  end_state.tail(3) = end_v;
  end_index = posToIndex(end_pt);
  cur_node->f_score = lambda_heu_ * estimateHeuristic(cur_node->state, end_state, time_to_goal);  //由于初始节点并没有相应的g代价 因此启发值就是初始节点的总代价 
  cur_node->node_state = IN_OPEN_SET;   //同时将当前节点加入到oplt中 等待扩展 状态量也进行修改 openlist是由一个优先队列进行维护的
  open_set_.push(cur_node);
  use_node_num_ += 1;           //对使用的节点数进行++

  //dynamic 是一个布尔参数，指示当前的搜索是否在动态环境中进行。
  if (dynamic)
  {
    time_origin_ = time_start;  //将 time_start（搜索的起始时间）赋值给成员变量 time_origin_，用于记录搜索的时间基准。
    cur_node->time = time_start;    //设置起点节点的时间属性 time 为搜索的起始时间。
    cur_node->time_idx = timeToIndex(time_start);   //调用函数 timeToIndex 将时间 time_start 离散化为时间索引 time_idx。
    expanded_nodes_.insert(cur_node->index, cur_node->time_idx, cur_node);    //将当前节点插入到已扩展节点内。
    // cout << "time start: " << time_start << endl;
  }
  else    //当 dynamic 为 false 时，表示当前搜索是在静态环境中进行，时间维度不需要考虑。 
    expanded_nodes_.insert(cur_node->index, cur_node);

  PathNodePtr neighbor = NULL;    //存储当前扩展节点的邻居节点的指针（初始化为 NULL）。
  PathNodePtr terminate_node = NULL;    //记录满足终止条件的节点（如果找到终止节点，将存储在这里）。
  bool init_search = init;      //标志变量，用于记录搜索是否为初始化状态，便于在不同情况下处理搜索逻辑。
  const int tolerance = ceil(1 / resolution_);  //取整得到容忍范围的网格单元数，用于判断当前节点是否接近目标点。

  //A*核心部分
  while (!open_set_.empty())
  {
    //后续将检查是否满足终止条件，并扩展其邻居节点。
    cur_node = open_set_.top(); //取出优先队列中的队头 也就是当前oplt中所有节点的最小f值的节点  当前正在处理的节点，

    // Terminate?
    //检查当前节点是否已经超出了搜索的范围（horizon_）。
    bool reach_horizon = (cur_node->state.head(3) - start_pt).norm() >= horizon_;    //计算的是当前节点与起点之间的欧氏距离是否超过了设定值
    //对每个坐标分量，比较当前节点的离散索引与目标点索引之间的差值，判断是否小于等于 tolerance。 如果在容忍范围内，认为已接近目标点。
    bool near_end = abs(cur_node->index(0) - end_index(0)) <= tolerance &&      
                    abs(cur_node->index(1) - end_index(1)) <= tolerance &&
                    abs(cur_node->index(2) - end_index(2)) <= tolerance;

    if (reach_horizon || near_end)    //如果满足在设定距离范围内 或者距离很近  则认为已经到达了目标点
    {
      terminate_node = cur_node;      // 满足条件的终止节点记为当前节点
      retrievePath(terminate_node);   //根据当前节点回溯其父节点链表，提取整个路径并存储在类的成员变量中path_nodes_。
      if (near_end)   //如果当前节点接近目标点（near_end），需要进一步检查是否存在one shot  也就是shot轨迹
      {
        // Check whether shot traj exist    进行shot尝试 one shot成功后 就无需再进行后续的扩展搜索
        estimateHeuristic(cur_node->state, end_state, time_to_goal);      //计算当前节点到目标点的启发式代价和估算时间 time_to_goal。
        computeShotTraj(cur_node->state, end_state, time_to_goal);
        if (init_search)      //如果一开始就找到了轨迹 认为不正确
          ROS_ERROR("Shot in first search loop!");
      }
    }
    if (reach_horizon)        //在设定距离外
    {
      if (is_shot_succ_)    //成功shot的话 表示当前搜索已经找到一条有效路径，返回状态 REACH_END。
      {
        std::cout << "reach end" << std::endl;
        return REACH_END;
      }
      else    //否则表示
      {
        std::cout << "reach horizon" << std::endl;        //如果未生成成功，表示搜索到达范围边界但未找到有效路径，返回状态 REACH_HORIZON。
        return REACH_HORIZON;
      }
    }

    if (near_end)   //如果到达目标点附近
    {
      if (is_shot_succ_)    //成功shot的话
      {
        std::cout << "reach end" << std::endl;
        return REACH_END;
      }
      else if (cur_node->parent != NULL)      //如果当前节点接近目标点，但未成功生成shot，可以返回状态 NEAR_END，表示搜索到接近目标点的区域。
      {
        std::cout << "near end" << std::endl;
        return NEAR_END;
      }
      else
      {
        std::cout << "no path" << std::endl;      //当前节点接近目标点，但由于没有父节点，搜索失败，返回状态 NO_PATH，表示无效路径。
        return NO_PATH;
      }
    }
    open_set_.pop();      //弹出当前节点 表示该节点不是要找的节点
    cur_node->node_state = IN_CLOSE_SET;    //改变当前节点的状态 表示已经访问过了
    iter_num_ += 1;       //迭代次数++
    /*
    这里的主要想法是 判断当前节点是否已经靠近了目标节点 以及shot到 如果没有 再进行扩展
    */
    //为扩展节点做参数准备
    /*
    res：空间分辨率，定义为 0.5 单位长度。
    用于将连续空间坐标映射到离散空间网格。
    time_res：时间分辨率 用于时间离散化的步长。
    time_res_init：初始化的时间分辨率，定义为 通常用于在初始化阶段生成更精细的时间步长
    */
    double res = 1 / 2.0, time_res = 1 / 1.0, time_res_init = 1 / 20.0;
    Eigen::Matrix<double, 6, 1> cur_state = cur_node->state;    
    Eigen::Matrix<double, 6, 1> pro_state;      //用于存储当前节点扩展出的候选邻居节点的状态（位置和速度）。
    vector<PathNodePtr> tmp_expand_nodes;   //临时存储当前节点扩展出的所有候选邻居节点。
    Eigen::Vector3d um;   //当前节点扩展时使用的控制输入（加速度）  用于扩展新节点的输入数据 在控制空间下进行节点的扩展 通过动力学扩展公式 进行子节点的扩展
    double pro_t;       //扩展时的预测时间 相当于扩展到某个节点的时间
    vector<Eigen::Vector3d> inputs;     //存储所有可能的控制输入（加速度）。  
    vector<double> durations;   //一个向量容器，存储扩展过程中所有可能的时间步长。
    if (init_search)    //如果是初始搜索
    {
      inputs.push_back(start_acc_);   //使用初始状态的控制量  压入容器
      for (double tau = time_res_init * init_max_tau_; tau <= init_max_tau_ + 1e-3;     //作用是通过给定的时间步长 在最大时间步长范围每 最小步长 
           tau += time_res_init * init_max_tau_)
        durations.push_back(tau);     //向容器中加入一次步长记录
      init_search = false;    //初始搜索准备工作结束
    }
    //非首次搜索 将控制输入也按照最大最小范围内 通过步长进行分割 最后压入容器内 
    else
    {
      for (double ax = -max_acc_; ax <= max_acc_ + 1e-3; ax += max_acc_ * res)
        for (double ay = -max_acc_; ay <= max_acc_ + 1e-3; ay += max_acc_ * res)
          for (double az = -max_acc_; az <= max_acc_ + 1e-3; az += max_acc_ * res)
          {
            um << ax, ay, az;
            inputs.push_back(um);
          }
      for (double tau = time_res * max_tau_; tau <= max_tau_; tau += time_res * max_tau_) //如果不是初始搜索 时间最小分辨率会相比于初始搜索的粗糙一点
        durations.push_back(tau);
    }

    // cout << "cur state:" << cur_state.head(3).transpose() << endl;  //组合不同的时间与加速度组合，得到不一样的下一状态集合
    for (int i = 0; i < inputs.size(); ++i)   //每种控制输入都进行扩展尝试
      for (int j = 0; j < durations.size(); ++j)  //尝试不同步长的扩展
      {
        um = inputs[i];
        double tau = durations[j];
        stateTransit(cur_state, pro_state, um, tau);    //进行扩展的计算 计算得到的扩展节点存储pro中 到也就是论文中带有指数函数以及积分的公式项
        pro_t = cur_node->time + tau;   //计算出的邻居节点的时间点为当前节点的时间加上步长

        Eigen::Vector3d pro_pos = pro_state.head(3);    //扩展出的节点的位置状态 

        // Check if in close set
        Eigen::Vector3i pro_id = posToIndex(pro_pos);   //记录空间索引
        int pro_t_id = timeToIndex(pro_t);      //时间索引
        //每轮 也就是确定的控制输入以及确定的步长 只会得到一个扩展点
        PathNodePtr pro_node = dynamic ? expanded_nodes_.find(pro_id, pro_t_id) : expanded_nodes_.find(pro_id);   //区分是否为动态环境 用来记录这一轮扩展得到的最佳节点
        //下面的几个if主要哦用来对当前节点是否可用进行检查 主要包括是否在openlist中 以及是否满足动力学（速度）要求 是否在同一网格内 是否经过障碍物
        if (pro_node != NULL && pro_node->node_state == IN_CLOSE_SET)   //如果该节点已经被扩展过 或者搜索过 直接跳过
        {
          if (init_search)        //首次搜索时提示
            std::cout << "close" << std::endl;
          continue;
        }

        // Check maximal velocity
        Eigen::Vector3d pro_v = pro_state.tail(3);
        if (fabs(pro_v(0)) > max_vel_ || fabs(pro_v(1)) > max_vel_ || fabs(pro_v(2)) > max_vel_)
        {
          if (init_search)
            std::cout << "vel" << std::endl;
          continue;
        }

        // Check not in the same voxel    对同一栅格内的节点进行剪枝
        Eigen::Vector3i diff = pro_id - cur_node->index;      //计算扩展的节点与当前节点的索引差
        int diff_time = pro_t_id - cur_node->time_idx;      //计算时间索引差值  表示扩展节点和当前节点处于同一个时间体素。
        //这里的剪枝策略与后面的会有所不同 由于这里的pro还只是根据控制量计算出来的节点 是根据当前节点扩展出来的节点 因此g值一定小于当前节点 所以直接进行剪枝
        if (diff.norm() == 0 && ((!dynamic) || diff_time == 0))       //如果处于同一体素网格 或者时间体素内 进行剪枝 
        {
          if (init_search)
            std::cout << "same" << std::endl;
          continue;
        }

        // Check safety 
        Eigen::Vector3d pos;    //用于存储预测路径上某点的位置（3D 坐标）。
        Eigen::Matrix<double, 6, 1> xt;   //用于存储预测路径上的状态（6维向量），包括位置和速度。
        bool is_occ = false;    //是否碰撞了
        for (int k = 1; k <= check_num_; ++k)   //成员变量，表示路径上需要检查的采样点数目。 对路径进行采样
        {
          double dt = tau * double(k) / double(check_num_);   //采样步长为之前计算扩展步长的分段
          stateTransit(cur_state, xt, um, dt);            //计算采样点的状态
          pos = xt.head(3);     //记录位置
          if (edt_environment_->sdf_map_->getInflateOccupancy(pos) == 1 )   //调用环境地图中的方法检查位置 pos 是否处在障碍物区域。 返回1为障碍物
          {
            is_occ = true;
            break;
          }
        }
        if (is_occ)   //如果碰撞了 也跳过这个点
        {
          if (init_search)
            std::cout << "safe" << std::endl;
          continue;
        }
        //到达目标的时间  以及临时变量
        /*
        um.squaredNorm()：控制输入um的平方范数，表示控制的力度。在规划中，将控制输入作为代价的一部分，鼓励选择较小的控制输入以节省能量。
        w_time_：时间代价的权重，用于平衡路径时间和空间控制输入的代价。
        这里计算的当前节点扩展到pro节点 主要考虑的代价是控制输入大小 以及耗费时间长短 f值计算正常
        */
        double time_to_goal, tmp_g_score, tmp_f_score;
        tmp_g_score = (um.squaredNorm() + w_time_) * tau + cur_node->g_score; //扩展节点g值等于当前节点g+当前节点到扩展结点的代价
        tmp_f_score = tmp_g_score + lambda_heu_ * estimateHeuristic(pro_state, end_state, time_to_goal);

        // Compare nodes expanded from the same parent      扩展节点的剪枝 留下最适合的扩展节点
        /*
        这里对应论文里的扩展结点剪枝策略 通过判断这个节点是否已经在临时扩展节点的序列内 
        如果在 判断是否f值小于序列内已经存在的 如果小于 更新序列内的节点状态 如果不小于 舍弃这个节点
        */
        bool prune = false;   
        for (int j = 0; j < tmp_expand_nodes.size(); ++j)   //遍历 tmp_expand_nodes 中的每个节点 expand_node，检查新扩展节点是否重复。 也就是pro是否已经在临时扩展的节点内了
        {
          PathNodePtr expand_node = tmp_expand_nodes[j];
          if ((pro_id - expand_node->index).norm() == 0 && ((!dynamic) || pro_t_id == expand_node->time_idx))   //通过空间索引以及时间索引判断是否在序列内
          {
            prune = true;       //表示该节点进行过了剪枝
            if (tmp_f_score < expand_node->f_score)
            {
              expand_node->f_score = tmp_f_score;
              expand_node->g_score = tmp_g_score;
              expand_node->state = pro_state;
              expand_node->input = um;
              expand_node->duration = tau;
              if (dynamic)
                expand_node->time = cur_node->time + tau;
            }
            break;
          }
        }

        // This node end up in a voxel different from others
        if (!prune)   //如果这个节点没有进行过剪枝
        {
          if (pro_node == NULL)   //如果扩展节点 pro_node 尚未被创建（pro_node == NULL），则初始化新节点。
          {
            pro_node = path_node_pool_[use_node_num_];
            pro_node->index = pro_id;
            pro_node->state = pro_state;
            pro_node->f_score = tmp_f_score;
            pro_node->g_score = tmp_g_score;
            pro_node->input = um;
            pro_node->duration = tau;
            pro_node->parent = cur_node;
            pro_node->node_state = IN_OPEN_SET;
            if (dynamic)
            {
              pro_node->time = cur_node->time + tau;
              pro_node->time_idx = timeToIndex(pro_node->time);
            }
            open_set_.push(pro_node);   //加入到openlist序列 

            if (dynamic)    //将这个节点记录为已经被扩展过的节点 加入序列
              expanded_nodes_.insert(pro_id, pro_node->time, pro_node);
            else
              expanded_nodes_.insert(pro_id, pro_node);

            tmp_expand_nodes.push_back(pro_node);   //临时存储当前扩展的节点，用于后续节点比较或代价更新。

            use_node_num_ += 1;   //每成功扩展一个节点，已使用的节点计数器增加 1。
            if (use_node_num_ == allocate_num_)   //节点池的最大分配数量。如果节点池用尽，打印错误信息并终止搜索。
            {
              cout << "run out of memory." << endl;
              return NO_PATH;
            }
          }
          else if (pro_node->node_state == IN_OPEN_SET)   //同样的 如果在openlist中已经有过这个节点  采用判断g值的方式来选择是否更新节点
          {
            if (tmp_g_score < pro_node->g_score)
            {
              // pro_node->index = pro_id;
              pro_node->state = pro_state;
              pro_node->f_score = tmp_f_score;
              pro_node->g_score = tmp_g_score;
              pro_node->input = um;
              pro_node->duration = tau;
              pro_node->parent = cur_node;
              if (dynamic)
                pro_node->time = cur_node->time + tau;
            }
          }
          else
          {
            cout << "error type in searching: " << pro_node->node_state << endl;
          }
        }
      }
    // init_search = false;
  }
  //到最后也没有节点可以加入openlist 或者openlist已经空了 返回没有路径
  cout << "open set empty, no path!" << endl;
  cout << "use node num: " << use_node_num_ << endl;
  cout << "iter num: " << iter_num_ << endl;
  return NO_PATH;
}
//参数设置 用于从 ROS 的参数服务器中读取配置参数并初始化类的成员变量。
void KinodynamicAstar::setParam(ros::NodeHandle& nh)
{
  nh.param("search/max_tau", max_tau_, -1.0);   //从参数服务器中读取指定参数。 参数路径为 "search/max_tau"。 如果参数未找到，则将默认值 -1.0 赋给变量 max_tau_。
  nh.param("search/init_max_tau", init_max_tau_, -1.0);
  nh.param("search/max_vel", max_vel_, -1.0);
  nh.param("search/max_acc", max_acc_, -1.0);
  nh.param("search/w_time", w_time_, -1.0);
  nh.param("search/horizon", horizon_, -1.0);
  nh.param("search/resolution_astar", resolution_, -1.0);
  nh.param("search/time_resolution", time_resolution_, -1.0);
  nh.param("search/lambda_heu", lambda_heu_, -1.0);
  nh.param("search/allocate_num", allocate_num_, -1);
  nh.param("search/check_num", check_num_, -1);
  nh.param("search/optimistic", optimistic_, true);
  tie_breaker_ = 1.0 + 1.0 / 10000;   //打破 A* 搜索中节点代价相同的平局情况。 默认值为增加一个很小的偏移量。

  double vel_margin;  
  nh.param("search/vel_margin", vel_margin, 0.0);
  max_vel_ += vel_margin;     //将 max_vel_ 增加一个裕量，提升路径规划的灵活性。
}
//从输入节点向前回溯 以此来获得整条路径
void KinodynamicAstar::retrievePath(PathNodePtr end_node)
{
  PathNodePtr cur_node = end_node;
  path_nodes_.push_back(cur_node);

  while (cur_node->parent != NULL)
  {
    cur_node = cur_node->parent;
    path_nodes_.push_back(cur_node);
  }

  reverse(path_nodes_.begin(), path_nodes_.end());
}
// 启发式函数设计 论文中采用了庞特里亚金极小值作为启发值
//这里并非是像论文中求代价的导数等于得出最小值 而是直接暴力比较所有解大小 得出最小的
double KinodynamicAstar::estimateHeuristic(Eigen::VectorXd x1, Eigen::VectorXd x2, double& optimal_time)
{
  const Vector3d dp = x2.head(3) - x1.head(3);    //位置差
  const Vector3d v0 = x1.segment(3, 3);   //segment 从索引3开始 提取长度为3的子集
  const Vector3d v1 = x2.segment(3, 3);

  //此处为代价函数的系数计算  此处前面系数可能是权重
  double c1 = -36 * dp.dot(dp); //dot为向量内积函数
  double c2 = 24 * (v0 + v1).dot(dp);   //速度和位置差的内积，用于体现速度在目标方向上的投影。
  double c3 = -4 * (v0.dot(v0) + v0.dot(v1) + v1.dot(v1));
  double c4 = 0;      //四次项系数，这里为 0，表示没有四次项。
  double c5 = w_time_;    //时间权重 w_time，影响代价函数中的时间项。

  std::vector<double> ts = quartic(c5, c4, c3, c2, c1);   //求一个关于时间的轨迹的解方程  这些时间值对应不同的轨迹可能性。
  //这里基于最大速度取一个速度值 然后计算一个不考虑动力学约束的最短时间 在后续计算中 只要小于这个的时间 都不考虑
  double v_max = max_vel_ * 0.5;
  double t_bar = (x1.head(3) - x2.head(3)).lpNorm<Infinity>() / v_max;
  ts.push_back(t_bar);    //将 t_bar 添加到时间候选值集合中，确保时间不会小于该值。

  double cost = 100000000;  //代价初始值
  double t_d = t_bar;   //初始化为 t_bar，表示默认最优时间为不受约束的时间下界。

  for (auto t : ts)  //直接遍历所有的t来自ts集合中
  {
    if (t < t_bar)
      continue;
    double c = -c1 / (3 * t * t * t) - c2 / (2 * t * t) - c3 / t + w_time_ * t;    //完整的代价函数计算
    if (c < cost)   //过滤掉 t<t_bar 的情况，因为这些时间违反了动力学约束。
    {
      cost = c;
      t_d = t;
    }
  }

  optimal_time = t_d;   // 记录最优时间

  return 1.0 * (1 + tie_breaker_) * cost;   //返回启发式代价，包含一个系数 tie_breaker_，用于打破搜索中代价相同节点的排序。
}
/////////////////////////
/*
    shot 思想   用于计算两个位置之间的平滑轨迹  该轨迹同样需要做动力学限制检测 碰撞检测
*/
bool KinodynamicAstar::computeShotTraj(Eigen::VectorXd state1, Eigen::VectorXd state2, double time_to_goal)
{
  /* ---------- get coefficient 获取传入的参数 ---------- */
  const Vector3d p0 = state1.head(3);
  const Vector3d dp = state2.head(3) - p0;
  const Vector3d v0 = state1.segment(3, 3);
  const Vector3d v1 = state2.segment(3, 3);
  const Vector3d dv = v1 - v0;
  double t_d = time_to_goal;
  MatrixXd coef(3, 4);
  end_vel_ = v1;
  // //计算轨迹系数 使用的是庞特里亚金的最优轨迹 三次多项式
  Vector3d a = 1.0 / 6.0 * (-12.0 / (t_d * t_d * t_d) * (dp - v0 * t_d) + 6 / (t_d * t_d) * dv);   
  Vector3d b = 0.5 * (6.0 / (t_d * t_d) * (dp - v0 * t_d) - 2 / t_d * dv);
  Vector3d c = v0;
  Vector3d d = p0;

  // 1/6 * alpha * t^3 + 1/2 * beta * t^2 + v0
  // a*t^3 + b*t^2 + v0*t + p0    每一列存储 a,b,c,d 的值，对应三次多项式的每个维度。
  coef.col(3) = a, coef.col(2) = b, coef.col(1) = c, coef.col(0) = d;

  Vector3d coord, vel, acc;   //coord三维向量，用于存储当前轨迹点的位置坐标
  VectorXd poly1d, t, polyv, polya;   //多项式系数的存储向量 t为t变量的0 1 2 3次形式向量 
  Vector3i index; //三维索引，用于离散化坐标，将连续位置转换为栅格地图中的离散索引。
  /*
  Tm是一个矩阵，用于计算多项式的一阶和二阶导数（速度和加速度）。 v=Tm*x  a=Tm*(Tm*x)
  */
  Eigen::MatrixXd Tm(4, 4);     
  Tm << 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0;

  /* ---------- forward checking of trajectory 
  
  检查是否满足物理学约束 也即速度 加速度 障碍物
  
  ---------- */
  double t_delta = t_d / 10;      //同样采取分段采样进行检测
  for (double time = t_delta; time <= t_d; time += t_delta)
  {
    //将t的各个阶次的计算结果 存入t 向量中 便于后续计算
    t = VectorXd::Zero(4);
    for (int j = 0; j < 4; j++)
      t(j) = pow(time, j);
    //针对 x,y,z 三个维度分别计算轨迹点的状态。
    for (int dim = 0; dim < 3; dim++)
    {
      poly1d = coef.row(dim);     //取出单独这个维度的系数
      coord(dim) = poly1d.dot(t);   //计算得到这一维度的多项式
      vel(dim) = (Tm * poly1d).dot(t);    //求导得到速度  
      acc(dim) = (Tm * Tm * poly1d).dot(t);   //求导得到加速度

      if (fabs(vel(dim)) > max_vel_ || fabs(acc(dim)) > max_acc_)   //做动力学检测
      {
        // cout << "vel:" << vel(dim) << ", acc:" << acc(dim) << endl;
        // return false;
      }
    }
    //做边界检测 不能超出地图的边界
    if (coord(0) < origin_(0) || coord(0) >= map_size_3d_(0) || coord(1) < origin_(1) || coord(1) >= map_size_3d_(1) ||
        coord(2) < origin_(2) || coord(2) >= map_size_3d_(2))
    {
      return false;
    }

    // if (edt_environment_->evaluateCoarseEDT(coord, -1.0) <= margin_) {
    //   return false;
    // }
    //碰撞检测
    if (edt_environment_->sdf_map_->getInflateOccupancy(coord) == 1)
    {
      return false;
    }
  }
  //如果没有问题 获取系数 表示成功生成了平滑的轨迹
  coef_shot_ = coef;
  t_shot_ = t_d;
  is_shot_succ_ = true;
  return true;
}
//这段代码实现了求解三次多项式实数根的函数 利用求根公式暴力求解
vector<double> KinodynamicAstar::cubic(double a, double b, double c, double d)
{
  vector<double> dts;

  double a2 = b / a;
  double a1 = c / a;
  double a0 = d / a;

  double Q = (3 * a1 - a2 * a2) / 9;
  double R = (9 * a1 * a2 - 27 * a0 - 2 * a2 * a2 * a2) / 54;
  double D = Q * Q * Q + R * R;
  if (D > 0)
  {
    double S = std::cbrt(R + sqrt(D));
    double T = std::cbrt(R - sqrt(D));
    dts.push_back(-a2 / 3 + (S + T));
    return dts;
  }
  else if (D == 0)
  {
    double S = std::cbrt(R);
    dts.push_back(-a2 / 3 + S + S);
    dts.push_back(-a2 / 3 - S);
    return dts;
  }
  else
  {
    double theta = acos(R / sqrt(-Q * Q * Q));
    dts.push_back(2 * sqrt(-Q) * cos(theta / 3) - a2 / 3);
    dts.push_back(2 * sqrt(-Q) * cos((theta + 2 * M_PI) / 3) - a2 / 3);
    dts.push_back(2 * sqrt(-Q) * cos((theta + 4 * M_PI) / 3) - a2 / 3);
    return dts;
  }
}
//用于求解四次多项式方程的实数根  通过调用cubic函数解决部分系数的关系，再进一步分解为两个二次方程
vector<double> KinodynamicAstar::quartic(double a, double b, double c, double d, double e)
{
  vector<double> dts;

  double a3 = b / a;
  double a2 = c / a;
  double a1 = d / a;
  double a0 = e / a;
  //四次多项式通过一个三次方程的根y1来进一步简化。
  vector<double> ys = cubic(1, -a2, a1 * a3 - 4 * a0, 4 * a2 * a0 - a1 * a1 - a3 * a3 * a0);
  double y1 = ys.front();
  double r = a3 * a3 / 4 - a2 + y1;
  if (r < 0)
    return dts;

  double R = sqrt(r);
  double D, E;
  if (R != 0)
  {
    D = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 + 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
    E = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 - 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
  }
  else
  {
    D = sqrt(0.75 * a3 * a3 - 2 * a2 + 2 * sqrt(y1 * y1 - 4 * a0));
    E = sqrt(0.75 * a3 * a3 - 2 * a2 - 2 * sqrt(y1 * y1 - 4 * a0));
  }

  if (!std::isnan(D))
  {
    dts.push_back(-a3 / 4 + R / 2 + D / 2);
    dts.push_back(-a3 / 4 + R / 2 - D / 2);
  }
  if (!std::isnan(E))
  {
    dts.push_back(-a3 / 4 - R / 2 + E / 2);
    dts.push_back(-a3 / 4 - R / 2 - E / 2);
  }

  return dts;
}
//初始化  
void KinodynamicAstar::init()
{
  /* ---------- map params ---------- */
  this->inv_resolution_ = 1.0 / resolution_;    //地图分辨率的倒数，用于快速进行空间坐标与离散网格索引之间的转换
  inv_time_resolution_ = 1.0 / time_resolution_;  //时间分辨率的导数
  edt_environment_->sdf_map_->getRegion(origin_, map_size_3d_); //获取地图原点以及尺寸等信息

  cout << "origin_: " << origin_.transpose() << endl;
  cout << "map size: " << map_size_3d_.transpose() << endl;

  /* ---------- pre-allocated node 
  提前分配的节点池，存储路径搜索中用到的节点，避免动态内存分配的开销。 主要用在扩展与搜索节点时---------- */
  path_node_pool_.resize(allocate_num_);
  for (int i = 0; i < allocate_num_; i++)
  {
    path_node_pool_[i] = new PathNode;
  }
  //状态转移矩阵，6x6 单位矩阵。 表示系统的状态转移关系。 在进行控制点扩展时结合步长 控制输入一起使用
  phi_ = Eigen::MatrixXd::Identity(6, 6);
  use_node_num_ = 0;
  iter_num_ = 0;
}
//在路径规划开始前，调用此函数将具体的环境数据与规划器关联起来，为后续搜索提供地图、障碍物等信息。
void KinodynamicAstar::setEnvironment(const EDTEnvironment::Ptr& env)
{
  this->edt_environment_ = env;
}
//重置路径规划器的所有内部状态，为新一轮搜索任务做准备。
void KinodynamicAstar::reset()
{
  expanded_nodes_.clear();
  path_nodes_.clear();
  //替换为一个空队列 因为 std::priority_queue 没有提供直接的清空函数。创建一个空队列并与当前队列交换，可以高效清空队列。
  std::priority_queue<PathNodePtr, std::vector<PathNodePtr>, NodeComparator> empty_queue;
  open_set_.swap(empty_queue);
  //遍历所有已使用的节点，将其状态重置为初始状态
  for (int i = 0; i < use_node_num_; i++)
  {
    PathNodePtr node = path_node_pool_[i];
    node->parent = NULL;
    node->node_state = NOT_EXPAND;
  }
  //清空计数变量以及标志位
  use_node_num_ = 0;
  iter_num_ = 0;
  is_shot_succ_ = false;
  has_path_ = false;
}
//用于生成从起点到目标点的完整轨迹，包括常规搜索轨迹和射击轨迹（shot trajectory）

std::vector<Eigen::Vector3d> KinodynamicAstar::getKinoTraj(double delta_t)  //delta_t决定轨迹点的采样间隔。
{
  vector<Vector3d> state_list;    //返回值 包含整个轨迹的所有离散位置点。

  /* ---------- get traj of searching ---------- */
  PathNodePtr node = path_nodes_.back();      //从 path_nodes_ 的最后一个节点（目标点）开始，逐级回溯父节点，直到起点。
  Matrix<double, 6, 1> x0, xt;

  while (node->parent != NULL)
  {
    Vector3d ut = node->input;    //获取控制输入
    double duration = node->duration; //当前节点到父节点的传播时间。
    x0 = node->parent->state;     //父节点状态

    //通过时间步长 delta_t 将节点间的轨迹离散化。
    for (double t = duration; t >= -1e-5; t -= delta_t) 
    {
      stateTransit(x0, xt, ut, t);
      state_list.push_back(xt.head(3));   //每计算一次 存储一次计算的中间采样点状态
    }
    node = node->parent;
  }
  reverse(state_list.begin(), state_list.end());    //计算完成后 翻转序列
  /* ---------- get traj of one shot ---------- */
  //如果是生成的shot轨迹 直接根据多项式的系数来进行间隔采样计算 假如搜索到一半 shot成功了 那先将上面的已经搜索过得到的节点压入list 
  //然后将后半段的shot再压入list shot的不需要反转序列
  if (is_shot_succ_)
  {
    Vector3d coord;
    VectorXd poly1d, time(4);

    for (double t = delta_t; t <= t_shot_; t += delta_t)
    {
      for (int j = 0; j < 4; j++)
        time(j) = pow(t, j);

      for (int dim = 0; dim < 3; dim++)
      {
        poly1d = coef_shot_.row(dim);
        coord(dim) = poly1d.dot(time);
      }
      state_list.push_back(coord);
    }
  }

  return state_list;
}
/*用于提取整个路径的时间采样点，同时计算路径的起点和终点的速度、加速度。以下是详细解析：
////ts:时间采样间隔，用于离散路径。 
point_set:传引用的容器，用于存储路径的采样点（位置向量）。
start_end_derivatives:传引用的容器，用于存储起点和终点的速度、加速度。

*/
void KinodynamicAstar::getSamples(double& ts, vector<Eigen::Vector3d>& point_set,
                                  vector<Eigen::Vector3d>& start_end_derivatives)
{
  /* ---------- path duration ---------- */
  double T_sum = 0.0;
  //如果射击轨迹（shot trajectory）成功，则加上射击轨迹的时长 t_shot_。
  if (is_shot_succ_)
    T_sum += t_shot_;
  PathNodePtr node = path_nodes_.back();
  //遍历回溯搜索路径节点，累加每段路径的持续时间 duration，得到路径总时长。
  while (node->parent != NULL)
  {
    T_sum += node->duration;
    node = node->parent;
  }
  // cout << "duration:" << T_sum << endl;

  // Calculate boundary vel and acc
  Eigen::Vector3d end_vel, end_acc;
  double t;
  //shot轨迹的速度与加速度直接使用多项式计算
  if (is_shot_succ_)
  {
    t = t_shot_;
    end_vel = end_vel_;
    for (int dim = 0; dim < 3; ++dim)
    {
      Vector4d coe = coef_shot_.row(dim);
      end_acc(dim) = 2 * coe(2) + 6 * coe(3) * t_shot_;
    }
  }
  //从轨迹的终点提取速度信息与加速度（输入）
  else
  {
    t = path_nodes_.back()->duration;
    end_vel = node->state.tail(3);
    end_acc = path_nodes_.back()->input;
  }

  // Get point samples
  int seg_num = floor(T_sum / ts);  //采样间隔时间 使用 floor(T_sum / ts) 确定段数，
  seg_num = max(8, seg_num);  //保证最少有8段
  ts = T_sum / double(seg_num); //更新采样时间步长 例如有少于8段的情况
  bool sample_shot_traj = is_shot_succ_;
  node = path_nodes_.back();

  for (double ti = T_sum; ti > -1e-5; ti -= ts)
  {
    if (sample_shot_traj)   //如果有shot轨迹，先进行采样。
    {
      // samples on shot traj
      Vector3d coord;
      Vector4d poly1d, time;

      for (int j = 0; j < 4; j++)
        time(j) = pow(t, j);

      for (int dim = 0; dim < 3; dim++)
      {
        poly1d = coef_shot_.row(dim);
        coord(dim) = poly1d.dot(time);
      }

      point_set.push_back(coord);   //尾插 t从大往小计算的
      t -= ts;

      /* end of segment */
      if (t < -1e-5)    //shot采样结束
      {
        sample_shot_traj = false;
        if (node->parent != NULL)
          t += node->duration;
      }
    }
    else
    {
      // samples on searched traj 同样是进行间隔采样
      Eigen::Matrix<double, 6, 1> x0 = node->parent->state;
      Eigen::Matrix<double, 6, 1> xt;
      Vector3d ut = node->input;

      stateTransit(x0, xt, ut, t);

      point_set.push_back(xt.head(3));
      t -= ts;

      // cout << "t: " << t << ", t acc: " << T_accumulate << endl;
      if (t < -1e-5 && node->parent->parent != NULL)
      {
        node = node->parent;
        t += node->duration;
      }
    }
  }
  reverse(point_set.begin(), point_set.end());//反转

  // calculate start acc
  //如果只有shot轨迹（没有搜索轨迹），通过射击轨迹的二次项系数计算起点加速度。  
  Eigen::Vector3d start_acc;
  if (path_nodes_.back()->parent == NULL)   //表示的是shot轨迹 因为shot轨迹是直接多项式计算出来的 没有父节点
  {
    // no searched traj, calculate by shot traj
    start_acc = 2 * coef_shot_.col(2);
  }
  else
  {
    // input of searched traj
    start_acc = node->input;
  }
  //保存变量
  start_end_derivatives.push_back(start_vel_);
  start_end_derivatives.push_back(end_vel);
  start_end_derivatives.push_back(start_acc);
  start_end_derivatives.push_back(end_acc);
}
//
std::vector<PathNodePtr> KinodynamicAstar::getVisitedNodes()
{
  vector<PathNodePtr> visited;
  visited.assign(path_node_pool_.begin(), path_node_pool_.begin() + use_node_num_ - 1);
  return visited;
}
//将位置坐标转换为索引
Eigen::Vector3i KinodynamicAstar::posToIndex(Eigen::Vector3d pt)
{
  Vector3i idx = ((pt - origin_) * inv_resolution_).array().floor().cast<int>();

  // idx << floor((pt(0) - origin_(0)) * inv_resolution_), floor((pt(1) -
  // origin_(1)) * inv_resolution_),
  //     floor((pt(2) - origin_(2)) * inv_resolution_);

  return idx;
}
//时间转换为索引
int KinodynamicAstar::timeToIndex(double time)
{
  int idx = floor((time - time_origin_) * inv_time_resolution_);
  return idx;
}
/// @brief 计算指定步长以及指定控制输入的状态转换
/// @param state0 当前状态
/// @param state1 要计算的状态
/// @param um 控制输入
/// @param tau 步长
void KinodynamicAstar::stateTransit(Eigen::Matrix<double, 6, 1>& state0, Eigen::Matrix<double, 6, 1>& state1,
                                    Eigen::Vector3d um, double tau)
{
  for (int i = 0; i < 3; ++i)
    phi_(i, i + 3) = tau;

  Eigen::Matrix<double, 6, 1> integral;
  integral.head(3) = 0.5 * pow(tau, 2) * um;
  integral.tail(3) = tau * um;

  state1 = phi_ * state0 + integral;
}

}  // namespace fast_planner
