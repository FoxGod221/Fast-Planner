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

/**
 * fast_planner_node.cpp文件是我们通过launch文件启动Fast Planner程序的接口，也是ROS中Fast Planner的主节点，
 * 该文件在完成节点初始化后，会根据参数planner的值选择进行kino_replan或者topo_replan的初始化。
 * 主要用于初始化和选择不同的规划算法，具体来说，是通过有限状态机（FSM）进行重新规划。代码的流程如下： 
 * 
 */

#include <ros/ros.h>
#include <visualization_msgs/Marker.h>

#include <plan_manage/kino_replan_fsm.h>
#include <plan_manage/topo_replan_fsm.h>

#include <plan_manage/backward.hpp>
namespace backward {
backward::SignalHandling sh;
}

using namespace fast_planner;

int main(int argc, char** argv) {
  ros::init(argc, argv, "fast_planner_node");
  ros::NodeHandle nh("~");    //NodeHandle 对象 nh 用于与 ROS 系统进行交互，~ 表示私有命名空间。

  int planner;    //选择规划模式
  nh.param("planner_node/planner", planner, -1);
  //创建了两个 FSM 对象：topo_replan 和 kino_replan，分别代表拓扑规划和运动学规划的有限状态机。 在FSM中初始化了定时器 定时器会调用规划器
  TopoReplanFSM topo_replan;
  KinoReplanFSM kino_replan;

  if (planner == 1) {
    kino_replan.init(nh);
  } else if (planner == 2) {
    topo_replan.init(nh);
  }

  ros::Duration(1.0).sleep();
  ros::spin();

  return 0;
}
