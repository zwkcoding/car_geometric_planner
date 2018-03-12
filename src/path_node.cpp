#include <ros/ros.h>
#include <nav_msgs/Path.h>
#include <tf/transform_datatypes.h>
#include <geometry_msgs/Pose.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <geometry_msgs/PoseStamped.h>

#include <chrono>
#include "r_s_planner/reeds_shepp.h"
#include "r_s_planner/dubins.h"

double start_state[3] = {0, 0, M_PI/4};
double final_state[3] = {3, 4, 0};
bool goal_set, start_set;

void startCb(const geometry_msgs::PoseWithCovarianceStampedConstPtr &msg) {
    start_state[0] = msg->pose.pose.position.x;
    start_state[1] = msg->pose.pose.position.y;
    start_state[2] = tf::getYaw(msg->pose.pose.orientation);// [-pi,pi]
    std::cout << "get the goal." << std::endl;
    start_set = true;
}

void goalCb(const geometry_msgs::PoseStampedConstPtr &goal) {
    final_state[0] = goal->pose.position.x;
    final_state[1] = goal->pose.position.y;
    final_state[2] = tf::getYaw(goal->pose.orientation);// [-pi,pi]
    std::cout << "get the goal." << std::endl;
    goal_set = true;
}
int main(int argc, char **argv) {
    // Initialize the node, publishers and subscribers.
    ros::init(argc, argv, "r_s_planner_node");
    ros::NodeHandle nh("~");
    // Create publishers and subscribers
    ros::Publisher r_s_pub = nh.advertise<nav_msgs::Path>("rs_path", 1, true);
    ros::Publisher db_pub = nh.advertise<nav_msgs::Path>("db_path", 1, true);
    ros::Subscriber start_sub = nh.subscribe("/initialpose", 1, startCb);
    ros::Subscriber end_sub = nh.subscribe("/move_base_simple/goal", 1, goalCb);
    double steer_radius = 5.8;
    double step_size = 0.5;
    ReedsSheppStateSpace rs_planner(steer_radius);
    DubinsStateSpace db_planner(steer_radius);
    start_set = goal_set = false;
    // Publisher in a loop.
    ros::Rate rate(10.0);
    while (nh.ok()) {
        ros::spinOnce();
//        if(!start_set || !goal_set) {
//            ROS_INFO("Waiting for input!");
//            rate.sleep();
//            continue;
//        }
        // Add data to grid map.
        auto start_time = std::chrono::steady_clock::now();
        std::vector<std::vector<double> > rs_path;
        double length = -1;
        rs_planner.sample(start_state, final_state, step_size, length, rs_path);
        auto end_time = std::chrono::steady_clock::now();
        double elapsed_msecondes = std::chrono::duration_cast<std::chrono::microseconds>
                                           (end_time - start_time).count() / 1000.0;
        ROS_INFO("R_S Path generation time:  %f msec.\n Path lenght : %f m",
                          elapsed_msecondes, length);

        nav_msgs::Path rs_path_msg;
        geometry_msgs::PoseStamped rs_pose;
        rs_path_msg.header.frame_id = "map";
        rs_path_msg.header.stamp = ros::Time::now();
        rs_pose.header = rs_path_msg.header;
        for (auto &point_itr : rs_path) {
            rs_pose.pose.position.x = point_itr[0];
            rs_pose.pose.position.y = point_itr[1];
            rs_pose.pose.position.z = 0;
            rs_pose.pose.orientation = tf::createQuaternionMsgFromYaw(point_itr[2]);
            rs_path_msg.poses.push_back(rs_pose);
        }
        r_s_pub.publish(rs_path_msg);

        std::vector<std::vector<double> > db_path;
        auto st = std::chrono::steady_clock::now();
        db_planner.sample(start_state, final_state, step_size, length, db_path);
        auto et = std::chrono::steady_clock::now();
        elapsed_msecondes = std::chrono::duration_cast<std::chrono::microseconds>
                                           (et - st).count() / 1000.0;
        ROS_INFO("DB Path generation time:  %f msec.\n Path lenght : %f m",
                 elapsed_msecondes, length);

        nav_msgs::Path db_path_msg;
        geometry_msgs::PoseStamped db_pose;
        db_path_msg.header.frame_id = "map";
        db_path_msg.header.stamp = ros::Time::now();
        db_pose.header = db_path_msg.header;
        for (auto &point_itr : db_path) {
            db_pose.pose.position.x = point_itr[0];
            db_pose.pose.position.y = point_itr[1];
            db_pose.pose.position.z = 0;
            db_pose.pose.orientation = tf::createQuaternionMsgFromYaw(point_itr[2]);
            db_path_msg.poses.push_back(db_pose);
        }
        db_pub.publish(db_path_msg);

        rate.sleep();
    }

    return 0;
}
