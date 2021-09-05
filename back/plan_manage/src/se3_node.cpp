#include <plan_manage/se3_planner.h>
ros::Subscriber waypoints_sub,_odom_sub, goal_sub;
Eigen::Vector3d start_pos,start_vel,start_acc;
Eigen::Vector3d endpos,endvel;
Eigen::Vector3d odom_position;
Eigen::Vector3d odom_velocity;
int run_type;
Eigen::Quaterniond odom_q;
Eigen::Vector3d Zero3d(0,0,0);
Eigen::MatrixXd initState,finState;
std::vector<Eigen::Vector3d> gate_list;
Eigen::Vector3d target_pt;
void rcvWaypointsCallback( nav_msgs::Path  wp);
void goalCallback(const geometry_msgs::PoseStamped::ConstPtr & msg);
void odom_callback(const nav_msgs::Odometry::ConstPtr & odom);
MavGlobalPlanner* glbPlanner;
std::vector<Eigen::Vector3d> wp_list;
int main(int argc, char **argv)
{
    ros::init(argc, argv, "map_se3_node");
    start_pos<<0,0,1;
    start_vel<<0,0,0;
    start_acc<<0,0,0;
    endvel<<0,0,0;
    initState.resize(3,3);
    finState.resize(3,3);
    // wp_list.push_back(Eigen::Vector3d())



    Config config;
    ros::NodeHandle nh_priv("~");
    nh_priv.param("se3_node/run_type", run_type, -1);
    goal_sub = nh_priv.subscribe("/goal", 1, goalCallback);
//    waypoints_sub= nh_priv.subscribe( "waypoints", 1, rcvWaypointsCallback );
    _odom_sub = nh_priv.subscribe("odom", 1, odom_callback);
    config.loadParameters(ros::NodeHandle("~"));
    glbPlanner = new MavGlobalPlanner(config, nh_priv);
    ros::Rate rate(100);
    bool status = ros::ok();
    while(status){
        ros::spinOnce();
        status = ros::ok();
        rate.sleep();
    }
    return 0;
}

void goalCallback(const geometry_msgs::PoseStamped::ConstPtr & msg){
    target_pt = Eigen::Vector3d(msg->pose.position.x,msg->pose.position.y,1.0);
    finState <<target_pt,Zero3d,Zero3d;
    initState << odom_position,Zero3d,Zero3d;
    ros::Time t1 = ros::Time::now();
    glbPlanner->plan(initState,finState);
    // glbPlanner->plan(initState,finState,&wp_list);

    ros::Time t2 = ros::Time::now();
    std::cout<<"total plan time is "<<(t2-t1).toSec()*1e3<<" ms"<<std::endl;
}

void rcvWaypointsCallback( nav_msgs::Path  wp)                                                                                      
{  
    if( wp.poses[0].pose.position.z < 0.0 )
        return;

    ROS_INFO_STREAM("TARGET="<<Eigen::Vector3d(wp.poses[0].pose.position.x,wp.poses[0].pose.position.y,1.5));
    ROS_INFO("[node] receive the planning target");
    target_pt = Eigen::Vector3d(wp.poses[0].pose.position.x,wp.poses[0].pose.position.y,1.0);
    // initState <<start_pos,Zero3d,Zero3d;
    // glbPlanner->plan(initState,finState,&wp_list);
    // glbPlanner->plan(initState,finState,&gate_list);
    // glbPlanner->plan(initState,finState);
    // wp.poses[0].pose.position.x = 0.619061;
    // wp.poses[0].pose.position.y = 61.0349;
    std::vector<Eigen::Vector3d> wp_list;
    //默认 yaw = 0;
    Eigen::Vector3d base_p = odom_position;
    // wp_list.push_back(Eigen::Vector3d(base_p[0]+2,base_p[1]-2,1.5));
    // base_p = wp_list.back();
    // wp_list.push_back(Eigen::Vector3d(base_p[0]+2,base_p[1]+2,1.5));
    // base_p = wp_list.back();
    // wp_list.push_back(Eigen::Vector3d(base_p[0]+4,base_p[1]-4,1.5)); 
    // base_p = wp_list.back();
    // wp_list.push_back(Eigen::Vector3d(base_p[0]-4,base_p[1]-2,1.5));
    // base_p = wp_list.back();
    // wp_list.push_back(Eigen::Vector3d(base_p[0]-2,base_p[1]+2,1.5));
    // base_p = wp_list.back();
    // wp_list.push_back(Eigen::Vector3d(base_p[0]-1,base_p[1]+1,1.5));
    // base_p = wp_list.back();
    // wp_list.push_back(Eigen::Vector3d(base_p[0]-3,base_p[1]+3,1.5));
    // base_p = wp_list.back();
    // wp_list.push_back(Eigen::Vector3d(base_p[0]-1,base_p[1]+3,1.5));
    // base_p = wp_list.back();
    // wp_list.push_back(Eigen::Vector3d(base_p[0],base_p[1],1.5));
    // wp_list.push_back(Eigen::Vector3d(base_p[0],base_p[1],1.5));

    //
    double r = 2.0;
    for(double t=1.0/12.0;t<=2;t+=1.0/12.0){
        double x = r * cos(2*M_PI*t);
        double y = r * sin(2*M_PI*t);
        // r=x*(1+t);
        Eigen::Vector3d tmp  = odom_q.toRotationMatrix() * Eigen::Vector3d(base_p[0]+x,base_p[1]+y,1.5);
        wp_list.push_back(tmp);
        // if(intt/0.25)
        r = 1.0*int(t/0.25)+4.0;
        std::cout<<"r";
    }
    base_p = wp_list.back();

    // initState<<odom_position,Zero3d,Zero3d;
    // finState<<Eigen::Vector3d(wp.poses[0].pose.position.x,wp.poses[0].pose.position.y,1.5),Zero3d,Zero3d;

    // finState<<odom_q.toRotationMatrix()*Eigen::Vector3d(base_p[0]-1,base_p[1]+1,1.5),Zero3d,Zero3d;

//    initState<<Eigen::Vector3d(0,0,1),Zero3d,Zero3d;
//
    finState <<target_pt,Zero3d,Zero3d;
     initState<<odom_position,Zero3d,Zero3d;
    ros::Time t1 = ros::Time::now();
    glbPlanner->plan(initState,finState);
    // glbPlanner->plan(initState,finState,&wp_list);

    ros::Time t2 = ros::Time::now();
    std::cout<<"total plan time is "<<(t2-t1).toSec()*1e3<<" ms"<<std::endl;
    

}
void odom_callback(const nav_msgs::Odometry::ConstPtr& odom)
{
  if(run_type == 0 || run_type == 1)
  {
    const Eigen::Vector3d position(odom->pose.pose.position.x,
                                    odom->pose.pose.position.y,
                                    odom->pose.pose.position.z);
    const Eigen::Vector3d velocity(odom->twist.twist.linear.x,
                                    odom->twist.twist.linear.y,
                                    odom->twist.twist.linear.z);
    odom_position = position;
    odom_velocity = velocity;                               
    Eigen::Quaterniond q;
    q.w() = odom->pose.pose.orientation.w;
    q.x() = odom->pose.pose.orientation.x;
    q.y() = odom->pose.pose.orientation.y;  
    q.z() = odom->pose.pose.orientation.z;
    q = q.normalized();  
    odom_q = q;
  }
  else if (run_type == 2)
  {
    const Eigen::Vector3d position(odom->pose.pose.position.x,
                                    -odom->pose.pose.position.y,
                                    -odom->pose.pose.position.z);
    const Eigen::Vector3d velocity(odom->twist.twist.linear.x,
                                    -odom->twist.twist.linear.y,
                                    -odom->twist.twist.linear.z);
    odom_position = position;
    odom_velocity = velocity;                               
    Eigen::Quaterniond q;
    q.w() = odom->pose.pose.orientation.w;
    q.x() = odom->pose.pose.orientation.x;
    q.y() = odom->pose.pose.orientation.y;  
    q.z() = odom->pose.pose.orientation.z;
    Eigen::Matrix3d R1 ;
    R1 << 1,0,0,0,-1,0,0,0,-1;
    q = R1*q*R1;
    q = q.normalized();  
    odom_q = q;
  }


}
