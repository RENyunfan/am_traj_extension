#include<plan_manage/se3_planner.h>
using namespace fmt;
MavGlobalPlanner::MavGlobalPlanner(Config &conf, ros::NodeHandle &nh_)
    : config(conf), nh(nh_),
      visualization(config, nh)
{
    print(fg(color::yellow_green)," -- [PLANNER] MavFlobalPlanner init process begin.\n");
    // point_cloud_sub_ = nh.subscribe("PointCloud_in", 10, &MavGlobalPlanner::point_cloud_cb, this,ros::TransportHints().tcpNoDelay());
    _traj_pub = nh.advertise<quadrotor_msgs::PolynomialTrajectory>("trajectory", 50);
    marker_pub = nh.advertise<visualization_msgs::MarkerArray>("/marker", 50);
    occ_map.reset(new OccMap);
    occ_map->init(nh);
    fmt::print(fg(fmt::color::yellow_green), " -- [PLANNER] OccMap init success.\n");
    pathfinder.setParam(nh);
    pathfinder.setEnvironment(occ_map);
    nh.param("trajplanner/trajbudgetTime",traj_budget_time,0.1);
    // jps_pathfinder.setParam(nh);
    // jps_pathfinder.setMapUtil(map_util);
    fmt::print(fg(color::yellow_green)," -- [PLANNER] Astart searcher init success.\n");

}
MavGlobalPlanner::~MavGlobalPlanner()
{
    delete occ_map->obs_pointer;
}
// void MavGlobalPlanner::point_cloud_cb(const sensor_msgs::PointCloud2 & pointcloud_map){
//     if(has_map) return;
//     sensor_msgs::PointCloud cloud;
//     sensor_msgs::convertPointCloud2ToPointCloud(pointcloud_map,cloud);
//     cloud.header.frame_id = config.odomFrame;
//     obs_pointer  = new vec_Vec3f();
//     *obs_pointer = DecompROS::cloud_to_vec(cloud);
//     has_map = true;
//     ROS_WARN("POINT BE BUILD!!!!!!!!!!!");
// }
void MavGlobalPlanner::plan(const Eigen::MatrixXd &iniState, const Eigen::MatrixXd &finState,vector<Eigen::Vector3d>* wp_list){

    // if(!has_map) return;
    ros::Time t1,t2,t3,t4,t5,t6,t7,t8;
    t1 = ros::Time::now();
    Eigen::Vector3d zeroVec(0.0,0.0,0.0);
    Eigen::Vector3d start_pt;
    Eigen::Vector3d end_pt;
    start_pt = iniState.col(0);
    end_pt  = finState.col(0);
    vec_Vec3f path;
    if(!wp_list){
        pathfinder.search(start_pt,end_pt);
        path = pathfinder.getPathAndClear();

        //jps_pathfinder.publishAll();
    }else{
        vec_Vec3f path_0;
        pathfinder.search(start_pt,(*wp_list)[0]);
        path_0 = pathfinder.getPathAndClear();       
        for(int i = 0;i< wp_list->size()-1;i++){
            vec_Vec3f tmp_path;
            pathfinder.search((*wp_list)[i],(*wp_list)[i+1]);
            tmp_path = pathfinder.getPathAndClear();
            path_0.pop_back();
            path_0.insert(path_0.end(),tmp_path.begin(),tmp_path.end());
        }
        vec_Vec3f end_path;
        pathfinder.search((*wp_list)[wp_list->size()-1],end_pt);
        end_path = pathfinder.getPathAndClear();
        path_0.pop_back();
        path_0.insert(path_0.end(),end_path.begin(),end_path.end());
        path = path_0;
    }
    t2 = ros::Time::now();
    visualization.visualizePath(path);
    EllipsoidDecomp3D decomp_util;
    decomp_util.set_obs(*(occ_map->obs_pointer));
    decomp_util.set_local_bbox(Eigen::Vector3d(config.polyhedronBox(0),
                                               config.polyhedronBox(1),
                                               config.polyhedronBox(2)));
    //!!
    t3 = ros::Time::now();
    vec_E<Polyhedron3D> decompPolys;
    for(int i = 0;i<path.size()-1;){
        //find the farest unblocked point
        int k;
        for(k = i+1;k<path.size();k++){
            if(occ_map->isBlocked(path[i],path[k],100)||((path[i]-path[k]).norm()>=1.2)){
                k--;
                break;
            }
        }
        if(k<i+1){
            k = i+1;
        }
        //hzchzc
        if(k==path.size())
            k--;
        vec_Vec3f line;
        line.push_back(path[i]);
        line.push_back(path[k]);
        //hzchzc
        // visualization.visualizePath(line);
        // sleep(4);

        decomp_util.dilate(line);
        Polyhedron3D poly = decomp_util.get_polyhedrons()[0];
        decompPolys.push_back(poly);
        //find the nearest one to the boundry of poly.
        int j;
        for(j=k;j<path.size();j++){
            Vec3f pt;
            pt[0] = path[j][0];
            pt[1] = path[j][1];
            pt[2] = path[j][2];
            if(!poly.inside(pt)){
                break;
            }
        }
        j--;
        if(j>=path.size()-1){
            break;
        }
        int wp;
        wp = int((1*i+4*j)/5);
        if(wp<=i){
            i++;
        }
        else{
            i = wp;
        }
        
    }
//    std::cout << "Number of polyhedra from map: " << decompPolys.size() << std::endl;

    t4 = ros::Time::now();
    for (size_t i = 0; i < decompPolys.size(); i++)
    {
        decompPolys[i].add(Hyperplane3D(Eigen::Vector3d(0.0, 0.0, config.mapHeight),
                                          Eigen::Vector3d(0.0, 0.0, 1.0)));
        decompPolys[i].add(Hyperplane3D(Eigen::Vector3d(0.0, 0.0, config.mapFloor),
                                          Eigen::Vector3d(0.0, 0.0, -1.0)));
    }
    visualization.visualizeInitPoly(decompPolys);
    // vec_E<Polyhedron3D> debug_polys;
    // debug_polys.push_back(decompPolys[decompPolys.size()-1]);
    //delete redundant convex poly hzc

    double dt1 = ros::Time::now().toSec();
    vec_E<Polyhedron3D> Modipolys;
    for(size_t i=0;i<decompPolys.size();){
        /* add more intersection*/
        Polyhedron3D Adpoly = decompPolys[i];
        if(Modipolys.size()>0){
            Eigen::MatrixXd Modipolym;
            Eigen::MatrixXd Adpolym;
            Polyhedron3D Modipoly = Modipolys[Modipolys.size()-1];
            vec_E<Hyperplane3D> Modihyper = Modipoly.hyperplanes();
            Modipolym.resize(6, Modihyper.size());
            for (uint j = 0; j < Modihyper.size(); j++)
            {
                Modipolym.col(j) << Modihyper[j].n_, Modihyper[j].p_;
            }
            vec_E<Hyperplane3D> Adhyper = Adpoly.hyperplanes();
            Adpolym.resize(6, Adhyper.size());
            for (uint j = 0; j < Adhyper.size(); j++)
            {
                Adpolym.col(j) << Adhyper[j].n_, Adhyper[j].p_;
            }
            if(geoutils::calInterVol(Modipolym,Adpolym)<0.5){
                // ROS_WARN("add one!");
                // int min_index = inf,max_index = -1;
                // for(int i = 0; i<path.size();i++){
                //     if(Adpoly.inside(path[i])&&Modipoly.inside(path[i])){
                //         min_index = min_index>i?i:min_index;
                //         max_index = max_index>i?max_index:i;
                //     }
                // }
                // vec_Vec3f line;
                // line.push_back(path[min_index]);
                // line.push_back(path[max_index]);
                // decomp_util.dilate(line);
                // Polyhedron3D poly = decomp_util.get_polyhedrons()[0];
                // Modipolys.push_back(poly);
            }
        }


        Modipolys.push_back(Adpoly);
        if(i == decompPolys.size()-1) break;
        Eigen::MatrixXd current_poly;
        vec_E<Hyperplane3D> current_hyperplanes = decompPolys[i].hyperplanes();
        current_poly.resize(6, current_hyperplanes.size());
        for (uint j = 0; j < current_hyperplanes.size(); j++)
        {
            current_poly.col(j) << current_hyperplanes[j].n_, current_hyperplanes[j].p_;
            //outside
        }
        //find the farthest poly
        size_t far_index =  decompPolys.size()-1;
        for(size_t k = i+1;k<decompPolys.size();k++){
            Eigen::MatrixXd spare_poly;
            vec_E<Hyperplane3D> spare_hyperplanes = decompPolys[k].hyperplanes();
            spare_poly.resize(6, spare_hyperplanes.size());
            for (uint j = 0; j < spare_hyperplanes.size(); j++)
            {
                spare_poly.col(j) << spare_hyperplanes[j].n_, spare_hyperplanes[j].p_;
                //outside
            }
            if(geoutils::calInterVol(current_poly,spare_poly)>3.3){
                // ROS_WARN("1111111111111111111: %d %d",k,i);
                //inter more, delete less 
                continue;
            }
            else{
                k--;
                if(k<=i){
                    k = i+1;
                }
                far_index = k;
                break;
            }
        }
        i = far_index;
    }
    decompPolys = Modipolys;
    double dt2 = ros::Time::now().toSec();


    visualization.visualizePolyH(decompPolys);   
    t5 = ros::Time::now();
    std::vector<Eigen::MatrixXd> hPolys;
    Eigen::MatrixXd current_poly;

    for (uint i = 0; i < decompPolys.size(); i++)
    {
        vec_E<Hyperplane3D> current_hyperplanes = decompPolys[i].hyperplanes();
        current_poly.resize(6, current_hyperplanes.size());
        for (uint j = 0; j < current_hyperplanes.size(); j++)
        {
            current_poly.col(j) << current_hyperplanes[j].n_, current_hyperplanes[j].p_;
            //outside
        }
        hPolys.push_back(current_poly);
    }
//    std::cout << "Number of polyhedra from map: " << hPolys.size() << std::endl;
//    std::cout<<"-----------------------------------------------\n";
//    for(int i = 0;i<hPolys.size();i++){
//        std::cout<<"poly id: "<<i<<endl;
//        std::cout<<hPolys[i]<<std::endl;
//    }
//    std::cout<<"--------------------------------------------\n";
    t6 = ros::Time::now();

    std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now();
    if (!nonlinOpt.setup(config.rho, config.totalT, iniState, finState, hPolys, INFINITY,
                            config.qdIntervals, config.horizHalfLen, config.vertHalfLen,
                            config.safeMargin, config.velMax, config.thrustAccMin, config.thrustAccMax,
                            config.bodyRateMax, config.gravAcc, config.penaltyPVTB, config.useC2Diffeo))
    {
        return;
    }
    double finalObj = nonlinOpt.optimize(traj, config.optRelTol);
    std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now();

    VisualUtils::PublishTrajJerkToMarkerArray(marker_pub,traj);
    VisualUtils::PublishTrajVelToMarkerArray(marker_pub,traj);
    double compTime = std::chrono::duration_cast<std::chrono::microseconds>(toc - tic).count() * 1.0e-3;
    // nonlinOpt.kill_kernel();
    t7 = ros::Time::now();
    fmt::print(fg(fmt::color::green), " -- [PLANNER] Plan success!\n");
    fmt::print(fg(fmt::color::sky_blue), "\t Final traj cost: {0:.4}\n", finalObj);
    fmt::print(fg(fmt::color::sky_blue), "\t Maximum Vel: {0:.4} m/s \n",  traj.getMaxVelRate());
    fmt::print(fg(fmt::color::sky_blue), "\t Maximum Acc: {0:.4} m/s^2 \n", traj.getMaxAccRate());
    fmt::print(fg(fmt::color::sky_blue), "\t Total duration: {0:.4} s\n", traj.getTotalDuration());
    fmt::print(fg(fmt::color::light_golden_rod_yellow), "\t Poly num: {}\n", hPolys.size()  );
    fmt::print(fg(fmt::color::light_golden_rod_yellow), "\t GCOPTER: {0:.4} ms \n", compTime);
    fmt::print(fg(fmt::color::light_golden_rod_yellow), "\t A* search: {0:.4} ms\n", (t2-t1).toSec()*1000  );
    fmt::print(fg(fmt::color::light_golden_rod_yellow), "\t SFC: {0:.4} ms\n", (t6-t2).toSec()*1000  );
    fmt::print(fg(fmt::color::light_golden_rod_yellow), "\t Delete SFC time: {0:.4} ms\n", (dt2-dt1)*1000  );
    fmt::print(fg(fmt::color::gold), "\t Total time: {0:.4} ms\n", compTime + ((t2-t1).toSec() + (t6-t2).toSec()) * 1000  );


    if (traj.getPieceNum() > 0)
    {
        lastTraj = traj;
        lastIniStamp = ros::Time::now();

        visualization.visualize(traj);
        visualization.visualizeEllipsoid(traj,100);
        // visualization.visualizeQuadrotor(traj, 70);
    }
    quadrotor_msgs::PolynomialTrajectory traj_msg;
    traj_msg = traj2msg(traj);
    _traj_pub.publish(traj_msg);

//    t8 = ros::Time::now();
//    std::cout<<"start: "<<start_pt.transpose()<<std::endl;
//    std::cout<<"end: "<<end_pt.transpose()<<std::endl;
//    std::cout<<"t2-t1: "<<(t2-t1).toSec()<<" t3-t2: "<<(t3-t2).toSec()<<"\n";
//    std::cout<<"t4-t3: "<<(t4-t3).toSec()<<" t5-t4: "<<(t5-t4).toSec()<<std::endl;
//    std::cout<<"t6-t5: "<<(t6-t5).toSec()<<" t7-t6: "<<(t7-t6).toSec()<<"\n";
//    std::cout<<"t8-t7: "<<(t8-t7).toSec()<<std::endl;
}

quadrotor_msgs::PolynomialTrajectory MavGlobalPlanner::traj2msg(Trajectory traj){
    static int count=0;
    quadrotor_msgs::PolynomialTrajectory traj_msg;
    traj_msg.header.seq = count;
    traj_msg.header.stamp = ros::Time::now();
    traj_msg.header.frame_id = config.odomFrame;
    traj_msg.trajectory_id = count;
    traj_msg.action = quadrotor_msgs::PolynomialTrajectory::ACTION_ADD;
    traj_msg.num_order = traj[0].getOrder(); // the order of polynomial
    traj_msg.num_segment = traj.getPieceNum();
    traj_msg.start_yaw = 0;
    traj_msg.final_yaw = 0;
    for(unsigned int i=0; i<traj_msg.num_segment; i++)
    {
        for (unsigned int j = 0; j <= traj[i].getOrder(); j++)
        {
          CoefficientMat coemat = traj[i].normalizePosCoeffMat();
          traj_msg.coef_x.push_back(coemat(0,j));
          traj_msg.coef_y.push_back(coemat(1,j));
          traj_msg.coef_z.push_back(coemat(2,j));
        }
        traj_msg.time.push_back(traj[i].getDuration());
        traj_msg.order.push_back(traj[i].getOrder());
    }
    traj_msg.mag_coeff = 1;
    count++;
    return traj_msg;
}
bool MavGlobalPlanner::trajplan(const Eigen::MatrixXd &iniState, const Eigen::MatrixXd &finState,Trajectory& optitraj,double & opticost,vec_Vec3f* ref_path,bool forcepub ){
    Eigen::Vector3d zeroVec(0.0,0.0,0.0);
    Eigen::Vector3d start_pt;
    Eigen::Vector3d end_pt;
    start_pt = iniState.col(0);
    end_pt  = finState.col(0); 
    vec_Vec3f path;
    /*
       
    if(!pathfinder.search(start_pt,end_pt,forcepub)){
        ROS_WARN("astar failed to find a path");
        pathfinder.reset();
        return false;
    }
    path = pathfinder.getPathAndClear();
    visualization.visualizePath(path);
    */
    if(ref_path==NULL){
        if(!pathfinder.search(start_pt,end_pt,forcepub)){
            ROS_WARN("astar failed to find a path");
            pathfinder.reset();
            return false;
        }
        path = pathfinder.getPathAndClear();
    }   
    else{
        path = *ref_path;
    }
    //line of sight
    



    visualization.visualizePath(path);
    EllipsoidDecomp3D decomp_util;
    decomp_util.set_obs(*(occ_map->obs_pointer));
    decomp_util.set_local_bbox(Eigen::Vector3d(config.polyhedronBox(0),
                                               config.polyhedronBox(1),
                                               config.polyhedronBox(2)));
    vec_E<Polyhedron3D> decompPolys;
    double sum = 0;
    vec_Vec3f refined_path;
    for(int i = 0;i<path.size()-1;){
        int k;
        double tm1 = ros::Time::now().toSec();
        for(k = i+1;k<path.size();k++){
            if(occ_map->isBlocked(path[i],path[k],100)||((path[i]-path[k]).norm()>=1.2)){
                k--;
                break;
            }
        }
        double tm2 = ros::Time::now().toSec();
        sum+=tm2-tm1;
        if(k<i+1){
            k = i+1;
        }
        //hzchzc
        if(k==path.size())
            k--;
        
        vec_Vec3f line;
        line.push_back(path[i]);
        line.push_back(path[k]);
        decomp_util.dilate(line);
        Polyhedron3D poly = decomp_util.get_polyhedrons()[0];
        decompPolys.push_back(poly);
        //find the nearest one to the boundry of poly.
        int j;
        for(j=k;j<path.size();j++){
            Vec3f pt;
            pt[0] = path[j][0];
            pt[1] = path[j][1];
            pt[2] = path[j][2];
            if(!poly.inside(pt)){
                break;
            }
        }
        j--;
        
        if(j>=path.size()-1){
            break;
        }
        //j is the end of the current line, while i is the front. 
        int wp;
        wp = int((1*i+4*j)/5);
        // wp = j;
        if(wp<=i){
            i++;
        }
        else{
            i = wp;
        }
        
        
    }
    for (size_t i = 0; i < decompPolys.size(); i++)
    {
        decompPolys[i].add(Hyperplane3D(Eigen::Vector3d(0.0, 0.0, config.mapHeight),
                                          Eigen::Vector3d(0.0, 0.0, 1.0)));
        decompPolys[i].add(Hyperplane3D(Eigen::Vector3d(0.0, 0.0, config.mapFloor),
                                          Eigen::Vector3d(0.0, 0.0, -1.0)));
    }
    visualization.visualizeInitPoly(decompPolys);
    //delete polys
    vec_E<Polyhedron3D> Modipolys;
    for(size_t i=0;i<decompPolys.size();){
        Polyhedron3D Adpoly = decompPolys[i];
        Modipolys.push_back(Adpoly);
        if(i == decompPolys.size()-1) break;
        Eigen::MatrixXd current_poly;
        vec_E<Hyperplane3D> current_hyperplanes = decompPolys[i].hyperplanes();
        current_poly.resize(6, current_hyperplanes.size());
        for (uint j = 0; j < current_hyperplanes.size(); j++)
        {
            current_poly.col(j) << current_hyperplanes[j].n_, current_hyperplanes[j].p_;
            //outside
        }
        //find the farthest poly
        size_t far_index =  decompPolys.size()-1;
        for(size_t k = i+1;k<decompPolys.size();k++){
            Eigen::MatrixXd spare_poly;
            vec_E<Hyperplane3D> spare_hyperplanes = decompPolys[k].hyperplanes();
            spare_poly.resize(6, spare_hyperplanes.size());
            for (uint j = 0; j < spare_hyperplanes.size(); j++)
            {
                spare_poly.col(j) << spare_hyperplanes[j].n_, spare_hyperplanes[j].p_;
                //outside
            }
            if(geoutils::calInterVol(current_poly,spare_poly)>0.25){
                //2.0
                continue;
            }
            else{
                k--;
                if(k<=i){
                    k = i+1;
                }
                far_index = k;
                break;
            }
        }
        i = far_index;
    }
    decompPolys = Modipolys;
    std::vector<Eigen::MatrixXd> hPolys;
    Eigen::MatrixXd current_poly;
    for (uint i = 0; i < decompPolys.size(); i++)
    {
        vec_E<Hyperplane3D> current_hyperplanes = decompPolys[i].hyperplanes();
        current_poly.resize(6, current_hyperplanes.size());
        for (uint j = 0; j < current_hyperplanes.size(); j++)
        {
            current_poly.col(j) << current_hyperplanes[j].n_, current_hyperplanes[j].p_;
            //outside
        }
        hPolys.push_back(current_poly);
    }
    std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now();
    if(!forcepub){
        if (!nonlinOpt.setup(config.rho, config.totalT, iniState, finState, hPolys, INFINITY,
                                config.qdIntervals, config.horizHalfLen, config.vertHalfLen,
                                config.safeMargin, config.velMax, config.thrustAccMin, config.thrustAccMax,
                                config.bodyRateMax, config.gravAcc, config.penaltyPVTB, config.useC2Diffeo,traj_budget_time))
        {
            ROS_WARN("traj generate failed");
            return false;
        }
    }
    else{
        if (!nonlinOpt.setup(config.rho, config.totalT, iniState, finState, hPolys, INFINITY,
                                config.qdIntervals, config.horizHalfLen, config.vertHalfLen,
                                config.safeMargin, config.velMax, config.thrustAccMin, config.thrustAccMax,
                                config.bodyRateMax, config.gravAcc, config.penaltyPVTB, config.useC2Diffeo))
        {
            ROS_WARN("traj generate failed");
            return false;
        }
    }
    double finalObj = nonlinOpt.optimize(traj, config.optRelTol);
    if(finalObj<0){
        ROS_WARN("traj planner exceed budget time");
        return false;
    }
    std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now();
    double compTime = std::chrono::duration_cast<std::chrono::microseconds>(toc - tic).count() * 1.0e-3;
    std::cout << "Optimization time usage: " << compTime << " ms" << std::endl;
    std::cout << "Final cost: " << finalObj << std::endl;
    std::cout << "Maximum Vel: " << traj.getMaxVelRate() << std::endl;
    std::cout << "Maximum Acc: " << traj.getMaxAccRate() << std::endl;
    std::cout << "Total Duation: " << traj.getTotalDuration() << std::endl;
    opticost = finalObj;
    visualization.visualizePolyH(decompPolys);   
    // visualization.visualize(traj);
    // visualization.visualizeEllipsoid(traj,100);
    optitraj = traj;
    return true;

}

Visualization::Visualization(Config &conf, ros::NodeHandle &nh_)
    : config(conf), nh(nh_)
{
    trajPub = nh.advertise<visualization_msgs::Marker>("/visualization/trajectory", 1);
    ellipsoidPub = nh.advertise<visualization_msgs::MarkerArray>("/visualization/ellipsoid", 1);
    quadrotorPub = nh.advertise<visualization_msgs::MarkerArray>("/visualization/quadrotor", 1);
    hPolyPub = nh.advertise<decomp_ros_msgs::PolyhedronArray>("/visualization/polyhedra", 1);
    tiltRatePub = nh.advertise<std_msgs::Float64>("/visualization/tilt_rate", 1);
    thrMagPub = nh.advertise<std_msgs::Float64>("/visualization/thrust_magnitute", 1);
    pathPub  = nh.advertise<visualization_msgs::Marker>("/visualization/path", 1);
    initPolyPub = nh.advertise<decomp_ros_msgs::PolyhedronArray>("/visualization/init_polyhedra", 1);
    cosEliPub = nh.advertise<visualization_msgs::MarkerArray>("/visualization/cos_ellipsoid", 1);
    DebugPub = nh.advertise<visualization_msgs::Marker>("/visualization/Debug_path", 1);
}

void Visualization::visualizeInitPoly(const vec_E<Polyhedron3D> &polyhedra){
    decomp_ros_msgs::PolyhedronArray poly_msg = DecompROS::polyhedron_array_to_ros(polyhedra);
    poly_msg.header.frame_id = config.odomFrame;
    poly_msg.header.stamp = ros::Time::now();
    initPolyPub.publish(poly_msg);
}

void Visualization::visualize(const Trajectory &traj, const int samples)
{
    visualization_msgs::Marker trajMarker;

    trajMarker.id = 0;
    trajMarker.header.stamp = ros::Time::now();
    trajMarker.header.frame_id = config.odomFrame;
    trajMarker.pose.orientation.w = 1.00;
    trajMarker.action = visualization_msgs::Marker::ADD;
    trajMarker.type = visualization_msgs::Marker::LINE_STRIP;
    trajMarker.header.frame_id = config.odomFrame;
    trajMarker.ns = "trajectory";
    trajMarker.color.r = config.trajVizRGB(0);
    trajMarker.color.g = config.trajVizRGB(1);
    trajMarker.color.b = config.trajVizRGB(2);
    trajMarker.color.a = 1.00;
    trajMarker.scale.x = config.trajVizWidth;

    double dt = traj.getTotalDuration() / samples;
    geometry_msgs::Point point;
    Eigen::Vector3d pos;
    for (int i = 0; i <= samples; i++)
    {
        pos = traj.getPos(dt * i);
        point.x = pos(0);
        point.y = pos(1);
        point.z = pos(2);
        trajMarker.points.push_back(point);
    }

    trajPub.publish(trajMarker);
}
void Visualization::visualizePath(const vec_Vec3f& path){
    visualization_msgs::Marker trajMarker;

    trajMarker.id = 0;
    trajMarker.header.stamp = ros::Time::now();
    trajMarker.header.frame_id = config.odomFrame;
    trajMarker.pose.orientation.w = 1.00;
    trajMarker.pose.orientation.x = 0.0;
    trajMarker.pose.orientation.y = 0.0;
    trajMarker.pose.orientation.z = 0.0;
    trajMarker.action = visualization_msgs::Marker::ADD;
    trajMarker.type = visualization_msgs::Marker::CUBE_LIST;
    trajMarker.ns = "path";
    trajMarker.color.r = 1.0;
    trajMarker.color.g = 0;
    trajMarker.color.b = 0;
    trajMarker.color.a = 1.00;
    trajMarker.scale.x = 0.1;
    trajMarker.scale.y = 0.1;
    trajMarker.scale.z = 0.1;


    geometry_msgs::Point point;
    Eigen::Vector3d pos;
    for (int i = 0; i < path.size(); i++)
    {
        pos = path[i];
        point.x = pos(0);
        point.y = pos(1);
        point.z = pos(2);
        trajMarker.points.push_back(point);
    }

    pathPub.publish(trajMarker);
}
void Visualization::visualizeDebugPath(const vec_Vec3f& path){
    visualization_msgs::Marker trajMarker;

    trajMarker.id = 0;
    trajMarker.header.stamp = ros::Time::now();
    trajMarker.header.frame_id = config.odomFrame;
    trajMarker.pose.orientation.w = 1.00;
    trajMarker.pose.orientation.x = 0.0;
    trajMarker.pose.orientation.y = 0.0;
    trajMarker.pose.orientation.z = 0.0;
    trajMarker.action = visualization_msgs::Marker::ADD;
    trajMarker.type = visualization_msgs::Marker::CUBE_LIST;

    trajMarker.ns = "debug_path";
    trajMarker.color.r = 0.0;
    trajMarker.color.g = 0;
    trajMarker.color.b = 1.0;
    trajMarker.color.a = 1.00;
    trajMarker.scale.x = 0.05;
    trajMarker.scale.y = 0.05;
    trajMarker.scale.z = 0.05;


    geometry_msgs::Point point;
    Eigen::Vector3d pos;
    // Eigen::Vector3d
    for (int i = 0; i < path.size(); i++)
    {
        pos = path[i];
        point.x = pos(0);
        point.y = pos(1);
        point.z = pos(2);
        trajMarker.points.push_back(point);
    }

    DebugPub.publish(trajMarker);
}
void Visualization::visualizeEllipsoid(const Trajectory &traj, const int samples)
{
    visualization_msgs::Marker ellipsoidMarker;
    visualization_msgs::MarkerArray ellipsoidMarkers;

    ellipsoidMarker.id = 0;
    ellipsoidMarker.type = visualization_msgs::Marker::MESH_RESOURCE;
    ellipsoidMarker.mesh_resource = config.ellipsoidPath;
    ellipsoidMarker.header.stamp = ros::Time::now();
    ellipsoidMarker.header.frame_id = config.odomFrame;
    ellipsoidMarker.pose.orientation.w = 1.00;
    ellipsoidMarker.action = visualization_msgs::Marker::ADD;
    ellipsoidMarker.ns = "ellipsoids";
    ellipsoidMarker.color.r = config.ellipsoidVizRGBA(0);
    ellipsoidMarker.color.g = config.ellipsoidVizRGBA(1);
    ellipsoidMarker.color.b = config.ellipsoidVizRGBA(2);
    ellipsoidMarker.color.a = config.ellipsoidVizRGBA(3);
    ellipsoidMarker.scale.x = config.horizHalfLen * 2.0;
    ellipsoidMarker.scale.y = config.horizHalfLen * 2.0;
    ellipsoidMarker.scale.z = config.vertHalfLen * 2.0;

    ellipsoidMarker.action = visualization_msgs::Marker::DELETEALL;
    ellipsoidMarkers.markers.push_back(ellipsoidMarker);
    ellipsoidPub.publish(ellipsoidMarkers);
    ellipsoidMarker.action = visualization_msgs::Marker::ADD;
    ellipsoidMarkers.markers.clear();

    double dt = traj.getTotalDuration() / samples;
    geometry_msgs::Point point;
    Eigen::Vector3d pos;
    Eigen::Matrix3d rotM;
    Eigen::Quaterniond quat;
    for (int i = 0; i <= samples; i++)
    {
        pos = traj.getPos(dt * i);
        ellipsoidMarker.pose.position.x = pos(0);
        ellipsoidMarker.pose.position.y = pos(1);
        ellipsoidMarker.pose.position.z = pos(2);
        traj.getRotation(dt * i, 0.0, config.gravAcc, rotM);
        quat = Eigen::Quaterniond(rotM);
        ellipsoidMarker.pose.orientation.w = quat.w();
        ellipsoidMarker.pose.orientation.x = quat.x();
        ellipsoidMarker.pose.orientation.y = quat.y();
        ellipsoidMarker.pose.orientation.z = quat.z();
        ellipsoidMarkers.markers.push_back(ellipsoidMarker);
        ellipsoidMarker.id++;
    }

    ellipsoidPub.publish(ellipsoidMarkers);
}
void Visualization::visualizeCosEllipsoid(std::vector<Eigen::Vector3d> poslist,std::vector<Eigen::Matrix3d> rotmlist){
    visualization_msgs::Marker ellipsoidMarker;
    visualization_msgs::MarkerArray ellipsoidMarkers;

    ellipsoidMarker.id = 0;
    ellipsoidMarker.type = visualization_msgs::Marker::MESH_RESOURCE;
    ellipsoidMarker.mesh_resource = config.ellipsoidPath;
    ellipsoidMarker.header.stamp = ros::Time::now();
    ellipsoidMarker.header.frame_id = config.odomFrame;
    ellipsoidMarker.pose.orientation.w = 1.00;
    ellipsoidMarker.action = visualization_msgs::Marker::ADD;
    ellipsoidMarker.ns = "cos_ellipsoids";
    ellipsoidMarker.color.r = 0.0;
    ellipsoidMarker.color.g = 0.0;
    ellipsoidMarker.color.b = 0.0;
    ellipsoidMarker.color.a = 0.8;
    ellipsoidMarker.scale.x = config.horizHalfLen * 2.0;
    ellipsoidMarker.scale.y = config.horizHalfLen * 2.0;
    ellipsoidMarker.scale.z = config.vertHalfLen * 2.0;
    // ellipsoidMarker.scale.x = 1;
    // ellipsoidMarker.scale.y = 1;
    ellipsoidMarker.scale.z = config.vertHalfLen * 2.0;



    ellipsoidMarker.action = visualization_msgs::Marker::DELETEALL;
    ellipsoidMarkers.markers.push_back(ellipsoidMarker);
    cosEliPub.publish(ellipsoidMarkers);
    ellipsoidMarker.action = visualization_msgs::Marker::ADD;
    ellipsoidMarkers.markers.clear();

    geometry_msgs::Point point;
    Eigen::Vector3d pos;
    Eigen::Matrix3d rotM;
    Eigen::Quaterniond quat;
    for (int i = 0; i < poslist.size(); i++)
    {
        pos = poslist[i];
        ellipsoidMarker.pose.position.x = pos(0);
        ellipsoidMarker.pose.position.y = pos(1);
        ellipsoidMarker.pose.position.z = pos(2);
        rotM = rotmlist[i];
        quat = Eigen::Quaterniond(rotM);
        ellipsoidMarker.pose.orientation.w = quat.w();
        ellipsoidMarker.pose.orientation.x = quat.x();
        ellipsoidMarker.pose.orientation.y = quat.y();
        ellipsoidMarker.pose.orientation.z = quat.z();
        ellipsoidMarkers.markers.push_back(ellipsoidMarker);
        ellipsoidMarker.id++;
    }

    cosEliPub.publish(ellipsoidMarkers);
}
void Visualization::visualizeQuadrotor(const Trajectory &traj, const int samples)
{
    visualization_msgs::Marker quadrotorMarker;
    visualization_msgs::MarkerArray quadrotorMarkers;

    quadrotorMarker.id = 0;
    quadrotorMarker.type = visualization_msgs::Marker::MESH_RESOURCE;
    quadrotorMarker.mesh_use_embedded_materials = true;
    quadrotorMarker.mesh_resource = config.quadrotorPath;
    quadrotorMarker.header.stamp = ros::Time::now();
    quadrotorMarker.header.frame_id = config.odomFrame;
    quadrotorMarker.pose.orientation.w = 1.00;
    quadrotorMarker.action = visualization_msgs::Marker::ADD;
    quadrotorMarker.ns = "quadrotor";
    quadrotorMarker.color.r = 0.0;
    quadrotorMarker.color.g = 0.0;
    quadrotorMarker.color.b = 0.0;
    quadrotorMarker.color.a = 0.0;
    quadrotorMarker.scale.x = (config.horizHalfLen+config.safeMargin*3.3+0.2) * sqrt(2.0);
    quadrotorMarker.scale.y = (config.horizHalfLen+config.safeMargin*3.3+0.2) * sqrt(2.0);
    quadrotorMarker.scale.z = config.vertHalfLen * 8.0;

    quadrotorMarker.action = visualization_msgs::Marker::DELETEALL;
    quadrotorMarkers.markers.push_back(quadrotorMarker);
    quadrotorPub.publish(quadrotorMarkers);
    quadrotorMarker.action = visualization_msgs::Marker::ADD;
    quadrotorMarkers.markers.clear();

    // double dt = traj.getTotalDuration() / samples;
    //hzc30
    double dt = 5.0/samples;
    geometry_msgs::Point point;
    Eigen::Vector3d pos;
    Eigen::Matrix3d rotM;
    Eigen::Quaterniond quat;
    for (int i = 0; i <= samples; i++)
    {
        pos = traj.getPos(6+dt * i);//4.5
        quadrotorMarker.pose.position.x = pos(0);
        quadrotorMarker.pose.position.y = pos(1);
        quadrotorMarker.pose.position.z = pos(2);
        // traj.getRotation(6+dt * i, M_PI_4, config.gravAcc, rotM);hzc
        traj.getRotation(6+dt * i, M_PI_2, config.gravAcc, rotM);
        quat = Eigen::Quaterniond(rotM);
        quadrotorMarker.pose.orientation.w = quat.w();
        quadrotorMarker.pose.orientation.x = quat.x();
        quadrotorMarker.pose.orientation.y = quat.y();
        quadrotorMarker.pose.orientation.z = quat.z();
        quadrotorMarkers.markers.push_back(quadrotorMarker);
        quadrotorMarker.id++;
    }

    quadrotorPub.publish(quadrotorMarkers);
}

void Visualization::visualizePolyH(const vec_E<Polyhedron3D> &polyhedra)
{
    decomp_ros_msgs::PolyhedronArray poly_msg = DecompROS::polyhedron_array_to_ros(polyhedra);
    poly_msg.header.frame_id = config.odomFrame;
    poly_msg.header.stamp = ros::Time::now();
    hPolyPub.publish(poly_msg);
}

void Visualization::visualizeProfile(const Trajectory &traj, const double &t)
{
    Eigen::Vector3d tiltRate = traj.getTiltRate(t, config.gravAcc);
    Eigen::Vector3d thrAcc = traj.getAcc(t);
    thrAcc(2) += config.gravAcc;

    std_msgs::Float64 tilt, thr;
    tilt.data = tiltRate.norm();
    thr.data = thrAcc.norm();

    tiltRatePub.publish(tilt);
    thrMagPub.publish(thr);
}
