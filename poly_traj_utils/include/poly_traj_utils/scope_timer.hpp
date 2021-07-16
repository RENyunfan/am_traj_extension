//
// Created by yunfan on 2021/3/19.
//


#ifndef SCROPE_TIMER_HPP//SRC_POLY_VISUAL_UTILS_HPP
#define SCROPE_TIMER_HPP
#include <chrono>
class TimeConsuming{
public:
    TimeConsuming();
    TimeConsuming(string msg, int repeat_time){
        repeat_time_ =repeat_time;
        msg_ = msg;
        tc_start = std::chrono::high_resolution_clock::now();
        has_shown = false;
    }
    TimeConsuming(string msg){
        msg_ = msg;
        repeat_time_ = 1;
        tc_start = std::chrono::high_resolution_clock::now();
        has_shown = false;
    }
    ~TimeConsuming(){
        if(!has_shown){
            tc_end = std::chrono::high_resolution_clock::now();
            double dt = std::chrono::duration_cast<std::chrono::duration<double>>(tc_end - tc_start).count();
            double t_us = (double)dt*1e6/repeat_time_;
            printf("%s time consuming \033[32m %lf us\033[0m\n",msg_.c_str(),t_us);
        }
    }
    void start(){
        tc_start = std::chrono::high_resolution_clock::now();
    }

    void stop(){
        has_shown=true;
        tc_end = std::chrono::high_resolution_clock::now();
        double dt = std::chrono::duration_cast<std::chrono::duration<double>>(tc_end - tc_start).count();
//            ROS_WARN("%s time consuming %lf us.",msg_.c_str(),(double)(end_t - start_t).toNSec()/ 1e3);
        printf("%s time consuming \033[32m %lf us\033[0m\n",msg_.c_str(),(double)dt*1e6/repeat_time_);
    }

private:
    std::chrono::high_resolution_clock::time_point tc_start, tc_end;
    string msg_;
    int repeat_time_;
    bool has_shown = false;
};

#endif //SRC_POLY_VISUAL_UTILS_HPP
