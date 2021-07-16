# Constrained BVP Test

# 1 运行时间分析

目前BVP有三个类型：

* 第一种是给定时间（或者启发式时间的bvp），计算速度为平均`0.45us`一次。
* 第二种是最优时间分配，计算速度平均$1.2us$一次
* 第三种是基于Bisection的约束最优生成算法，目前速度最慢，为$12us$一次。

分析约束BVP算法时间

| 函数                   | 时间   | 功能                             |
| ---------------------- | ------ | -------------------------------- |
| EnforceFeasibility     | 可忽略 | 起点和重点状态的饱和函数         |
| GenFixStateMinSnapOptT | 1.6us  | 计算global minimum               |
| CheckPiece             | 0.5us  | 检测轨迹速度和加速度是否超出上限 |
| GenFixStateMinSnap     | 0.45us | 生成最优控制                     |
| EvaluateSnapCost       | 0.26us | 评估最终轨迹的cost               |

时间消耗基本为3：2：1

* 3： 生成最优时间轨迹+1次Check【但是理论上只需要2us】
* 2：搜索可行解【两次Check+bvp】【理论上2us】
* 2：Bisection跑【5次check+bvp】【理论上5us】

```cpp
BVP with t time consuming  0.472259 us
BVP + Opt t time consuming  1.209608 us
CBVP time consuming  11.362139 us
```

