[toc]

# Sample-based Motion Planning Algorithm

# 1 Overview

## Graph-search

对于图搜索问题，我们往往需要构建图。这种构建图的方法通常是以一定分辨率的离散空间，描述环境。很多方法保证分辨率完备（Resolution complete）和分辨率最优（Resolution Optimal）。意思是足够大的分辨率下，保证路径搜索的完备性和最优性。其中还有一些方法，比如A*算法是最优效率的。

图搜索也有一些缺点。

*  搜索结果的质量（例如连续性）受到分辨率的影响。
* 维度诅咒：随着维度的上升，图搜索算法的复杂度指数爆炸。

## Sampling-based 

与图搜索算法相比，则是基于采样的算法。基于采样的算法有这些有点

* 无须对工作空间进行离散化，而是从连续的空间进行随机采样。
* 算法复杂度随着维度的增加，增加没那么可怕。
* 很多算法是概率完备的（Probabilistically complete），意思就是说，随着采样的不断进行，一定能找到解。
* 渐进精度（Anytime resolution），意思是随着采样的进行，轨迹的精度和分辨率也不断上升。
* 有些算法是渐进最优的（Asymptotically optimal），例如RRT*，意思是随着采样的进行，最终轨迹会收敛到这个规划问题的最优轨迹。

当然基于采样的算法也有一些缺点，最大的问题就是计算复杂度比较高。

常见的基于采样的算法可以分为单查询（single-query）和多查询（Multiple-query）。

**单查询**

* 输入为一个起点和重点。
* 每一次采样一个随机点并进行相应的算法。
* 代表算法有RRT，RDT，RRT*，EST，SRT，RRM。

**多查询**

* 构建一个图（Road map）。
* 可以有很多的起点和重点状态对。（可以理解为构建好图之后可以在任意节点间进行搜索）。
* 代表算法有PRM，Lazy-PRM，dynamic PRM， PRM*。



## 发展史

* Dijkstra's Algorithm 1956
* A* 1968
* RRT 1998
* RRT* 2010
* CHOMP 2013 （基于梯度的方法）

# 1 FMT*

> Fast Marching Tree: a Fast Marching Sampling-Based Method for Optimal Motion Planning in Many Dimensions [Janson15]

快速扩展树。

FMT*是多查询和单查询的组合，采用懒动态规划（Lazy-Dynamic programming）

>In [programming language theory](https://en.wikipedia.org/wiki/Programming_language_theory), **lazy evaluation**, or **call-by-need**,[[1\]](https://en.wikipedia.org/wiki/Lazy_evaluation#cite_note-1) is an [evaluation strategy](https://en.wikipedia.org/wiki/Evaluation_strategy) which delays the evaluation of an [expression](https://en.wikipedia.org/wiki/Expression_(computer_science)) until its value is needed ([non-strict evaluation](https://en.wikipedia.org/wiki/Non-strict_evaluation)) 
>
>https://en.wikipedia.org/wiki/Lazy_evaluation

FMT*同样是渐进最优的，其收敛速度高于RRT\*和PRM\*。并且搜索路径是从期间增量式构建的。得益于扩展方式，他不需要重写父节点。

FMT*的优点有

* FMT*是比RRT\*更好的算法，因为在搜索树上的连接是接近最优的。
* 不需要重写父节点
* 比PRM*更好在于，FMT\*构建的是树状的地图，有很好的微分约束，逐渐的添加最优边直到连接到终点。

![image-20210712142533166](Sample-based%20Motion%20Planning%20Algorithm.assets/image-20210712142533166.png)

首先是Multiple-query的部分，随机采样一组点V。随后将会把V分成三类，$V_{open},V_{unvisited},V_{closed}$。首先在当前状态下搜索最近邻（可以是距离，也可以是最近k个点），并将这些点设置为$V_{open}$。随后在其中选在cost最小的一个点进行第步。同样对于这个点（如图中的$z$）进行最近邻搜索，找到最近的点如b中的$x$。随后在$x$上对$V_{open}$中的点进行连接，并进行障碍检测，删除掉无法连接的点。这样就可以得到一个局部最优的边。随后将无障碍物的边加入到树中，并将改点设为closed。进行下一次迭代，直到找到终点，或者所有点都访问玩。

![image-20210712142642866](Sample-based%20Motion%20Planning%20Algorithm.assets/image-20210712142642866.png)

简单的说FMT*就是一种具有树状结构的PRM。

![image-20210712144330211](Sample-based%20Motion%20Planning%20Algorithm.assets/image-20210712144330211.png)

甚至可以用广义的cost来代替举例等星系，实现一个类似距离场的算法。