[toc]

# Polynomial Trajectory Closed-form Optimization

> Yunfan REN
>
> renyunfan@outlook.com

本文算法的为了适应G. Lu, W. Xu, F. Zhang, "Model Predictive Control for Trajectory Tracking on Differentiable Manifolds", [https://arxiv.org/abs/2106.15233](https://www.youtube.com/redirect?event=video_description&redir_token=QUFFLUhqbXRlYWZSQkxuaXVkTVJQM3ljdTJqQjkxU2RNd3xBQ3Jtc0ttYXpGMTBZUlhWNkJNdG1CbG5ENHFvazVVVmlzZ3N6RmYxWDlkRnNyV3BVLWJMaXdfcW15cFZfeVNVSlBKZEU3MTFWZzc5LW5sby1pRzJDTTRtQmg2VTBkbzV0Q24tZ2RkRjJhVTZ1M3UxeXVoODZDNA&q=https%3A%2F%2Farxiv.org%2Fabs%2F2106.15233)中的MPC控制器，完全采用7阶多项式。通过微分平坦特性将7阶多项式映射到最高角速度连续的轨迹。

# 1 BVP问题

BVP问题致力于解决固定边界条件的单段多项式轨迹最优求解。本节算法全部是现在文件`obvp_solver.hpp`中

## 1.0 极小值原理求解

首先我们考虑单段多项式轨迹，主要参考文献【1】。

对于单段多项式轨迹，定义状态为$S = (p,v,a,j)$，那么状态转移方程为
$$
\dot S = (a,v,j,s) \tag{1-1}
$$
代价函数可以写成
$$
J = \int_0^T [\rho + s^2(t)]\mathrm dt \tag{1-2}
$$
其中$s$表示snap。根据庞特里亚金极小值原理，设$\vec\lambda = (\lambda_1, \lambda_2,\lambda_3,\lambda_4)$，并定义哈密顿函数$H(s,j,\lambda)$:
$$
\begin{aligned}
H(s,j,\lambda) &= s^2+\rho + \vec\lambda^Tf_s\\
&=s^2+\rho + \lambda_1v+\lambda_2a+\lambda_3j+\lambda_4s

\end{aligned} \tag{1-3}
$$
那么由极小值原理可以得到
$$
\dot{\vec\lambda} = -\nabla_SH (S^*,s^*,\vec\lambda)= (0,-\lambda_1,-\lambda_2,-\lambda_3) \tag{1-4}
$$
其中$s^*$表示最优的控制输入（snap）$S^*(t)$则表示最优状态轨迹。上述公式求解可以得到
$$
\vec\lambda(t) = 
\left[
\begin{matrix}
 6\alpha\\
-6\alpha t-4\beta\\
 3\alpha t^2 + 4\beta t + 2\gamma\\
 -\alpha t^3 - 2\beta^2 - 2\gamma t- 2\omega
\end{matrix}
\right] \tag{1-5}
$$
那么零哈密顿函数对于输入的偏导数为0：$\frac{\partial H}{\partial s} = 2s+\lambda_4 = 0$最优输入轨迹可以求解得到
$$
\begin{aligned}
s^*(t)& = \arg\min_JH(S^*(t),s^*(t),\lambda(t))\\
&= \frac{1}{2}\lambda_4\\
& = \frac{1}{2}\alpha t^3 + \beta t^2+\gamma t+\omega

\end{aligned} \tag{1-6}
$$
对最优状态轨迹进行多次积分，即可得到各阶的状态轨迹。
$$
S^*(t) = \left(\begin{array}{c}
\frac{\alpha \,t^4 }{8}+\frac{\beta \,t^3 }{3}+\frac{\gamma \,t^2 }{2}+\omega \,t+j_0 \\
\frac{\alpha \,t^5 }{40}+\frac{\beta \,t^4 }{12}+\frac{\gamma \,t^3 }{6}+\frac{\omega \,t^2 }{2}+j_0 \,t+a_0 \\
\frac{\alpha \,t^6 }{240}+\frac{\beta \,t^5 }{60}+\frac{\gamma \,t^4 }{24}+\frac{\omega \,t^3 }{6}+\frac{j_0 \,t^2 }{2}+a_0 \,t+v_0 \\
\frac{\alpha \,t^7 }{1680}+\frac{\beta \,t^6 }{360}+\frac{\gamma \,t^5 }{120}+\frac{\omega \,t^4 }{24}+\frac{j_0 \,t^3 }{6}+\frac{a_0 \,t^2 }{2}+v_0 \,t+p_0 
\end{array}\right) \tag{1-7}
$$
因此对于minsnap+mintime的轨迹来说，最优的状态就取决于参数$\alpha,\beta,\gamma,\omega$。由于我们固定了一部分的边界条件，我们可以将以上四个参数表示为边界条件的函数。

```matlab
clear
syms p v a j s alpha beta gamma omega t rho real
syms p_0 v_0 a_0 j_0 p_f v_f a_f j_f real
% 1） 写出最优状态轨迹
s_star = 0.5 * alpha * t^3 +  beta * t^2+ gamma *t + omega
j_star = int(s_star, t);    
a_star = int(j_star, t);    
v_star = int(a_star, t);    
p_star = int(v_star, t);    
j_star = j_star + j_0;
a_star = a_star + a_0 + j_0 * t;    
v_star = v_star + v_0 + a_0*t + j_0 * t^2/2;    
p_star = p_star + p_0+ v_0*t + a_0/2*t^2 + j_0/6*t^3;    

S_star = [j_star;a_star;v_star;p_star;]
```

## 1.1 完全固定边界条件的解

### 1.1.1 求解最优状态轨迹

根据极小值原理，对于固定的边界条件，可以直接求解。对于自由的边界条件，则要引入终止代价带来的横截条件（Transversal Condition）。对于完全固定边界条件的情况，我们只需要将边界条件代入公式（1-7）
$$
\left(\begin{array}{c}
j_f \\
a_f \\
v_f \\
p_f 
\end{array}\right) = \left(\begin{array}{c}
\frac{\alpha \,t^4 }{8}+\frac{\beta \,t^3 }{3}+\frac{\gamma \,t^2 }{2}+\omega \,t+j_0 \\
\frac{\alpha \,t^5 }{40}+\frac{\beta \,t^4 }{12}+\frac{\gamma \,t^3 }{6}+\frac{\omega \,t^2 }{2}+j_0 \,t+a_0 \\
\frac{\alpha \,t^6 }{240}+\frac{\beta \,t^5 }{60}+\frac{\gamma \,t^4 }{24}+\frac{\omega \,t^3 }{6}+\frac{j_0 \,t^2 }{2}+a_0 \,t+v_0 \\
\frac{\alpha \,t^7 }{1680}+\frac{\beta \,t^6 }{360}+\frac{\gamma \,t^5 }{120}+\frac{\omega \,t^4 }{24}+\frac{j_0 \,t^3 }{6}+\frac{a_0 \,t^2 }{2}+v_0 \,t+p_0 
\end{array}\right)
\tag{1-8}
$$
联立方程组求解得到
$$
\left(\begin{array}{c}
\alpha \\
\beta \\
\gamma \\
\omega
\end{array}\right) =  
\left(\begin{array}{c}
\frac{3360\,a_0 }{t^5 }-\frac{3360\,a_f }{t^5 }+\frac{280\,j_0 }{t^4 }+\frac{280\,j_f }{t^4 }+\frac{33600\,p_0 }{t^7 }-\frac{33600\,p_f }{t^7 }+\frac{16800\,v_0 }{t^6 }+\frac{16800\,v_f }{t^6 }\\
\frac{2340\,a_f }{t^4 }-\frac{2700\,a_0 }{t^4 }-\frac{240\,j_0 }{t^3 }-\frac{180\,j_f }{t^3 }-\frac{25200\,p_0 }{t^6 }+\frac{25200\,p_f }{t^6 }-\frac{12960\,v_0 }{t^5 }-\frac{12240\,v_f }{t^5 }\\
\frac{1200\,a_0 }{t^3 }-\frac{840\,a_f }{t^3 }+\frac{120\,j_0 }{t^2 }+\frac{60\,j_f }{t^2 }+\frac{10080\,p_0 }{t^5 }-\frac{10080\,p_f }{t^5 }+\frac{5400\,v_0 }{t^4 }+\frac{4680\,v_f }{t^4 }\\
\frac{60\,a_f }{t^2 }-\frac{120\,a_0 }{t^2 }-\frac{16\,j_0 }{t}-\frac{4\,j_f }{t}-\frac{840\,p_0 }{t^4 }+\frac{840\,p_f }{t^4 }-\frac{480\,v_0 }{t^3 }-\frac{360\,v_f }{t^3 }
\end{array}\right)\tag{1-9}
$$
将公式（1-9）代入公式（1-7）即可得到最优多项式轨迹的多项式系数

```matlab
syms T real
subs(S_star, t, T)
end_state = [j_f a_f, v_f, p_f]'
eqns=[];
for i = 1:4
    cur_eqn = S_star(i,:) == end_state(i);
    eqns = [eqns ;cur_eqn];
end
vars = [alpha;beta;gamma;omega];
anss = solve(eqns,vars);
aa = expand(anss.alpha );
bb = expand(anss.beta );
cc = expand(anss.gamma );
dd = expand(anss.omega );
abyw =expand( [aa;bb;cc;dd])
```

### 1.1.2 求解最优轨迹cost和最优时间分配

我们还可以获得代价函数关于$\alpha,\beta,\gamma,\omega$的表达式
$$
J = \frac{\alpha^2 \,t^7 }{28}+\frac{\alpha \,\beta \,t^6 }{6}+\frac{\alpha \,\gamma \,t^5 }{5}+\frac{\alpha \,\omega \,t^4 }{4}+\frac{\beta^2 \,t^5 }{5}+\frac{\beta \,\gamma \,t^4 }{2}+\frac{2\,\beta \,\omega \,t^3 }{3}+\frac{\gamma^2 \,t^3 }{3}+\gamma \,\omega \,t^2 +\omega^2 \,t+\rho \,t \tag{1-10}
$$
代入公式（1-9）即可获得最优状态的代价。
$$
J = \left(\begin{array}{c}
\rho \\
16\,{j_0 }^2 +8\,j_0 \,j_f +16\,{j_f }^2 \\
240\,a_0 \,j_0 +120\,a_0 \,j_f -120\,a_f \,j_0 -240\,a_f \,j_f \\
1200\,{a_0 }^2 -1680\,a_0 \,a_f +1200\,{a_f }^2 +960\,j_0 \,v_0 +720\,j_0 \,v_f +720\,j_f \,v_0 +960\,j_f \,v_f \\
10800\,a_0 \,v_0 +9360\,a_0 \,v_f -9360\,a_f \,v_0 -10800\,a_f \,v_f +1680\,j_0 \,p_0 -1680\,j_0 \,p_f +1680\,j_f \,p_0 -1680\,j_f \,p_f \\
25920\,{v_0 }^2 +48960\,v_0 \,v_f +25920\,{v_f }^2 +20160\,a_0 \,p_0 -20160\,a_0 \,p_f -20160\,a_f \,p_0 +20160\,a_f \,p_f \\
100800\,p_0 \,v_0 +100800\,p_0 \,v_f -100800\,p_f \,v_0 -100800\,p_f \,v_f \\
100800\,{p_0 }^2 -201600\,p_0 \,p_f +100800\,{p_f }^2 
\end{array}\right)^T
\left(\begin{array}{cccccccc}
t^8  \\ t^6  \\ t^5  \\ t^4  \\ t^3  \\ t^2  \\ t \\ 1
\end{array}\right)
$$
除此之外我们可以对$J$求时间的倒数，以得到最优的时间分配
$$
\frac{dJ}{dt} = 
\left(\begin{array}{c}
\rho \\
-16\,{j_0 }^2 -8\,j_0 \,j_f -16\,{j_f }^2 \\
240\,a_f \,j_0 -240\,a_0 \,j_f -480\,a_0 \,j_0 +480\,a_f \,j_f \\
-3600\,{a_0 }^2 +5040\,a_0 \,a_f -3600\,{a_f }^2 -2880\,j_0 \,v_0 -2160\,j_0 \,v_f -2160\,j_f \,v_0 -2880\,j_f \,v_f \\
37440\,a_f \,v_0 -37440\,a_0 \,v_f -43200\,a_0 \,v_0 +43200\,a_f \,v_f -6720\,j_0 \,p_0 +6720\,j_0 \,p_f -6720\,j_f \,p_0 +6720\,j_f \,p_f \\
-129600\,{v_0 }^2 -244800\,v_0 \,v_f -129600\,{v_f }^2 -100800\,a_0 \,p_0 +100800\,a_0 \,p_f +100800\,a_f \,p_0 -100800\,a_f \,p_f \\
604800\,p_f \,v_0 -604800\,p_0 \,v_f -604800\,p_0 \,v_0 +604800\,p_f \,v_f \\
-705600\,{p_0 }^2 +1411200\,p_0 \,p_f -705600\,{p_f }^2 
\end{array}\right)\left(\begin{array}{c}
t^8 \\
t^6 \\
t^5 \\
t^4 \\
t^3 \\
t^2 \\
t\\
1
\end{array}\right)
$$
以上多项式求根，再对每个根的cost进行evaluation，就可以得到最优的时间分配了。

```matlab
J = expand( int((s_star^2+rho),t))
eqn = expand(subs(J,alpha,aa));
eqn = expand(subs(eqn,beta,bb));
eqn = expand(subs(eqn,gamma,cc));
eqn = expand(subs(eqn,omega,dd));
% 注意这里为了形式简洁，乘以了t^7，需要后续除掉
J = expand(eqn*t^7)
[cs,ts] = coeffs(J,t)

dj = expand(diff(J/t^7,t)*t^8)
[cs,ts] = coeffs(dj,t)
```

![image-20210726234615779](Polynomial%20Trajectory%20Closed-form%20Optimization.assets/image-20210726234615779.png)



## 1.2 只固定位置和速度的解

### 1.2.1 最优状态轨迹

对于只固定终点位置和速度的解，在极小值原理中相当于释放了终点状态中的a和j。因此需要引入横截条件。极小值原理中，横截条件表示将代价函数中的终止代价对自由变量求导需要满足为0。在我们考虑的问题中，代价函数（1-2）只包含了过程代价，故终止代价$h = 0$，因此各个偏导数也为0。因此对于自由变量，只需要满足
$$
\lambda_i(T) = \frac{\part h}{\part x_i}(T) = 0 \tag{1-11}
$$
那么对于自由的$a,j$，需要满足方程$\lambda_3 = \lambda_4 = 0$，根据公式（1-5）可以得到两个方程
$$
\left(\begin{array}{c}

3\,\alpha \,t^2 +4\,\beta \,t+2\,\gamma =0\\
-\alpha \,t^3 -2\beta \,t^2 -2\,\gamma \,t-2\,\omega =0
\end{array}\right)\tag{1-12}
$$
再更根据固定的终止条件$p_f,v_f$，可以联立方程组
$$
\left(\begin{array}{c}
\frac{\alpha \,t^6 }{240}+\frac{\beta \,t^5 }{60}+\frac{\gamma \,t^4 }{24}+\frac{\omega \,t^3 }{6}+\frac{j_0 \,t^2 }{2}+a_0 \,t+v_0 =v_f \\
\frac{\alpha \,t^7 }{1680}+\frac{\beta \,t^6 }{360}+\frac{\gamma \,t^5 }{120}+\frac{\omega \,t^4 }{24}+\frac{j_0 \,t^3 }{6}+\frac{a_0 \,t^2 }{2}+v_0 \,t+p_0 =p_f \\
3\,\alpha \,t^2 +4\,\beta \,t+2\,\gamma =0\\
-\alpha \,t^3 -2\,\beta \,t^2 -2\,\gamma \,t-2\,\omega =0
\end{array}\right)
$$
求解得到
$$
\left(\begin{array}{c}
\alpha \\
\beta \\
\gamma \\
\omega
\end{array}\right) =\left(\begin{array}{c}
\frac{672\,a_0 }{t^5 }+\frac{84\,j_0 }{t^4 }+\frac{3024\,p_0 }{t^7 }-\frac{3024\,p_f }{t^7 }+\frac{2184\,v_0 }{t^6 }+\frac{840\,v_f }{t^6 }\\
\frac{3276\,p_f }{t^6 }-\frac{96\,j_0 }{t^3 }-\frac{3276\,p_0 }{t^6 }-\frac{738\,a_0 }{t^4 }-\frac{2376\,v_0 }{t^5 }-\frac{900\,v_f }{t^5 }\\
\frac{468\,a_0 }{t^3 }+\frac{66\,j_0 }{t^2 }+\frac{2016\,p_0 }{t^5 }-\frac{2016\,p_f }{t^5 }+\frac{1476\,v_0 }{t^4 }+\frac{540\,v_f }{t^4 }\\
\frac{252\,p_f }{t^4 }-\frac{12\,j_0 }{t}-\frac{252\,p_0 }{t^4 }-\frac{66\,a_0 }{t^2 }-\frac{192\,v_0 }{t^3 }-\frac{60\,v_f }{t^3 }
\end{array}\right)
$$
在得到四个中间变量后，我们可以依照同样的方法求得最优代价和最优时间

### 1.2.2 最优代价和最优时间

同样采样(1.1.2)中的方法，我们可以计算得到代价函数的表达式
$$
J =\left(\begin{array}{c}
t^6 \\
t^5 \\
t^4 \\
t^3 \\
t^2 \\
t\\
1 
\end{array}\right)/ t^7
$$
同样，代价函数对时间的偏导数为
$$
\left(\begin{array}{c}
-12\,{j_0 }^2 \\
-264\,a_0 \,j_0 \\
-1404\,{a_0 }^2 -1152\,j_0 \,v_0 -360\,j_0 \,v_f \\
2016\,j_0 \,p_f -4320\,a_0 \,v_f -2016\,j_0 \,p_0 -11808\,a_0 \,v_0 \\
-23760\,{v_0 }^2 -18000\,v_0 \,v_f -3600\,{v_f }^2 -20160\,a_0 \,p_0 +20160\,a_0 \,p_f \\
78624\,p_f \,v_0 -30240\,p_0 \,v_f -78624\,p_0 \,v_0 +30240\,p_f \,v_f \\
-63504\,{p_0 }^2 +127008\,p_0 \,p_f -63504\,{p_f }^2 
\end{array}\right)
\left(\begin{array}{c}
t^6 \\
t^5 \\
t^4 \\
t^3 \\
t^2 \\
t\\
1
\end{array}\right)
$$


最终可以看到算法估计的代价和最有时间是完全正确的，并且与实际轨迹的离散积分完全一致。

![image-20210726232100572](Polynomial%20Trajectory%20Closed-form%20Optimization.assets/image-20210726232100572.png)

## 1.3 BVP C++代码

```cpp
        if ((type == FULLFIX)) {
            coeffsSnapObjective(0) = (8 * j_0 * j_f + 16 * j_0.square() + 16 * j_f.square()).sum();
            coeffsSnapObjective(1) = (240 * a_0 * j_0 + 120 * a_0 * j_f - 120 * a_f * j_0 - 240 * a_f * j_f).sum();
            coeffsSnapObjective(2) = (960 * j_0 * v_0 - 1680 * a_0 * a_f + 720 * j_0 * v_f + 720 * j_f * v_0 +
                                      960 * j_f * v_f + 1200 * a_0.square() + 1200 * a_f.square()).sum();
            coeffsSnapObjective(3) = (10800 * a_0 * v_0 + 9360 * a_0 * v_f - 9360 * a_f * v_0 - 10800 * a_f * v_f +
                                      1680 * j_0 * p_0 - 1680 * j_0 * p_f + 1680 * j_f * p_0 - 1680 * j_f * p_f).sum();
            coeffsSnapObjective(4) = (20160 * a_0 * p_0 - 20160 * a_0 * p_f - 20160 * a_f * p_0 + 20160 * a_f * p_f +
                                      48960 * v_0 * v_f + 25920 * v_0.square() + 25920 * v_f.square()).sum();
            coeffsSnapObjective(5) = (100800 * p_0 * v_0 + 100800 * p_0 * v_f - 100800 * p_f * v_0 -
                                      100800 * p_f * v_f).sum();
            coeffsSnapObjective(6) = (100800 * p_0.square() - 201600 * p_0 * p_f + 100800 * p_f.square()).sum();
        } else if ((type == FIXPV)) {

            coeffsSnapObjective(0) = (12 * j_0.square()).sum();
            coeffsSnapObjective(1) = (132 * a_0 * j_0).sum();
            coeffsSnapObjective(2) = (468 * a_0.square() + 384 * j_0 * v_0 + 120 * j_0 * v_f).sum();
            coeffsSnapObjective(3) = (2952 * a_0 * v_0 + 1080 * a_0 * v_f + 504 * j_0 * p_0 -
                                      504 * j_0 * p_f).sum();
            coeffsSnapObjective(4) = (4752 * v_0.square() + 3600 * v_0 * v_f + 720 * v_f.square() +
                                      4032 * a_0 * p_0 - 4032 * a_0 * p_f).sum();
            coeffsSnapObjective(5) = (13104 * p_0 * v_0 + 5040 * p_0 * v_f - 13104 * p_f * v_0 -
                                      5040 * p_f * v_f).sum();
            coeffsSnapObjective(6) = (9072 * p_0.square() - 18144 * p_0 * p_f + 9072 * p_f.square()).sum();
        }else {
            fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "Error, undefined bvp type!\n");
        }

```



## 1.4 Constrained BVP

使用二分法搜索得到最优的约束时间分配

![image-20210727000540062](Polynomial%20Trajectory%20Closed-form%20Optimization.assets/image-20210727000540062.png)

## 1.5 时间分析

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

采用保守时间分配，省略掉第二步，可以把速度提升到9.7us

# 2 AM分段多项式优化

# 参考文献

1. “A Computationally Efficient Motion Primitive for Quadrocopter Trajectory Generation | IEEE Journals & Magazine | IEEE Xplore.” n.d. Accessed March 29, 2021. https://ieeexplore.ieee.org/document/7299672?reload=true.
2. Wang, Z., X. Zhou, C. Xu, J. Chu, and F. Gao. 2020. “Alternating Minimization Based Trajectory Generation for Quadrotor Aggressive Flight.” *IEEE Robotics and Automation Letters* 5 (3): 4836–43. https://doi.org/10.1109/LRA.2020.3003871.
3. Richter, Charles, Adam Bry, and Nicholas Roy. 2016. “Polynomial Trajectory Planning for Aggressive Quadrotor Flight in Dense Indoor Environments.” In *Robotics Research*, edited by Masayuki Inaba and Peter Corke, 114:649–66. Springer Tracts in Advanced Robotics. Cham: Springer International Publishing. https://doi.org/10.1007/978-3-319-28872-7_37.

