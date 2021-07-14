# Optimal Boundary Value Problem for Polynomial of order 7

# 1 完全定义边界条件

> 在参考文献 A Computationally Efficient Motion Primitive for Quadrocopter Trajectory的基础上加入了minimum time的cost。

对于单段多项式轨迹，定义状态为$S = (p,v,a,j)$，那么状态转移方程为
$$
\dot S = (a,v,j,s)
$$
代价函数可以写成
$$
J = \int_0^T [\rho + s^2(t)]\mathrm dt
$$
其中$s$表示snap。根据庞特里亚金极小值原理，设$\vec\lambda = (\lambda_1, \lambda_2,\lambda_3,\lambda_4)$，并定义哈密顿函数$H(s,j,\lambda)$:
$$
\begin{aligned}
H(s,j,\lambda) &= s^2+\rho + \vec\lambda^Tf_s\\
&=s^2+\rho + \lambda_1v+\lambda_2a+\lambda_3j+\lambda_4s

\end{aligned}
$$
那么由极小值原理可以得到
$$
\dot{\vec\lambda} = -\nabla_SH (S^*,s^*,\vec\lambda)= (0,-\lambda_1,-\lambda_2,-\lambda_3)
$$
其中$s^*$表示最优的控制输入（snap）$S^*(t)$则表示最优状态轨迹。上述公式求解可以得到
$$
\vec\lambda(t) = 
\left[
\begin{matrix}
 -6\alpha\\
6\alpha t+2\beta\\
 -3\alpha t^2 - 2\beta t - 2\gamma\\
 \alpha t^3 + \beta^2 + 2\gamma t+ 2\omega
\end{matrix}
\right]
$$
那么零哈密顿函数对于输入的偏导数为0：$\frac{\partial H}{\partial s} = 2s+\lambda_4 = 0$最优输入轨迹可以求解得到
$$
\begin{aligned}
s^*(t)& = \arg\min_JH(S^*(t),s^*(t),\lambda(t))\\
&= \frac{-1}{2}\lambda_4\\
& = -\frac{1}{2}\alpha t^3 - \frac{1}{2}\beta t^2-\gamma t-\omega

\end{aligned}
$$
对最优状态轨迹进行多次积分，即可得到各阶的状态轨迹。
$$
S^*(t) = \left(\begin{array}{c} \frac{\alpha \,t^4}{8}+\frac{\beta \,t^3}{6}+\frac{\gamma \,t^2}{2}+\omega \,t+j_{0}\\ \frac{\alpha \,t^5}{40}+\frac{\beta \,t^4}{24}+\frac{\gamma \,t^3}{6}+\frac{\omega \,t^2}{2}+j_{0}\,t+a_{0}\\ \frac{\alpha \,t^6}{240}+\frac{\beta \,t^5}{120}+\frac{\gamma \,t^4}{24}+\frac{\omega \,t^3}{6}+\frac{j_{0}\,t^2}{2}+a_{0}\,t+v_{0}\\ \frac{\alpha \,t^7}{1680}+\frac{\beta \,t^6}{720}+\frac{\gamma \,t^5}{120}+\frac{\omega \,t^4}{24}+\frac{j_{0}\,t^3}{6}+\frac{a_{0}\,t^2}{2}+v_{0}\,t+p_{0} \end{array}\right)
$$


由此我们只需要求解$\alpha,\beta,\gamma,\omega$就可以得到最优状态轨迹。利用边界条件的对应关系。
$$
\left(\begin{array}{c} 
\alpha\\\beta\\\gamma\\\omega
\end{array}\right)
=
\left(\begin{array}{c} 33600\,p_{0}-33600\,j_{f}+16800\,a_{f}\,t+16800\,t\,v_{0}+3360\,a_{0}\,t^2+280\,j_{0}\,t^3+280\,p_{f}\,t^3-3360\,t^2\,v_{f}\\ 50400\,j_{f}\,t-50400\,p_{0}\,t-5400\,a_{0}\,t^3-24480\,a_{f}\,t^2-480\,j_{0}\,t^4-360\,p_{f}\,t^4-25920\,t^2\,v_{0}+4680\,t^3\,v_{f}\\ 1200\,a_{0}\,t^4+4680\,a_{f}\,t^3+120\,j_{0}\,t^5-10080\,j_{f}\,t^2+10080\,p_{0}\,t^2+60\,p_{f}\,t^5+5400\,t^3\,v_{0}-840\,t^4\,v_{f}\\ 840\,j_{f}\,t^3-360\,a_{f}\,t^4-16\,j_{0}\,t^6-120\,a_{0}\,t^5-840\,p_{0}\,t^3-4\,p_{f}\,t^6-480\,t^4\,v_{0}+60\,t^5\,v_{f} \end{array}\right)
$$

```matlab
syms p v a j s alpha beta gamma omega  rho real
由极小值原理，获得最优控制输入状态轨迹
syms t real
s_star = 0.5 * alpha * t^3 + 0.5 * beta * t^2+ gamma *t + omega
对控制输入进行多次积分，并在考虑并加入边界条件，即可得到各阶最优状态轨迹
syms p_0 v_0 a_0 j_0 p_f v_f a_f j_f real
j_star = int(s_star, t);    
a_star = int(j_star, t);    
v_star = int(a_star, t);    
p_star = int(v_star, t);    
j_star = j_star + j_0;
a_star = a_star + a_0 + j_0 * t;    
v_star = v_star + v_0 + a_0*t + j_0 * t^2/2;    
p_star = p_star + p_0+ v_0*t + a_0/2*t^2 + j_0/6*t^3;    
S_star = [
j_star;a_star;v_star;p_star;
]
由于初始条件以及隐含在积分中，我们将终止条件T代入得
syms T real
subs(S_star, t, T)
end_state = [p_f v_f, a_f, j_f]'
eqns=[];
for i = 1:4
    cur_eqn = S_star(i,:) == end_state(i);
    eqns = [eqns ;cur_eqn];
end
vars = [alpha;beta;gamma;omega];
anss = solve(eqns,vars);
aa = expand(anss.alpha * t^7)
bb = expand(anss.beta * t^7)
cc = expand(anss.gamma * t^7)
dd = expand(anss.omega * t^7)
ans = [aa;bb;cc;dd]
latex(ans)
最后我们还可以获得代价函数的表达式
J =  expand(s_star^2+rho);
subs(J,alpha,aa);
subs(J,beta,bb);
subs(J,gamma,cc);
subs(J,omega,dd);

我们可以通过求导获得cost最小的时间分配
eqn = expand(diff(J,t))

```

# 2 只固定起点全状态和终点位置

只固定位置和速度的最优控制问题也是我们采样算法中常常用到的。

前面部分的推导基本是一致的

对于单段多项式轨迹，定义状态为$S = (p,v,a,j)$，那么状态转移方程为
$$
\dot S = (a,v,j,s)
$$
代价函数可以写成
$$
J = \int_0^T [\rho + s^2(t)]\mathrm dt
$$
需要注意的是，我们的代价函数通常由两部分组成
$$
J  = cost_p + cost_f
$$
其中$cost_p$表示2过程代价，例如我们的代价函数就是过程代价。除此之外通常还有终端代价，在我们的costfunc中为0。

其中$s$表示snap。根据庞特里亚金极小值原理，设$\vec\lambda = (\lambda_1, \lambda_2,\lambda_3,\lambda_4)$，并定义哈密顿函数$H(s,j,\lambda)$:
$$
\begin{aligned}
H(s,j,\lambda) &= s^2+\rho + \vec\lambda^Tf_s\\
&=s^2+\rho + \lambda_1v+\lambda_2a+\lambda_3j+\lambda_4s

\end{aligned}
$$
那么由极小值原理可以得到变量满足
$$
\dot{\vec\lambda} = -\nabla_SH (S^*,s^*,\vec\lambda)= (0,-\lambda_1,-\lambda_2,-\lambda_3)
$$
对于自由变量，庞特里亚金极小值原理有一个新的约束，称为横截条件（Transversal Condition）【注意横截条件只针对终止时间T成立】，其中$h$表示终端代价，因此在我们的问题中均为0
$$
 \lambda_i(T) = \frac{\part h}{\part x_i}(T) = 0
$$
我们需要将哈密顿函数对状态求导
$$
\lambda_2(T) = \frac{\partial H}{\part v} = 0\\
\lambda_3(T) = \frac{\partial H}{\part a} = 0\\
\lambda_4(T) = \frac{\partial H}{\part j} = 0\\
$$
其中$s^*$表示最优的控制输入（snap）$S^*(t)$则表示最优状态轨迹。上述公式求解可以得到
$$
\vec\lambda(t) = 
\left[
\begin{matrix}
 -6\alpha\\
6\alpha t+2\beta\\
 -3\alpha t^2 - 2\beta t - 2\gamma\\
 \alpha t^3 + \beta^2 + 2\gamma t+ 2\omega
\end{matrix}
\right]
$$
那么零哈密顿函数对于输入的偏导数为0：$\frac{\partial H}{\partial s} = 2s+\lambda_4 = 0$最优输入轨迹可以求解得到
$$
\begin{aligned}
s^*(t)& = \arg\min_JH(S^*(t),s^*(t),\lambda(t))\\
&= \frac{-1}{2}\lambda_4\\
& = -\frac{1}{2}\alpha t^3 - \frac{1}{2}\beta t^2-\gamma t-\omega

\end{aligned}
$$
对最优状态轨迹进行多次积分，即可得到各阶的状态轨迹。
$$
S^*(t) = \left(\begin{array}{c} \frac{\alpha \,t^4}{8}+\frac{\beta \,t^3}{6}+\frac{\gamma \,t^2}{2}+\omega \,t+j_{0}\\ \frac{\alpha \,t^5}{40}+\frac{\beta \,t^4}{24}+\frac{\gamma \,t^3}{6}+\frac{\omega \,t^2}{2}+j_{0}\,t+a_{0}\\ \frac{\alpha \,t^6}{240}+\frac{\beta \,t^5}{120}+\frac{\gamma \,t^4}{24}+\frac{\omega \,t^3}{6}+\frac{j_{0}\,t^2}{2}+a_{0}\,t+v_{0}\\ \frac{\alpha \,t^7}{1680}+\frac{\beta \,t^6}{720}+\frac{\gamma \,t^5}{120}+\frac{\omega \,t^4}{24}+\frac{j_{0}\,t^3}{6}+\frac{a_{0}\,t^2}{2}+v_{0}\,t+p_{0} \end{array}\right)
$$

$$
\left(\begin{array}{c} 
\alpha\\\beta\\\gamma\\\omega
\end{array}\right)
=\left(\begin{array}{c}
\frac{8\,p_f }{t^4 }-\frac{8\,j_0 }{t^4 }\\
\frac{24\,j_0 }{t^3 }-\frac{24\,p_f }{t^3 }\\
\frac{12\,p_f }{t^2 }-\frac{12\,j_0 }{t^2 }\\
\frac{4\,j_0 }{t}-\frac{4\,p_f }{t}
\end{array}\right)
$$

# 3 固定位置和速度

只固定位置和速度的最优控制问题也是我们采样算法中常常用到的。

前面部分的推导基本是一致的

对于单段多项式轨迹，定义状态为$S = (p,v,a,j)$，那么状态转移方程为
$$
\dot S = (a,v,j,s)
$$
代价函数可以写成
$$
J = \int_0^T [\rho + s^2(t)]\mathrm dt
$$
需要注意的是，我们的代价函数通常由两部分组成
$$
J  = cost_p + cost_f
$$
其中$cost_p$表示2过程代价，例如我们的代价函数就是过程代价。除此之外通常还有终端代价，在我们的costfunc中为0。

其中$s$表示snap。根据庞特里亚金极小值原理，设$\vec\lambda = (\lambda_1, \lambda_2,\lambda_3,\lambda_4)$，并定义哈密顿函数$H(s,j,\lambda)$:
$$
\begin{aligned}
H(s,j,\lambda) &= s^2+\rho + \vec\lambda^Tf_s\\
&=s^2+\rho + \lambda_1v+\lambda_2a+\lambda_3j+\lambda_4s

\end{aligned}
$$
那么由极小值原理可以得到变量满足
$$
\dot{\vec\lambda} = -\nabla_SH (S^*,s^*,\vec\lambda)= (0,-\lambda_1,-\lambda_2,-\lambda_3)
$$
对于自由变量，庞特里亚金极小值原理有一个新的约束，称为横截条件（Transversal Condition）【注意横截条件只针对终止时间T成立】，其中$h$表示终端代价，因此在我们的问题中均为0
$$
 \lambda_i(T) = \frac{\part h}{\part x_i}(T) = 0
$$
我们需要将哈密顿函数对状态求导
$$
\lambda_2(T) = \frac{\partial H}{\part v} = 0\\
\lambda_3(T) = \frac{\partial H}{\part a} = 0\\
\lambda_4(T) = \frac{\partial H}{\part j} = 0\\
$$
其中$s^*$表示最优的控制输入（snap）$S^*(t)$则表示最优状态轨迹。上述公式求解可以得到
$$
\vec\lambda(t) = 
\left[
\begin{matrix}
 -6\alpha\\
6\alpha t+2\beta\\
 -3\alpha t^2 - 2\beta t - 2\gamma\\
 \alpha t^3 + \beta^2 + 2\gamma t+ 2\omega
\end{matrix}
\right]
$$
那么零哈密顿函数对于输入的偏导数为0：$\frac{\partial H}{\partial s} = 2s+\lambda_4 = 0$最优输入轨迹可以求解得到
$$
\begin{aligned}
s^*(t)& = \arg\min_JH(S^*(t),s^*(t),\lambda(t))\\
&= \frac{-1}{2}\lambda_4\\
& = -\frac{1}{2}\alpha t^3 - \frac{1}{2}\beta t^2-\gamma t-\omega

\end{aligned}
$$
对最优状态轨迹进行多次积分，即可得到各阶的状态轨迹。
$$
S^*(t) = \left(\begin{array}{c} \frac{\alpha \,t^4}{8}+\frac{\beta \,t^3}{6}+\frac{\gamma \,t^2}{2}+\omega \,t+j_{0}\\ \frac{\alpha \,t^5}{40}+\frac{\beta \,t^4}{24}+\frac{\gamma \,t^3}{6}+\frac{\omega \,t^2}{2}+j_{0}\,t+a_{0}\\ \frac{\alpha \,t^6}{240}+\frac{\beta \,t^5}{120}+\frac{\gamma \,t^4}{24}+\frac{\omega \,t^3}{6}+\frac{j_{0}\,t^2}{2}+a_{0}\,t+v_{0}\\ \frac{\alpha \,t^7}{1680}+\frac{\beta \,t^6}{720}+\frac{\gamma \,t^5}{120}+\frac{\omega \,t^4}{24}+\frac{j_{0}\,t^3}{6}+\frac{a_{0}\,t^2}{2}+v_{0}\,t+p_{0} \end{array}\right)
$$

$$
\left(\begin{array}{c} 
\alpha\\\beta\\\gamma\\\omega
\end{array}\right)
=
\left(\begin{array}{c} \frac{160\,v_{f}}{t^5}-\frac{40\,j_{0}}{t^4}-\frac{120\,p_{f}}{t^4}-\frac{160\,a_{0}}{t^5}\\ \frac{360\,a_{0}}{t^4}+\frac{96\,j_{0}}{t^3}+\frac{264\,p_{f}}{t^3}-\frac{360\,v_{f}}{t^4}\\ \frac{120\,v_{f}}{t^3}-\frac{36\,j_{0}}{t^2}-\frac{84\,p_{f}}{t^2}-\frac{120\,a_{0}}{t^3}\\ \frac{20\,a_{0}}{t^2}+\frac{8\,j_{0}}{t}+\frac{12\,p_{f}}{t}-\frac{20\,v_{f}}{t^2} \end{array}\right)
$$

# 4 只固定位置、速度和加速度

$$
\left(\begin{array}{c} 
\alpha\\\beta\\\gamma\\\omega
\end{array}\right)
=
\left(\begin{array}{c} \frac{360\,p_{f}}{t^4}-\frac{120\,j_{0}}{t^4}-\frac{960\,a_{0}}{t^5}-\frac{2400\,v_{0}}{t^6}-\frac{1440\,v_{f}}{t^5}+\frac{2400\,v_{f}}{t^6}\\ \frac{1800\,a_{0}}{t^4}+\frac{240\,j_{0}}{t^3}-\frac{600\,p_{f}}{t^3}+\frac{4320\,v_{0}}{t^5}+\frac{2520\,v_{f}}{t^4}-\frac{4320\,v_{f}}{t^5}\\ \frac{132\,p_{f}}{t^2}-\frac{72\,j_{0}}{t^2}-\frac{480\,a_{0}}{t^3}-\frac{1080\,v_{0}}{t^4}-\frac{600\,v_{f}}{t^3}+\frac{1080\,v_{f}}{t^4}\\ \frac{60\,a_{0}}{t^2}+\frac{12\,j_{0}}{t}-\frac{12\,p_{f}}{t}+\frac{120\,v_{0}}{t^3}+\frac{60\,v_{f}}{t^2}-\frac{120\,v_{f}}{t^3} \end{array}\right)
$$

