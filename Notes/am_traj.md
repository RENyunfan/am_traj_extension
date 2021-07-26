# AM-Traj

仍然以7阶多项式为例，实现交替优化算法。首先进行问题描述：

给定起点和终点的全状态
$$
s_0 = [p_0,v_0,a_0,j_0]\\
s_f = [p_f,v_f,a_f,j_f]
$$
和一组路标点（waypoint）
$$
p_1,p_2,\dots,p_{N-2}
$$
求解最优的分段多项式轨迹。设总共有$N$个点，则有$N+6$个固定状态和总共$4N$个状态。

因此我们可以写出状态向量
$$
\mathbf D = \left[\begin{array}{c}
p_0\\v_0\\a_0\\j_0\\
p_1\\v_1\\a_1\\j_1\\

\vdots\\p_{N-2}\\
p_f\\v_f\\a_f\\j_f
 \end{array}\right]\in \R^{N\times3}
$$
通过设计选择矩阵$C$对状态进行重新排列
$$
\left[\begin{array}{c}
\mathbf d_1\\\mathbf d_2\\\vdots \\\mathbf d_M\\
 \end{array}\right] = C
 \left[\begin{array}{c}
\mathbf d_F\\\mathbf d_p
 \end{array}\right]
$$

```matlab
% clear ;clc
% syms p0 v0 a0 j0 pf vf af jf p1 v1 a1 j1 p2 v2 a2 j2 real
% d=[p0 v0 a0 j0 p1 v1 a1 j1  p2 v2 a2 j2 pf vf af jf]';
% C = GetMatrixCT(length(d)/4);
% dfdp = inv(C) *d
function C = GetMatrixC(N)
    C = zeros(N*4,N*4);
    C(1:4,1:4) = eye(4);
    
    for i = 1:N-2
       C(4*i+1,i+4) = 1;
    end
    
    for i = 1:4
       C(4*N-4+i,N+2+i) = 1;
    end
    
    for i = 1:(N-2)
       C(4*i+2, 3*i-2 + N+6) =  1;
       C(4*i+3, 3*i-1 + N+6) =  1;
       C(4*i+4, 3*i + N+6) =  1;
    end
end
```

接下来要计算从边界条件到多项式系数的转换矩阵A。考虑到轨迹可以表示为
$$
p = c_0+c_1t+c_2t^2+c_3t^3+...+c_7t^7\\ 
= \mathbf p^T\mathbf T(t)
$$

分别取起点状态和终点状态，单段轨迹可以表示为
$$
\left(\begin{array}{cccccccc}
p_0\\v_0\\a_0\\j_0\\
p_f\\v_f\\a_f\\j_f
\end{array}\right)= 
\left(\begin{array}{cccccccc}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & 0 & 0 & 0 & 2 & 0 & 0\\
0 & 0 & 0 & 0 & 6 & 0 & 0 & 0\\
T^7  & T^6  & T^5  & T^4  & T^3  & T^2  & T & 1\\
7\,T^6  & 6\,T^5  & 5\,T^4  & 4\,T^3  & 3\,T^2  & 2\,T & 1 & 0\\
42\,T^5  & 30\,T^4  & 20\,T^3  & 12\,T^2  & 6\,T & 2 & 0 & 0\\
210\,T^4  & 120\,T^3  & 60\,T^2  & 24\,T & 6 & 0 & 0 & 0
\end{array}\right)
\left(\begin{array}{cccccccc}
c_7\\c_6\\c_5\\c_4\\c_3\\c_2\\c_1\\c_0\\
\end{array}\right)\\
\mathbf d = M\mathbf p
$$
考虑到我们需要求的代价函数为
$$
J = \int_o^T p^{(4)}(t)^2\mathrm d t
$$
因此可以将代价函数写成矩阵的形式考虑到核函数
$$
\begin{aligned}
&f(t)=\sum_{i} p_{i} t^{i} \\
&\Rightarrow f^{(4)}(t)=\sum_{i \geq 4} i(i-1)(i-2)(i-3) t^{i-4} p_{i} \\
&\Rightarrow\left(f^{(4)}(t)\right)^{2}=\sum_{i \geq 4, l \geq 4} i(i-1)(i-2)(i-3) l(l-1)(l-2)(l-3) t^{i+l-8} p_{i} p_{l} \\
&\Rightarrow J(T)=\int_{T_{j-1}}^{T_{j}}\left(f^{4}(t)\right)^{2} d t=\sum_{i \geq 4, l \geq 4} \frac{i(i-1)(i-2)(i-3) j(l-1)(l-2)(l-3)}{i+l-7}\left(T_{j}^{i+l-7}-T_{j-1}^{i+l-7}\right) p_{i} p_{l} \\
&\Rightarrow J(T)=\int_{T_{j-1}}^{T_{j}}\left(f^{4}(t)\right)^{2} d t
\end{aligned}
$$
因此Smoothness的代价函数可以写成
$$
J = p^TQp
$$
其中Q可以表示为
$$
\left[\begin{array}{cc}
&\vdots \\
\cdots & \frac{i(i-1)(i-2)(i-3) l(l-1)(l-2)(l-3)}{i+l-7} T^{i+l-7} & \ldots
\\
&\vdots&
\end{array}\right]
$$
对于7阶多项式，Minsnap的Q可以表示为
$$
Q=\left(\begin{array}{cccccccc}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 576\,T & 1440\,T^2  & 2880\,T^3  & 5040\,T^4 \\
0 & 0 & 0 & 0 & 1440\,T^2  & 4800\,T^3  & 10800\,T^4  & 20160\,T^5 \\
0 & 0 & 0 & 0 & 2880\,T^3  & 10800\,T^4  & 25920\,T^5  & 50400\,T^6 \\
0 & 0 & 0 & 0 & 5040\,T^4  & 20160\,T^5  & 50400\,T^6  & 100800\,T^7 
\end{array}\right)
$$


因此代价函数可以表示为
$$
J = \mathbf d_{fp}C^TM^{-T}QM^{-1}C\mathbf d_{fp}+\rho \sum T_i
$$
得到选择矩阵之后，我们就可以把整段轨迹的代价函数写成
$$
J=\left[\begin{array}{l}
\mathbf{d}_{F} \\
\mathbf{d}_{P}
\end{array}\right]^{T} \underbrace{C A^{-T} Q A^{-1} C^{T}}_{R}\left[\begin{array}{l}
\mathbf{d}_{F} \\
\mathbf{d}_{P}
\end{array}\right]=\left[\begin{array}{l}
\mathbf{d}_{F} \\
\mathbf{d}_{P}
\end{array}\right]^{T}\left[\begin{array}{ll}
R_{F F} & R_{F P} \\
R_{P F} & R_{P P}
\end{array}\right]\left[\begin{array}{l}
\mathbf{d}_{F} \\
\mathbf{d}_{P}
\end{array}\right]
$$
解得
$$
\mathbf d_P^*  =R_{PP}^{-1}R_{FP}^T\mathbf d_F
$$
