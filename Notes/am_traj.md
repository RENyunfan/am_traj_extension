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
接下来要计算从边界条件到多项式系数的转换矩阵A。考虑到轨迹可以表示为
$$
p = 
$$




