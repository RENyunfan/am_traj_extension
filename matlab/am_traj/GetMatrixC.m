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

