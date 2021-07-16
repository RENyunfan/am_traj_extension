function M = GetMatrixM(t_vec)
syms T real;
M_0 = [0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 2 0 0;
     0 0 0 0 6 0 0 0;
     T^7 T^6 T^5 T^4 T^3 T^2 T 1;
     7*T^6 6*T^5 5*T^4 4*T^3 3*T^2 2*T 1 0 ;
     42*T^5 30*T^4 20*T^3 12*T^2 6*T 2 0 0;
     210*T^4 120*T^3 60*T^2 24*T 6 0 0 0];
 num_t = length(t_vec)
 M = zeros(8*num_t, 8*num_t);
 for i = 1 : length(t_vec)
     M(8*(i-1)+1:8*(i) , 8*(i-1)+1:8*(i)) = subs(M_0,T,t_vec(i));
 end

end

