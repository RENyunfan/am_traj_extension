function Q = GetMatrixQ(T)

for i = 4:7
    for l = 4:7
       Q(i+1,l+1) = i*(i-1)*(i-2)*(i-3)*l*(l-1)*(l-2)*(l-3)/(i+l-7)*T^(i+l-7) ;
    end 
end

end

