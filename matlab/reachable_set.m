clc;clear;close all;
max_a = 3;
max_v = 0;
init_a = 0;
init_v = [1,0,0];

end_p=[0,0,0];
a0=[0,2,0];
end_ps=[end_p];
reso=1;
for int_t = 0.1:0.05:1
    for ax = -max_a:reso:max_a
        for ay = -max_a:reso:max_a
            az=0;
%             for az = -max_a:reso:max_a
                if(ax^2+ay^2+az^2>max_a^2)
                    continue
                end
                end_p(1) = init_v(1) * int_t + 1/2*(ax+a0(1))*int_t^2;
                end_p(2) = init_v(2) * int_t + 1/2*(ay+a0(2))*int_t^2;
                end_p(3) = init_v(3) * int_t + 1/2*(az+a0(3))*int_t^2;
                end_ps = [end_ps;end_p];
%             end
        end
    end
end
plot(end_ps(:,1),end_ps(:,2),'.')
hold on
axis('equal')
plot(0,0,'r*')

init_a = 0;
init_v = [2,0,0];
end_p=[0,0,0];
a0=[0,2,0];
end_ps=[end_p];
reso=0.5;
for int_t = 0.1:0.05:1
    for ax = -max_a:reso:max_a
        for ay = -max_a:reso:max_a
            az=0;
%             for az = -max_a:reso:max_a
                if(ax^2+ay^2+az^2>max_a^2)
                    continue
                end
                end_p(1) = init_v(1) * int_t + 1/2*(ax+a0(1))*int_t^2;
                end_p(2) = init_v(2) * int_t + 1/2*(ay+a0(2))*int_t^2;
                end_p(3) = init_v(3) * int_t + 1/2*(az+a0(3))*int_t^2;
                end_ps = [end_ps;end_p];
%             end
        end
    end
end
plot(end_ps(:,1),end_ps(:,2),'r.')


init_a = 0;
init_v = [3,0,0];
end_p=[0,0,0];
a0=[0,2,0];
end_ps=[end_p];
reso=0.5;
for int_t = 0.1:0.05:1
    for ax = -max_a:reso:max_a
        for ay = -max_a:reso:max_a
            az=0;
%             for az = -max_a:reso:max_a
                if(ax^2+ay^2+az^2>max_a^2)
                    continue
                end
                end_p(1) = init_v(1) * int_t + 1/2*(ax+a0(1))*int_t^2;
                end_p(2) = init_v(2) * int_t + 1/2*(ay+a0(2))*int_t^2;
                end_p(3) = init_v(3) * int_t + 1/2*(az+a0(3))*int_t^2;
                end_ps = [end_ps;end_p];
%             end
        end
    end
end
plot(end_ps(:,1),end_ps(:,2),'o')