
%% cordic
function sina = sin_dds_cordic(angle,Iterate)
 
if angle < -pi/2 || angle > pi/2
    if angle < 0
        sina = sin_dds_cordic(angle+pi,Iterate);
    else
        sina = sin_dds_cordic(angle-pi,Iterate);
    end
   sina = -sina;
%     v = -v; % flip the sign for second or third quadrant
   return
end
% Iterate = 16;%迭代次数
x = zeros(Iterate+1,1);
y = zeros(Iterate+1,1);
z = zeros(Iterate+1,1);
x(1) = 0.607253;%初始设置

z(1) = angle;%待求角度θ
for i = 1:Iterate
    if z(i) >= 0
        d = 1;
    else
        d = -1;
    end
    x(i+1) = x(i) - d*y(i)*(2^(-(i-1)));
    y(i+1) = y(i) + d*x(i)*(2^(-(i-1)));
    z(i+1) = z(i) - d*atan(2^(-(i-1)));
end
sina = y(Iterate+1);

end
 
