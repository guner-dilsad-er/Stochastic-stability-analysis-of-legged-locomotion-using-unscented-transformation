function [u] = u96Func(x,a,tp,tm)

h = hFunc(x,a,tp,tm);
hDot = hDotFunc(x,a,tp,tm);

LgLfh = LgLfhFunc(x,a,tp,tm);
L2fh = L2fhFunc(x,a,tp,tm);
global cont_num
Kp=[];
Kd=[];
switch cont_num
    case 1       
        Kp = diag([60 90 90 50]);
        Kd = diag([10 20 20 10]);
    case 2
        Kp = diag([20 20 20 20])/4;
        Kd = diag([20 20 20 20])/4;
    case 3
        %P controller
        Kp = diag([20 20 20 20])*2;
        Kd = diag([20 20 20 20])*0;
    case 4
        
        Kp = diag([20 20 20 20])*2;
        Kd = diag([20 20 20 20])/20;
        
    case 5
        Kp = diag([10 89 83 50]);
        Kd = diag([5.4 21 21 9]);
        % %
end
%
v = Kp*h + Kd*hDot;
u = -LgLfh\(L2fh + v);

end