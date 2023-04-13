function [xp, deltaqsDot] = impactModel(x)
De = DeFunc(x);
E2 = E2Func(x);
Upsilon = zeros(2,5);
R = RFunc();

deltaF2 = -(E2/De*transpose(E2))\E2*[eye(5,5);Upsilon];
deltaqeDot = De\transpose(E2)*deltaF2 + [eye(5,5);Upsilon];
deltaqsDot = [R, zeros(5,2)] * deltaqeDot;
xp = [R*x(1:5); deltaqsDot*x(6:10)];
end