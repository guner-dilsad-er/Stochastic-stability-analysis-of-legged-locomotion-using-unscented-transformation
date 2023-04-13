function [x0m,a,tp,tm] = decode(stack)
M=MFunc();

x0m = stack(end-9:end);

a = zeros(4,M+1);
a(1,3:5) = stack(1:3);
a(2,3:5) = stack(4:6);
a(3,3:5) = stack(7:9);
a(4,3:5) = stack(10:12);

[x0p, ~] = impactModel(x0m);

qp = x0p(1:5);
qm = x0m(1:5);

qpDot = x0p(6:10);
qmDot = x0m(6:10);

H0 = H0Func();
c = cFunc();
H = [H0;c];

tp = c*qp;
tm = c*qm;

tpDot = c*qpDot;
tmDot = c*qmDot;

cond1 = H*qp;
a(:,1) = cond1(1:4);

cond2 = H*qm;
a(:,end) = cond2(1:4);

cond3 = H*(qpDot/tpDot);
cond3 = cond3(1:4);
cond3 = cond3*(tm-tp)/M;
a(:,2) = cond3 + a(:,1);

cond4 = H*(qmDot/tmDot);
cond4 = cond4(1:4);
cond4 = cond4*(tm-tp)/M;
a(:,end-1) = a(:,end)-cond4;

end