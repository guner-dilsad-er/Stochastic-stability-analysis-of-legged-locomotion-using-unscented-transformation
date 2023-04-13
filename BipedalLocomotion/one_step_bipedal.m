function [xe]=one_step_bipedal(x0in,disturbance)

walk = load("periodicWalking_0.4m_0.5s_noPD.mat");
addpath('../code')
[x0m,a,tp,tm] = decode(walk.optimizedStack);
x0 = x0m';
x0(5) = x0in;
stepCounter = 1;
uAll = [];
hAll = [];
tAll = [];
xAll = [];
FnAll = [];
FtAll = [];
vAll = [];
tTotal = 0;
tInc = 0;
frameNo = 0;

for k = 1:stepCounter
    x0p = (impactModel(x0'))';
    x0p(5) = x0p(5)+disturbance;
    options = odeset('RelTol',1e-8,'MaxStep',1e-2, 'Events', @(t,x)eventFunc(t,x));
    [t,x,te,xe,ie] = ode45(@(t,x)xDotFunc(t,x,a,tp,tm), 0:0.020:10, x0p, options);
    x0 = xe;
    %disp(xe);
    tInc = tInc + te;
end
if isempty(xe)
    xe=[];
else
    xe=xe(5);
end
