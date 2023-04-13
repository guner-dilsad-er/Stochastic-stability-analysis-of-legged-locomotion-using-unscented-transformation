function [xe]=one_step_bipedal_for_returnmap(x0in)

walk = load("periodicWalking_0.4m_0.5s_noPD.mat");

[x0m,a,tp,tm] = decode(walk.optimizedStack);

x0 = x0in;

x0p = (impactModel(x0'))';

options = odeset('RelTol',1e-8,'MaxStep',1e-2, 'Events', @(t,x)eventFunc(t,x));

[t,x,te,xe,ie] = ode45(@(t,x)xDotFunc(t,x,a,tp,tm), 0:0.020:10, x0p, options);



