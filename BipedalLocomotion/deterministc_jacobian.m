%% Deterministic Jacobian
clear all
clc
close all
%%
%% 
walk = load("periodicWalking_0.4m_0.5s_noPD.mat");
[x0m,a,tp,tm] = decode(walk.optimizedStack);
x0in= x0m';
f=@one_step_bipedal_for_returnmap;
x0=x0in;
[max_eigval,J,V,D]=jcbn_linearized(f,x0);
%%
figure()
plot(eig(J),'*r')
grid on