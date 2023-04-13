%% Solving Approximate Analytical Map for Sine force input
% SLIP with Constant Forcing (F-SLIP) Template
% Error in initial condition and all solutions

function y_next=fslip_map_sine_force_fzero(y_0,disturbance,output_disturbance)
%% System Paramaters
%   clear all
%   clc, close all
% y_0=1;
% vel_disturbance=0;
% output_disturbance=0;
% leg_params.L_0 = 0.2;
% leg_params.g = 9.8;
% leg_params.m = 2;
% leg_params.k = 2000;
% leg_params.d = 10;
% leg_params.force_mag=50;
% leg_params.force_periode=0.05/(2*pi);

global leg_params
L_0=leg_params.L_0 ;
g=leg_params.g ;
m=leg_params.m ;
k=leg_params.k;%+500*disturbance ;
d=leg_params.d;%+output_disturbance ;
L=leg_params.force_mag;%+3*disturbance;
T=leg_params.force_periode;
omega=2*pi/T;
T=[];
%%
w_n = sqrt(k/m);
xi = d/(2*sqrt(k*m));
w_d = sqrt(1-xi^2)*w_n;
y_td = L_0;
%%
if y_0<y_td
    disp('Error in heigth')
    y_next=[];
   % y_next=y_0;
    return
else
y_dot_td=-sqrt(2*g*(y_0-y_td))+disturbance;

D=-m*g/k+y_td;

alpha=k-m*omega^2;
beta=d*omega;

A=beta*L/(alpha^2+beta^2);
B=-alpha*L/(alpha^2+beta^2);

phi=-atan2((-y_dot_td+B*omega-w_n*xi*(m*g/k-A))/w_d,m*g/k-A);
phi_2=atan2(w_d,xi*w_n);
C=sqrt(((-y_dot_td+B*omega-w_n*xi*(m*g/k-A))/w_d)^2+(m*g/k-A)^2);

%%
%y_sln=@(t)(C*exp(-xi*w_n*t)*cos(w_d*t-phi)+A*cos(t*omega)+B*sin(t*omega)+D);
y_sln=@(t)(C*exp(-xi*w_n*t)*cos(w_d*t-phi)+A*cos(t*omega)+B*sin(t*omega)+D);
y_dot_sln=@(t)(-C*exp(-xi*w_n*t)*(w_n*xi*cos(w_d*t-phi)+w_d*sin(w_d*t-phi))-A*omega*sin(t*omega)+B*omega*cos(t*omega));

%%
t_b1=(pi/2+phi+phi_2)/(w_d);
t_b=fzero(y_dot_sln,t_b1);
counter=0;
while t_b<0 && counter<10
    t_b=fzero(y_dot_sln,-t_b) ;
    counter=counter+1;
end
% t_b
%%
%t_lo1

cos_lo1=(m*g/k-A)/(C*k*exp(-2*xi*w_n*t_b));
 if abs(cos_lo1)>1
    t_lo1=10e5;
 else
     sin_lo1=sqrt(1-cos_lo1^2);
     t_lo1=(-atan2(sin_lo1,cos_lo1)+phi+2*pi)/(w_d);
 end

cond1=@(x)(y_sln(x)-y_td);
% disp('condition tlo1 by analytical')
% cond1(t_lo1)
% disp('condition tlo1 by optimization')
 t_lo1=fzero(cond1,t_lo1);
% cond1(t_lo1c)
%% t_lo2
a=k-d*xi*w_n;
b=-d*w_d;
R_2=sqrt(a^2+b^2);
phi_3=atan2(b,a);
cos_lo2=(-k*A+m*g+d*B*omega)/(R_2*C*exp(-2*xi*w_n*t_b));
if abs(cos_lo2)>1
    t_lo2=10e5;
else
    sin_lo2=sqrt(1-cos_lo2^2);
    t_lo2=(-atan2(sin_lo2,cos_lo2)+phi+phi_3+2*pi)/(w_d);
   % t_lo2=(atan2(sin_lo2,cos_lo2)+phi+phi_3+2*pi)/(w_d)
end  

cond2=@(x)(k*(y_sln(x)-y_td)-d*y_dot_sln(x)+L*sin(x*omega));
% disp('condition tlo2 by analytical')
% cond2(t_lo2)
 t_lo2=fzero(cond2,t_lo2);
% disp('condition tlo2 by optimization')
% cond2(t_lo2c)
%%

t_lo=min(t_lo1,t_lo2);

if t_lo>100
    disp('Error in liftoff time')
    y_next=[];
    %y_next=y_0;
    %y_next=[];
    %t_lo=2*t_b;
    return
end
leg_params.t_lo=t_lo;

%%
y_dot_lo=y_dot_sln(t_lo);%+output_disturbance;%+disturbance;
y_lo=y_sln(t_lo);

% if t_lo1<t_lo2
%     y_lo=y_td
% else
% end

leg_params.y_lo=y_lo;
y_reach=y_dot_lo^2/(2*g)+y_lo;
y_next=y_reach;
end