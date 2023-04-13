%% Unscented Transform (Augmented)
% Comparison of Katie Byl's exhaustive calculation, markovian chain calculation by
% monte carlo simulation and transition matrix by unscented transform
%clear all, close all, clc
%% Define the system
% global leg_params
% leg_params.L_0 = 0.2;
% leg_params.g = 9.8;
% leg_params.m = 2;
% leg_params.k = 2000;
% leg_params.d = 5;
% leg_params.force_mag=20;
% w_n = sqrt(leg_params.k/leg_params.m);
% xi = leg_params.d/(2*sqrt(leg_params.k*leg_params.m));
% w_d = sqrt(1-xi^2)*w_n;
% leg_params.force_periode=pi/w_n;
% f=@fslip_map_sine_force_fzero;
% load('Tsys_CLT.mat')
% load('Tmontecarlo.mat')
%% Define the process noise
% noise1_sweep=-1.5:1e-3:(1.5-1e-3);
% noise1_mean=0;
% noise1_var=0.1;
% discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
% discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);
% %% Define the measurement noise
% noise2_sweep=0; % -0.1:1e-3:(0.1-1e-3);
% noise2_mean=0;
% noise2_var=eps; % 0.001;
% discrete_noise2_pdf=normpdf(noise2_sweep,noise2_mean,sqrt(noise2_var));
% discrete_noise2_pdf=discrete_noise2_pdf/sum(discrete_noise2_pdf);
% %% Define states
% increment=0.005;
% states=0.4:increment:1.5;
% all_states=0:increment:5;
% height_sweep=states+increment/2;
% height_sweep(end)=[];
%force_period_controller()
%% Initialize Transition matrix
T=zeros(length(states)-1,length(states)-1);
%%
tic
% Unscented Kalman Filter for nonlinear dynamic systems
% x1: x[k-1], w1: w[k-1], v1: v[k-1]
%           x = f(x1, w1)
%           z = h(x, v1)
idx=1;
UTpts=[];
mean_list=[];
var_list=[];

for y_0=height_sweep(2:end)
    idx=idx+1;
    % x_aug=[x, w_k1, v_k1];
    P0=eps;
    Rw=noise1_var;
    Rv=noise2_var+eps;
    P=diag([P0,Rw,Rv]);
    
    x_aug=[y_0, 0, 0];
    
    %Sigma points around reference point
    %Inputs:
    %       x: reference point
    %       P: covariance
    %       c: coefficient
    %Output:
    %       X: Sigma points
    
    n= length(x_aug);
    
    W0=-0.5;
    L=numel(x_aug);
    alpha=2e-3;                                 %default, tunable
    ki=0;                                       %default, tunable
    beta=2;                                     %default, tunable
    lambda=0.5;%alpha^2*(L+ki)-L;                    %scaling factor
    c=n+lambda;                                 %scaling factor
    %Wm=[lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
    Wm=[0.5,0.25,0.25];%
    %Wm=[0.66 0.17 0.17]
    %Wm(1:3)/sum(Wm(1:3));
    Wc=Wm(1:3);
    %Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance
    Wc=Wc/sum(Wc);
    A = diag(sqrt(n/(1-W0))*chol(P)');
    x_a1=x_aug';
    X = [x_a1, x_a1+A ,x_a1-A];
    
    y=zeros(n,1);
    Y=zeros(n,L);
    for k=1:L
        Y(:,k)=[f(X(1,k),X(2,k),X(3,k)); 0 ;0];
        y=y+Wm(k)*Y(:,k);
    end
    %x=f(x1, w_k, v_k);
    Y;
    y;
    Y1=Y-y(:,ones(1,L));
    P=Y1*diag(Wc)*Y1';
    P1=P(1);%*sqrt(7); %%%%%%%%%%%%%%%%%%%%%%%%%
    z1=Y(1);%* 0.9969;
    sqrt(P(1));
    mean_list=[mean_list z1];
    var_list=[var_list P1];
%     if idx==110
%         sdfsdf=5;
%         Y;
%         y;
%     end
    
    height_pdf=normpdf(height_sweep,z1,sqrt(P1));
    all_height_pdf=normpdf(all_states,z1,sqrt(P1));
    height_pdf=height_pdf/sum(all_height_pdf);
    
    height_pdf(1)=height_pdf(1)+normcdf(height_sweep(1),z1,sqrt(P1))+normcdf(all_states(end),z1,sqrt(P1))-normcdf(height_sweep(end)+increment,z1,sqrt(P1));
    
    % pd = makedist('Weibul','A',z1,'B',5);
    %b=P1/z1;
    %a=z1/b;
    %pd = makedist('Normal','mu',z1,'sigma',sqrt(P1));
    %pd = makedist('ExtremeValue','mu',z1,'sigma',sqrt(P1));
    % pd = makedist('InverseGaussian','mu',z1,'lambda',);
    %    pd = makedist('gamma','A',a,'B',b);
   % height_pdf=pdf(pd,height_sweep);
    %height_pdf(1)=height_pdf(1)+cdf(pd,height_sweep(1))+cdf(pd,all_states(end))-cdf(pd,height_sweep(end)+increment);
    
    
    % height_pdf = chi2pdf(height_sweep,sqrt(P1));
    
    T(:,idx)=(height_pdf'/sum(height_pdf));
    UTpts=[UTpts;X(1,:)];
    
    % var(T_sys(:,end))=  2.6078e-04;
end
%%
T(1)=1;
sum(T,1);
Tnorm=(T)./sum(T,1);
Tnorm(isnan(Tnorm))=0;
cols_with_all_zeros = find(all(Tnorm==0));
Tnorm(1,cols_with_all_zeros)=1;
S=Tnorm(1,2:end);
R=Tnorm(2:end,2:end);
T_unsc=Tnorm;
disp('----------------------------------------')
disp('Unscented Transform is done')
toc

Reig=sort(eig(R),'descend');
fprintf('Number of jumps: %d \n',idx*3)
disp('First 2 Eigenvalues of Tbar')
disp(Reig(1)),disp(Reig(2))

discrete_noise_p_pdf=normpdf(noise1_sweep,0,sqrt(Rw));
discrete_noise_p_pdf=discrete_noise_p_pdf/sum(discrete_noise_p_pdf);
sample=round(length(height_sweep)/2);
%%
% figure(3)
% 
% subplot(3,1,1)
% plot(noise1_sweep,discrete_noise_p_pdf')
% title('Known Noise PDF')
% 
% subplot(3,1,2)
% plot(height_sweep,T_unsc(:,sample)),hold on
% UTpts(sample,:);
% xline(UTpts(sample,1))
% xline(UTpts(sample,2))
% xline(UTpts(sample,3))
% titlestr=strcat('Output PDF (x_{',string(sample),'}='+string(height_sweep(sample))+')');
% title(titlestr)
% xlim([min(states) max(states)])
% xlabel('States')
% 
% subplot(3,1,3)
% [H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
% s=surf(H1,H2,T_unsc(2:end,2:end));
% s.EdgeColor = 'none';
% title('Transition matrix excluding the absorbing state')
% view(0,90)
% axis tight
% suptitle('Unscented Transform Method')
% 
% figure()
% [H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
% s=surf(H1,H2,T_sys(2:end,2:end)-T_unsc(2:end,2:end));
% s.EdgeColor = 'none';
% title('Transition matrix excluding the absorbing state')
% view(0,90)
% axis tight
% suptitle('Unscented Transform Method')
