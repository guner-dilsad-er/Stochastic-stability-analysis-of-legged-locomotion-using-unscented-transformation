%% Leg
% Comparison of Katie Byl's exhaustive calculation, markovian chain calculation by
% monte carlo simulation and transition matrix by unscented transform
clear all, close all, clc
%% Define the system
global leg_params
leg_params.L_0 = 0.2;
leg_params.g = 9.8;
leg_params.m = 2;
leg_params.k = 2000;
leg_params.d = 5;
leg_params.force_mag=20;
w_n = sqrt(leg_params.k/leg_params.m);
xi = leg_params.d/(2*sqrt(leg_params.k*leg_params.m));
w_d = sqrt(1-xi^2)*w_n;
leg_params.force_periode=pi/w_n;
f=@fslip_map_sine_force_fzero;
%% Define Experiment parameters
MC_exp=6600;
%% Define the process noise
noise1_sweep=-1.5:1e-4:(1.5-1e-4);
noise1_mean=0;
noise1_var=0.1;
discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);
%% Define the measurement noise
noise2_sweep=0; % -0.1:1e-3:(0.1-1e-3);
noise2_mean=0;
noise2_var=eps; % 0.001;
discrete_noise2_pdf=normpdf(noise2_sweep,noise2_mean,sqrt(noise2_var));
discrete_noise2_pdf=discrete_noise2_pdf/sum(discrete_noise2_pdf);
%% Define states
increment=0.005;
states=0.4:increment:1.5;
all_states=0:increment:5;
height_sweep=states+increment/2;
height_sweep(end)=[];
%force_period_controller()
%% Initialize Transition matrix
T=zeros(length(states)-1,length(states)-1);
%% Fill Transition matrix for each height and noise value
tic
for i=2:length(height_sweep)
    for j=1:length(noise1_sweep)
        % for k=1:length(noise2_sweep)
        x_k=height_sweep(i);
        w_k=noise1_sweep(j);
        v_k=0;
        h_next=f(x_k,w_k,v_k);
        out1=i;
        if h_next>max(states) || h_next<min(states) || isempty(h_next)
            out2=1;
        else
            out2=find(min(abs(height_sweep - h_next)) == abs(height_sweep - h_next),1);
        end
        T(out2,out1)=T(out2,out1)+1*discrete_noise1_pdf(j);
        %  end
    end
end
disp('Byl computation is done')
toc
%%
fprintf('Number of jumps: %d \n',(length(noise1_sweep)*length(height_sweep)))
%%
T(1)=1;
Tnorm=(T)./sum(T,1);
S=Tnorm(1,2:end);
R=Tnorm(2:end,2:end);
T_sys=Tnorm;
Reig=sort(eig(R),'descend');
disp('First 2 Eigenvalues of Tbar')
disp(Reig(1)),disp(Reig(2))
disp('--------------------------------')
%%
figure(1)

subplot(3,1,1)
plot(noise1_sweep,discrete_noise1_pdf)
title('Known Noise PDF')

subplot(3,1,2)
sample=round(length(height_sweep)/2);
plot(height_sweep,T_sys(:,sample)),hold on
titlestr=strcat('Output PDF (x_{',string(sample),'}='+string(height_sweep(sample))+')');
title(titlestr)
xlim([min(states) max(states)])
xlabel('States')

subplot(3,1,3)
[H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
s=surf(H1,H2,T_sys(2:end,2:end));
s.EdgeColor = 'none';
title('Transition matrix excluding the absorbing state')
view(0,90)
axis tight
suptitle('Exhaustive Method')
%%
save('Tsys.mat','T_sys')
%% Monte Carlo Transition Matrix;
MC_cell=zeros(1,length(height_sweep));
MC_all= zeros(length(states)-1,length(states)-1);

for mc_batch=1:50
    %% Initialize Transition matrix
    T=zeros(length(states)-1,length(states)-1);
    %% Fill Transition matrix for random height and noise value
    tic
    noise1_mem=zeros(1,MC_exp);
    noise2_mem=zeros(1,MC_exp);
    for exp=1:MC_exp
        i=randi(length(height_sweep)-1)+1;
        w_k= noise1_mean + randn(1)*sqrt(noise1_var);
        v_k= 0;%noise2_mean + randn(1)*sqrt(noise2_var);
        x_k=height_sweep(i);
        h_next=f(x_k,w_k,v_k);
        out1=i;
        if h_next>max(states) || h_next<min(states) || isempty(h_next)
            out2=1;
        else
            out2=find(min(abs(height_sweep - h_next)) == abs(height_sweep - h_next),1);
        end
        T(out2,out1)=T(out2,out1)+1;
        noise1_mem(exp)=w_k;
        noise2_mem(exp)=v_k;
    end
    disp('Monte Carlo computation is done')
    toc
    %%
    fprintf('Number of jumps: %d \n',exp)
    %%
    T(1)=1;
    Tnorm=(T)./sum(T,1);
    S=Tnorm(1,2:end);
    R=Tnorm(2:end,2:end);
    T_montecarlo=Tnorm;
    Reig=sort(eig(R),'descend');
    disp('First 2 Eigenvalues of Tbar')
    disp(Reig(1)),disp(Reig(2))
    KL_X=[];
    for idx=1:length(height_sweep)
        kl_x=kldiv(height_sweep',T_montecarlo(:,idx),T_sys(:,idx));
        KL_X=[KL_X;kl_x];
    end
    MC_cell(mc_batch,:)=KL_X';
    MC_all=MC_all+T_montecarlo;
end
%%
T_montecarlo= MC_all/mc_batch;
save('Tmontecarlo.mat','T_montecarlo','MC_cell')
figure(2)
subplot(3,1,1)
histogram(noise1_mem);
title('Measured Noise Histogram')
%xlim([-0.2 0.2])
subplot(3,1,2)
plot(height_sweep,T_montecarlo(:,sample)),hold on
titlestr=strcat('Output PDF (x_{',string(sample),'}='+string(height_sweep(sample))+')');
title(titlestr)
xlim([min(states) max(states)])
xlabel('States')
%ylim([0 0.3])
subplot(3,1,3)
[H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
s=surf(H1,H2,T_montecarlo(2:end,2:end));
s.EdgeColor = 'none';
title('Transition matrix excluding the absorbing state')
%zlim([0 0.1])
view(0,90)
axis tight
suptitle('Monte Carlo Method')
%% Unscented Transformation
close all
Augmented_UT();

Matrix_Comparison(T_sys, T_montecarlo,T_unsc,height_sweep,MC_cell); %MC_cell

%%
var(T_sys(2:end,110))

var(T_unsc(2:end,110))

% %%
% 
% T=zeros(length(states)-1,length(states)-1);
% 
% tic
% % Unscented Kalman Filter for nonlinear dynamic systems
% % [x, P] = ukf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P
% % for nonlinear dynamic system (for simplicity, noises are assumed as additive):
% %           x_k+1 = f(x_k) + w_k
% %           z_k   = h(x_k) + v_k
% % where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
% %       v ~ N(0,R) meaning v is gaussian noise with covariance R
% % Inputs:   f: function handle for f(x)
% fstate=@impact_map;
% %fstate=@fslip_map_sine_force_fzero;
% n=1;
% P=noise1_var;
% %           P: known process noise covariance
% %           R: measurement noise covariance
% R=eps;%noise2_var;
% %           h: function handle for h(x)
% hmeas=@ascent_map;
% idx=1;
% UTpts=[];
% for z=height_sweep(2:end)
%     idx=idx+1;
%     %           x: "a priori" state estimate
%     y_dot_td=-sqrt(2*leg_params.g*(z-leg_params.L_0));
%     
%     y_td= leg_params.L_0;
%     
%     x=z;
%     % Output:   x: "a posteriori" state estimate
%     %           P: "a posteriori" state covariance
%     % Reference: Julier, SJ. and Uhlmann, J.K., Unscented Filtering and
%     % Nonlinear Estimation, Proceedings of the IEEE, Vol. 92, No. 3,
%     % pp.401-422, 2004.
%     %
%     % By Yi Cao at Cranfield University, 04/01/2008
%     % x_n=fstate(x)
%     % hmeas(x_n)
%     
%     L=numel(x);                                 %number of states
%     m=numel(z);                                 %number of measurements
%     alpha=1e-3;                                 %default, tunable
%     ki=0;                                       %default, tunable
%     beta=2;                                     %default, tunable
%     lambda=alpha^2*(L+ki)-L;                    %scaling factor
%     c=n+lambda;                                 %scaling factor
%     Wm=[lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
%     Wc=Wm;
%     Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance
%     c=sqrt(c);
%     
%     % sigma points around x
%     X=sigmas(x,P,c);
%     
%     % unscented transformation of impact map
%     [x1,X1,P1,X2]=ut(fstate,X,Wm,Wc,L,P);
%     
%     % x1 : transformed mean
%     % X1 : transformed sampling points f(X)
%     % P1 : transformed covariance
%     % X2 : transformed deviations X1-x1
%     
%     % Unscented transformation of measurements
%     % sigma points for the second UT are taken as
%     % the transformed sampling points X1
%     
%     %[z1,Z1,P2,Z2]=ut(hmeas,X1,Wm,Wc,m,R);
%     
%     % z1: transformed mean
%     % Z1: transformed sampling points
%     % P2: transformed covariance
%     % Z2: transformed deviations
%     
%     % P12=X2*diag(Wc)*Z2';                        %transformed cross-covariance
%     % K=P12*inv(P2);
%     % x=x1+K*(z-z1);                              %state update
%     % P=P1-K*P12';                                %covariance update
%     
%     %     y_lo=z1(1);
%     % y_dot_lo=-z1(2);
%     %
%     % y_reach=y_dot_lo^2/(2*g)+y_lo;
%     % y_next=y_reach;
%     %
%     % z1=y_next;
%     z1=X1(1);P2=P1;
%     height_pdf=normpdf(height_sweep,z1,sqrt(P2));
%     all_height_pdf=normpdf(all_states,z1,sqrt(P2));
%     height_pdf=height_pdf/sum(all_height_pdf);
%     
%     height_pdf(1)=height_pdf(1)+normcdf(height_sweep(1),z1,sqrt(P2))+normcdf(all_states(end),z1,sqrt(P2))-normcdf(height_sweep(end)+increment,z1,sqrt(P2));
%     
%     T(:,idx)=(height_pdf');
%     UTpts=[UTpts;X1];
% end
% 
% %%
% T(1)=1;
% sum(T,1);
% Tnorm=(T)./sum(T,1);
% Tnorm(isnan(Tnorm))=0;
% cols_with_all_zeros = find(all(Tnorm==0));
% Tnorm(1,cols_with_all_zeros)=1;
% S=Tnorm(1,2:end);
% R=Tnorm(2:end,2:end);
% T_unsc=Tnorm;
% disp('----------------------------------------')
% disp('Unscented Transform is done')
% toc
% 
% Reig=sort(eig(R),'descend');
% fprintf('Number of jumps: %d \n',idx*3)
% disp('First 2 Eigenvalues of Tbar')
% disp(Reig(1)),disp(Reig(2))
% %%
% discrete_noise_p_pdf=normpdf(noise1_sweep,0,sqrt(P));
% discrete_noise_p_pdf=discrete_noise_p_pdf/sum(discrete_noise_p_pdf);
% sample=round(length(height_sweep)/2);
% 
% figure(5)
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
% function [y,Y,P,Y1]=ut(f,X,Wm,Wc,n,R)
% %Unscented Transformation
% %Input:
% %        f: nonlinear map
% %        X: sigma points
% %       Wm: weights for mean
% %       Wc: weights for covraiance
% %        n: numer of outputs of f
% %        R: additive covariance
% %Output:
% %        y: transformed mean
% %        Y: transformed sampling points
% %        P: transformed covariance
% %       Y1: transformed deviations
% 
% L=size(X,2);
% y=zeros(n,1);
% Y=zeros(n,L);
% for k=1:L
%     Y(:,k)=f(X(:,k),0,0);
%     y=y+Wm(k)*Y(:,k);
% end
% Y1=Y-y(:,ones(1,L));
% P=Y1*diag(Wc)*Y1'+R;
% end
% 
% function X=sigmas(x,P,c)
% %Sigma points around reference point
% %Inputs:
% %       x: reference point
% %       P: covariance
% %       c: coefficient
% %Output:
% %       X: Sigma points
% 
% A = c*chol(P)';
% Y = x(:,ones(1,numel(x)));
% X = [x Y+A Y-A]; end