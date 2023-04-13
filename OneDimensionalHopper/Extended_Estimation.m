%% Extended Kalman Filter Approach

clear all
close all
%% Define the system
% global leg_params
% leg_params.L_0 = 0.2;
% leg_params.g = 9.8;
% leg_params.m = 2;
% leg_params.k = 50;
% leg_params.d = 0.2;
% % As damping increases, eigenvalues
% leg_params.f = 20;
% %f=@fslip_map;
% %% Define Experiment parameters
% %% Define the process noise
% noise1_sweep=-0.3:1e-4:(0.3-1e-4);
% noise1_mean=0;
% noise1_var=0.01;
% discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
% discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);
% %% Define the measurement noise
% noise2_sweep=0;%-0.1:1e-3:(0.1-1e-3);
% noise2_mean=0;
% noise2_var=1e-7;%0.001;
% discrete_noise2_pdf=normpdf(noise2_sweep,noise2_mean,sqrt(noise2_var));
% discrete_noise2_pdf=discrete_noise2_pdf/sum(discrete_noise2_pdf);
%% Define states
% increment=0.01;
% states=0.2:increment:2;
% all_states=0:increment:5;
% height_sweep=states+increment/2;
% height_sweep(end)=[];

global leg_params
leg_params.L_0 = 0.2;
leg_params.g = 9.8;
leg_params.m = 2;
leg_params.k = 2000;
leg_params.d = 5;%5;
%leg_params.f = 20;
%f=@fslip_map;
leg_params.force_mag=20;
w_n = sqrt(leg_params.k/leg_params.m);
xi = leg_params.d/(2*sqrt(leg_params.k*leg_params.m));
w_d = sqrt(1-xi^2)*w_n;
leg_params.force_periode=pi/w_n;
f=@fslip_map_sine_force;
%% Define Experiment parameters
MC_exp=660000;
%% Define the process noise
noise1_sweep=-1.5:1e-4:(1.5-1e-4);
noise1_mean=0;
noise1_var=0.05;
discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);
%% Define the measurement noise
noise2_sweep=0;%-0.1:1e-3:(0.1-1e-3);
noise2_mean=0;
noise2_var=eps;%0.001;
discrete_noise2_pdf=normpdf(noise2_sweep,noise2_mean,sqrt(noise2_var));
discrete_noise2_pdf=discrete_noise2_pdf/sum(discrete_noise2_pdf);
%% Define states
increment=0.005;
states=0.4:increment:1.5;
all_states=0:increment:5;
height_sweep=states+increment/2;
height_sweep(end)=[]; 
sample=round(length(height_sweep)/2);
%%
T=zeros(length(states)-1,length(states)-1);
%%
tic

% Extended Kalman filter
%
% -------------------------------------------------------------------------
%
% State space model is
% X_k+1 = f_k(X_k) + V_k+1   -->  state update
% Y_k = h_k(X_k) + W_k       -->  measurement
%
% V_k+1 zero mean uncorrelated gaussian, cov(V_k) = Q_k
% W_k zero mean uncorrelated gaussian, cov(W_k) = R_k
% V_k & W_j are uncorrelated for every k,j
%
% -------------------------------------------------------------------------
%
% Inputs:
% f = f_k
%fstate=@sine_impact_map;

% Q = Q_k+1
Q=noise1_var;
% h = h_k
%hmeas=@ascent_map;
hmeas=@(x)(x);
% y = y_k % loop

% R = R_k
R=noise2_var;
% del_f = gradient of f_k % loop
% del_h = gradient of h_k % loop
% x_hat = current state prediction
% P_hat = current error covariance (predicted)
%
% -------------------------------------------------------------------------
idx=1;
ext_mean_list=[];
ext_var_list=[];
for z=height_sweep(2:end)
    idx=idx+1;
    
    %y_dot_td=-sqrt(2*leg_params.g*(z-leg_params.L_0));
    x_hat= z;%y_dot_td;
    P_hat=Q;
    % J=jcbn_finite_diff(f,x0)
    fstate=@(x)(fslip_map_sine_force_fzero(z,x,0));
    del_f=jcbn_finite_diff(fstate,0);
    del_h=1;%jcbn_finite_diff(hmeas,x_hat);
    %
    % Outputs:
    % x_next = next state prediction
    % P_next = next error covariance (predicted)
    % x_dgr = current state estimate
    % P_dgr = current estimated error covariance
    %
    % -------------------------------------------------------------------------
  
    
% % %     y_hat = hmeas(x_hat);
% % %     y_tilde = z - y_hat;
% % %     t = del_h;
% % %     tmp = P_hat*t;
% % %     M = inv(t'*tmp+R+eps);
% % %     K = tmp*M;
% % %     p = del_f;
% % %     x_dgr = x_hat + K*y_tilde;
% % %     x_next = fstate(x_dgr);
% % %     P_dgr = P_hat - tmp*K';
% % %     P_next = p* P_dgr* p' + Q;
    
    % -------------------------------------------------------------------------

    x_k_f=fslip_map_sine_force_fzero(z,0,0);
    P_k_f=del_f*Q*del_f;
    
    % -------------------------------------------------------------------------
    
    %
    %z1=x_next^2/(2*leg_params.g)+leg_params.L_0;
    z1=x_k_f;
    P2=P_k_f;
    ext_mean_list=[ext_mean_list z1];
    ext_var_list=[ext_var_list P2];
    % -------------------------------------------------------------------------
    %
    height_pdf=normpdf(height_sweep,z1,sqrt(P2));
    all_height_pdf=normpdf(all_states,z1,sqrt(P2));
    height_pdf=height_pdf/sum(all_height_pdf);
    
    %height_pdf(1)=height_pdf(1)+normcdf(height_sweep(1),z1,sqrt(P2))+normcdf(all_states(end),z1,sqrt(P2))-normcdf(height_sweep(end)+increment,z1,sqrt(P2));
     height_pdf(1)=height_pdf(1)+normcdf(height_sweep(1),z1,sqrt(P2))+normcdf(all_states(end),z1,sqrt(P2))-normcdf(height_sweep(end)+increment,z1,sqrt(P2));

    T(:,idx)=height_pdf'/sum(height_pdf);
end
%P12=X2*diag(Wc)*Z2';                        %transformed cross-covariance
%K=P12*inv(P2);
%x=x1+K*(z-z1);                              %state update
%P=P1-K*P12';                                %covariance update
%%
T(1)=1;
sum(T,1);
Tnorm=(T)./sum(T,1);
Tnorm(isnan(Tnorm))=0;
cols_with_all_zeros = find(all(Tnorm==0));
Tnorm(1,cols_with_all_zeros)=1;
S=Tnorm(1,2:end);
R=Tnorm(2:end,2:end);
T_ext=Tnorm;
disp('----------------------------------------')
disp('Linear Transform is done')
toc

Reig=sort(eig(R),'descend');
fprintf('Number of jumps: %d \n',idx*2)
disp('First 2 Eigenvalues of Tbar')
disp(Reig(1)),disp(Reig(2))
%%
savename='T_ext_'+string(noise1_var)+'.mat';
save(savename,'T_ext','ext_mean_list','ext_var_list')
%%
discrete_noise_p_pdf=normpdf(noise1_sweep,0,sqrt(noise1_var));
discrete_noise_p_pdf=discrete_noise_p_pdf/sum(discrete_noise_p_pdf);

figure(8)
subplot(3,1,1)
plot(noise1_sweep,discrete_noise_p_pdf')
title('Known Noise PDF')
subplot(3,1,2)
plot(height_sweep,T_ext(:,sample))
hold on
titlestr=strcat('Output PDF ($$x_{',string(sample),'}$$='+string(height_sweep(sample))+')');
title(titlestr,'Interpreter','latex')
xlim([min(states) max(states)])
xlabel('States','Interpreter','latex')
subplot(3,1,3)
[H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
s=surf(H1,H2,T_ext(2:end,2:end));
s.EdgeColor = 'none';
title('Transition matrix excluding the absorbing state')
view(0,90)
axis tight
suptitle('Numerical Linearization Method')

%% Comparison
disp('----------------------------------------')
% SSE=sum(sum((T_sys-T_ext).^2));
% Fnorm=sqrt(SSE);
% disp('Frobenius Norm of difference between Transition Matrices Byl and UF')
% disp(Fnorm)
% disp(Fnorm/Fnorm_byl*100)
% 
% figure(3)
% e=eig(T_montecarlo);
% plot(real(e),imag(e),'r*')
% hold on
% e=eig(T_sys);
% plot(real(e),imag(e),'b*')
% e=eig(T_ext);
% plot(real(e),imag(e),'k*')
% grid on
% legend('Monte Carlo','Byl','UT')
% title('Eigenvalue comparison')
%%
[H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
figure()
s=surf(H1,H2,T_ext(2:end,2:end));
s(1).EdgeColor = 'none';
s(1).FaceAlpha=1;
xlabel('$$h_k$$','Interpreter','latex')
ylabel('$$h_{k+1}$$','Interpreter','latex')
title('Numerical Linearization Method','Interpreter','latex')
view(0,90)
xlim([min(height_sweep(2:end)) max(height_sweep)])
ylim([min(height_sweep(2:end)) max(height_sweep)])
axis square