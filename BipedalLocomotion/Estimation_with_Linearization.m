%% Extended Kalman Filter Approach
addpath('code/')

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
% fstate=@sine_impact_map;

% Q = Q_k+1

    %Rw=eps;
Rw=noise_var;
Rwp=eps;
Q=diag([Rwp,Rwp,Rwp,Rwp,Rwp,Rw,Rw,Rw,Rw,Rw]);
% del_f = gradient of f_k % loop
% del_h = gradient of h_k % loop
% x_hat = current state prediction
% P_hat = current error covariance (predicted)
%
% -------------------------------------------------------------------------
%f=@one_step_bipedal_for_noise_experiments(x0in,distruabnce);
idx=1;
ext_mean_list=[];
ext_var_list=[];
for y_0=sweep(2:end)
    idx=idx+1;
    disp(idx)
    x0in(10)=y_0;
    fstate=@(x)(one_step_bipedal_for_noise_experiments(x0in,x));
    del_f=jcbn_finite_diff(fstate,zeros(1,10));
%     del_h=1;%jcbn_finite_diff(hmeas,x_hat);
    %
    % Outputs:
    % x_next = next state prediction
    % P_next = next error covariance (predicted)
    % x_dgr = current state estimate
    % P_dgr = current estimated error covariance
    %
    % -------------------------------------------------------------------------

    x_k_f=one_step_bipedal_for_noise_experiments(x0in,0);
    P_k_f=del_f*Q*del_f';
    
    % -------------------------------------------------------------------------
    
    %
    z1=x_k_f(10);
    P2=P_k_f(10,10);
    ext_mean_list=[ext_mean_list z1];
    ext_var_list=[ext_var_list P2];
    % -------------------------------------------------------------------------
    %
    height_pdf=normpdf(sweep,z1,sqrt(P2));
    all_height_pdf=normpdf(all_states,z1,sqrt(P2));
    height_pdf=height_pdf/sum(all_height_pdf);
    
     height_pdf(1)=height_pdf(1)+normcdf(sweep(1),z1,sqrt(P2))+normcdf(all_states(end),z1,sqrt(P2))-normcdf(sweep(end)+increment,z1,sqrt(P2));

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
savename='T_ext_'+string(noise_var)+'.mat';
save(savename,'T_ext','ext_mean_list','ext_var_list')
%%
discrete_noise_p_pdf=normpdf(noise_sweep,0,sqrt(noise_var));
discrete_noise_p_pdf=discrete_noise_p_pdf/sum(discrete_noise_p_pdf);
sample=50;
figure(8)
subplot(3,1,1)
plot(noise_sweep,discrete_noise_p_pdf')
title('Known Noise PDF')
subplot(3,1,2)
plot(sweep,T_ext(:,sample))
hold on
titlestr=strcat('Output PDF ($$x_{',string(sample),'}$$='+string(sweep(sample))+')');
title(titlestr,'Interpreter','latex')
xlim([min(states) max(states)])
xlabel('States','Interpreter','latex')
subplot(3,1,3)
[H1,H2]=meshgrid(sweep(2:end),sweep(2:end));
s=surf(H1,H2,T_ext(2:end,2:end));
s.EdgeColor = 'none';
title('Transition matrix excluding the absorbing state')
view(0,90)
axis tight
suptitle('Numerical Linearization Method')

%% Comparison
disp('----------------------------------------')
%%
[H1,H2]=meshgrid(sweep(2:end),sweep(2:end));
figure()
%T_ext=T_ext(1:100,1:100);
s=surf(H1,H2,T_ext(2:end,2:end));
s(1).EdgeColor = 'none';
s(1).FaceAlpha = 1;
xlabel('$$\dot{q}^5_k$$','Interpreter','latex')
ylabel('$$\dot{q}^5_{k+1}$$','Interpreter','latex')
title('Numerical Linearization Method')
view(0,90)
xlim([min(sweep(2:end)) max(sweep)])
ylim([min(sweep(2:end)) max(sweep)])
caxis([0 0.1])
colorbar()

axis square
