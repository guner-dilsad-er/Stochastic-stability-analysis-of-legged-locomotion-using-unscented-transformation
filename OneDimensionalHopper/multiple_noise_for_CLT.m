%% Leg
% Comparison of Systematic calculation, markovian chain calculation by
% monte carlo simulation and transition matrix by unscented transform
clear all, close all, clc

%% Define the system
global leg_params
leg_params.L_0 = 0.2;
leg_params.g = 9.8;
leg_params.m = 2;
leg_params.k = 2000;
leg_params.k = 3000;
leg_params.d = 5;
leg_params.force_mag=20;
w_n = sqrt(leg_params.k/leg_params.m);
xi = leg_params.d/(2*sqrt(leg_params.k*leg_params.m));
w_d = sqrt(1-xi^2)*w_n;
leg_params.force_periode=pi/w_n;
f=@fslip_map_sine_force_fzero;
%% Define Experiment parameters
MC_exp=6600;
%% Define the Touchdown Velocity noise
noise1_sweep=-1:1e-4:(1-1e-4);
noise1_mean=0;
noise1_var=0.05;
discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);

figure()
plot(noise1_sweep,discrete_noise1_pdf,'LineWidth',5)

hold on
noise1_sweep=-1:5e-2:(1-5e-2);
discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);
bar(noise1_sweep,discrete_noise1_pdf/500)
xlabel('Disturbance','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
title('Zero mean Gaussian with variance 0.05','Interpreter','latex')
noise1_sweep=-1:1e-3:(1-1e-3);
noise1_mean=0;
noise1_var=0.05;
discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);

%% Define the Lift-off Velocity noise
noise2_sweep=0;%-0.5:5e-3:(0.5-5e-3);
noise2_mean=0;
noise2_var= 0;%0.01;
discrete_noise2_pdf=normpdf(noise2_sweep,noise2_mean,sqrt(noise2_var));
discrete_noise2_pdf=discrete_noise2_pdf/sum(discrete_noise2_pdf);
%% Define states
increment=0.005;
states=0.4:increment:1.5;
all_states=0:increment:5;
height_sweep=states+increment/2;
height_sweep(end)=[];
%% Deterministic Return Map

height_sweep_tight=0.4:increment/100:1.5+increment/200;

Tdet=zeros(length(height_sweep_tight),1);
for i=1:length(height_sweep_tight)
            x_k=height_sweep_tight(i);
            h_next=f(x_k,0,0);
            Tdet(i)=h_next;
            %  end
end

figure(15)
plot(height_sweep_tight,Tdet,'LineWidth',2)
hold on
plot(height_sweep_tight,height_sweep_tight,'LineWidth',2)
title('Deterministic Return Map')
axis tight
grid on
plot(0.614675,0.614675,'bo')
legend('Return Map','Unity Slope','Fixed Point','Location','Southeast')
xlabel('$$h_k$$','Interpreter','latex')
ylabel('$$h_{k+1}$$','Interpreter','latex')
%%
%% Systematic Experiments
for noise1_var=0.02:0.01:0.07
discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);
%% Initialize Transition matrix

T=zeros(length(states)-1,length(states)-1);

%% Fill Transition matrix for each height and noise value
tic
experiment_mem_byl=[];
for i=2:length(height_sweep)
    for j=1:length(noise1_sweep)
       % for k=1:length(noise2_sweep)
            x_k=height_sweep(i);
            w_k=noise1_sweep(j);
            v_k=0;%noise2_sweep(k);
            h_next=f(x_k,w_k,v_k);
            out1=i;
            if h_next>max(states) || h_next<min(states) || isempty(h_next)
                out2=1;
            else
                out2=find(min(abs(height_sweep - h_next)) == abs(height_sweep - h_next),1);
            end
            T(out2,out1)=T(out2,out1)+1*discrete_noise1_pdf(j);%*discrete_noise2_pdf(k);
            experiment_mem_byl=[experiment_mem_byl; x_k,w_k,v_k, h_next,discrete_noise1_pdf(j)];

        %end         
     end   
   % toc
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
colormap(parula)
title('Transition matrix excluding the absorbing state')
view(0,90)
axis tight
suptitle('Exhaustive Method')
%%
savename='Tsys_CLT_'+string(noise1_var)+'_expmem.mat';
save(savename,'T_sys','experiment_mem_byl')

%%
%% Monte Carlo Transition Matrix;
MC_cell=zeros(1,length(height_sweep));
MC_all= zeros(length(states)-1,length(states)-1);
MC_stack= zeros(length(states)-1,length(states)-1);
experiment_mem=zeros(1,4);
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
        v_k= noise2_mean + randn(1)*sqrt(noise2_var);
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
        experiment_mem=[experiment_mem; x_k,w_k,v_k, h_next];
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
    MC_stack(:,:,mc_batch)=T_montecarlo;
end
%%
T_montecarlo= MC_all/mc_batch;

savename='Tmontecarlo_'+string(noise1_var)+'.mat';
save(savename,'T_montecarlo','MC_cell','noise1_mem','MC_all','MC_stack','experiment_mem')

%%
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
end
%%
close all
for noise1_var=0.02:0.01:0.07
Augmented_UT();
savename='T_unsc_high_stiff3000_'+string(noise1_var)+'.mat';
save(savename,'T_unsc','mean_list','var_list')
end
%%
different_noise_figs();
%%
close all
%%
Matrix_Comparison(noise1_var,height_sweep)