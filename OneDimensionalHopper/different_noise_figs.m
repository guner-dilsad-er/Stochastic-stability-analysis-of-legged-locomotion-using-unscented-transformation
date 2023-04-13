close all
clear all
my_fig_defaults()
noise1_sweep=-1:1e-4:(1-1e-4);
noise1_mean=0;
increment=0.005;
states=0.4:increment:1.5;
all_states=0:increment:5;
height_sweep=states+increment/2;
height_sweep(end)=[];
%%
T_b_muhat=[];  
figure()
for noise1_var=0.02:0.01:0.07
    discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
    discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);
    
    savename='Tsys_CLT_'+string(noise1_var)+'.mat';
    load(savename)
    

    
    subplot(2,1,1)
    plot(noise1_sweep,discrete_noise1_pdf)
    xlabel('Disturbance')
    ylabel('PDF')
    hold on
    title('Known Noise PDF')
    
    subplot(2,1,2)
    sample=43;
    plot(height_sweep,T_sys(:,sample)),hold on
    titlestr=strcat('Output PDF ($$x_{',string(sample),'}$$='+string(height_sweep(sample))+')');
    title(titlestr,'Interpreter','latex')
    xlim([min(states) max(states)])
    xlabel('States')
    ylabel('Probability')
    % subplot(3,1,3)
    % [H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
    % s=surf(H1,H2,T_sys(2:end,2:end));
    % s.EdgeColor = 'none';
    % colormap(parula)
    % title('Transition matrix excluding the absorbing state')
    % view(0,90)
    % axis tight
    
    [Vb,D]=eig(T_sys);
    Md_b=Vb(2:end,2);
    mu_hat_mfpt=1/(1-max(real(D(2,2)))) ;
    T_b_muhat=[ T_b_muhat mu_hat_mfpt];
end
legend(string(0.02:0.01:0.07))
% subplot(3,1,3)
% 
% plot(0.02:0.01:0.07,T_b_muhat,'-o')
% hold on
% set(gca, 'YScale', 'log')


sgtitle('Exhaustive Method')


figure()
plot(0.02:0.01:0.07,T_b_muhat,'-o')
hold on
grid on
set(gca, 'YScale', 'log')
xlabel('Variance of Noise')
ylabel('MFPT of stochastic one leg hopper')
title('Exhaustive Method')
%%
T_mc_muhat=[];
figure()
for noise1_var=0.02:0.01:0.07
    
    discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
    discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);
    savename='Tmontecarlo_'+string(noise1_var)+'.mat';
    load(savename)
    
    
    subplot(2,1,1)
    h1=histogram(noise1_mem,'FaceAlpha',0.2,'EdgeAlpha',0.2);
    h1.Normalization = 'probability';
    hold on
    title('Measured Noise Histogram')
    %xlim([-0.2 0.2])
    subplot(2,1,2)
    plot(height_sweep,T_montecarlo(:,sample)),hold on
    titlestr=strcat('Output PDF ($$x_{',string(sample),'}$$='+string(height_sweep(sample))+')');
    title(titlestr,'Interpreter','latex')
    xlim([min(states) max(states)])
    xlabel('States')
    
    %ylim([0 0.3])
    % subplot(3,1,3)
    % [H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
    % s=surf(H1,H2,T_montecarlo(2:end,2:end));
    % s.EdgeColor = 'none';
    % title('Transition matrix excluding the absorbing state')
    % %zlim([0 0.1])
    % view(0,90)
    % axis tight
    
    [Vm,D]=eig(T_montecarlo);
    Md_m=Vm(2:end,2);
    mu_hat_mfpt=1/(1-max(real(D(2,2)))) ;
    T_mc_muhat=[ T_mc_muhat mu_hat_mfpt];
end
legend(string(0.02:0.01:0.07))
sgtitle('Monte Carlo Method')


figure()
plot(0.02:0.01:0.07,T_mc_muhat,'-o')
hold on
grid on
set(gca, 'YScale', 'log')
xlabel('Variance of Noise')
ylabel('MFPT of stochastic one leg hopper')
title('Monte Carlo Method')


%%
T_u_muhat=[];
T_u_hd_muhat=[];
T_u_hd2_muhat=[];

  figure()
for noise1_var=0.02:0.01:0.07
    discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
    discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);
    
    savename='T_unsc_'+string(noise1_var)+'.mat';
    load(savename)
    
  
    
    subplot(2,1,1)
    plot(noise1_sweep,discrete_noise1_pdf)
    xlabel('Disturbance')
    ylabel('PDF')
    hold on
    title('Known Noise PDF')
    
    subplot(2,1,2)
    plot(height_sweep,T_unsc(:,sample)),hold on
    titlestr=strcat('Output PDF ($$x_{',string(sample),'}$$='+string(height_sweep(sample))+')');
    title(titlestr,'Interpreter','latex')
    xlim([min(states) max(states)])
    xlabel('States')
    ylabel('Probability')
    xline(height_sweep(sample))
    % subplot(3,1,3)
    % [H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
    % s=surf(H1,H2,T_sys(2:end,2:end));
    % s.EdgeColor = 'none';
    % colormap(parula)
    % title('Transition matrix excluding the absorbing state')
    % view(0,90)
    % axis tight
    
    [Vu,D]=eig(T_unsc);
    Md_u=Vu(2:end,2);
    mu_hat_mfpt=1/(1-max(diag(D(2:end,2:end)))) ;
    T_u_muhat=[ T_u_muhat mu_hat_mfpt];
    savename='T_unsc_high_stiff3000_'+string(noise1_var)+'.mat';
    load(savename)
    [Vu,D]=eig(T_unsc);
    Md_u=Vu(2:end,2);
    mu_hat_mfpt=1/(1-max(diag(D(2:end,2:end)))) ;
    T_u_hd_muhat= [T_u_hd_muhat mu_hat_mfpt];
       
    savename='T_unsc_high_stiff5000_'+string(noise1_var)+'.mat';
    load(savename)
    [Vu,D]=eig(T_unsc);
    Md_u=Vu(2:end,2);
    mu_hat_mfpt=1/(1-max(diag(D(2:end,2:end)))) ;
    T_u_hd2_muhat= [T_u_hd2_muhat mu_hat_mfpt];
    
    subplot(2,1,1)
    plot(noise1_sweep,discrete_noise1_pdf)
    xlabel('Disturbance')
    ylabel('PDF')
    hold on
    title('Known Noise PDF')
    
    subplot(2,1,2)
    plot(height_sweep,T_unsc(:,sample)),hold on
    titlestr=strcat('Output PDF ($$x_{',string(sample),'}$$='+string(height_sweep(sample))+')');
    title(titlestr,'Interpreter','latex')
    xlim([min(states) max(states)])
    xlabel('States')
    ylabel('Probability')
    xline(height_sweep(sample))
end
legend(string(0.02:0.01:0.07))

T_unsc(2:end,2:end);

% subplot(3,1,3)
% 
% plot(0.02:0.01:0.07,T_u_muhat,'-o')
% hold on
% set(gca, 'YScale', 'log')
% grid on
sgtitle('Unscented Transform')

%%
figure()
plot(0.02:0.01:0.07,T_u_muhat,'-o')
hold on
plot(0.02:0.01:0.07,T_u_hd2_muhat,'-o')
plot(0.02:0.01:0.07,T_u_hd_muhat,'-o')
%plot(0.02:0.01:0.07,T_b_muhat,'-o')

set(gca, 'YScale', 'log')
grid on
%legend('Unscented Transformation','Systematic Experiments')
xlabel('Variance of Noise')
ylabel('MFPT of stochastic one leg hopper')
title('Mean First Passsage Time')
legend('k=2000','k=3000','k=5000')
%title('Comparison of MFPT Values')

