% Comparison on Metastable Distributions of Matrices
% close all
% noise1_sweep=-1:1e-4:(1-1e-4);
% noise1_mean=0;
% increment=0.005;
% states=0.4:increment:1.5;
% all_states=0:increment:5;
% height_sweep=states+increment/2;
% height_sweep(end)=[];
function Matrix_Comparison(noise1_var,height_sweep)
%dilsad_thesis_fig_defaults()
savename='T_ext_'+string(noise1_var)+'.mat';
load(savename)
%T_ext,ext_mean_list,ext_var_list
savename='Tsys_CLT_'+string(noise1_var)+'_expmem.mat';
load(savename)
%T_sys
savename='Tmontecarlo_'+string(noise1_var)+'.mat';
load(savename)
%T_montecarlo,MC_cell,noise1_mem,MC_all,MC_stack,experiment_mem
savename='T_unsc_'+string(noise1_var)+'.mat';
load(savename)
%T_unsc,mean_list,var_list

sample=round(length(height_sweep)/2);
%T_montecarlo=sum(MC_cell{:})/length(MC_cell);
%% Comparison
disp('----------------------------------------')
Fnorm_byl=sqrt(sum(sum((T_sys.^2))));
disp('Frobenius norm of the Truth matrix')
disp(Fnorm_byl)
disp('----------------------------------------')
SSE=sum(sum((T_sys-T_montecarlo).^2));
Fnorm=sqrt(SSE);
disp('Frobenius Norm of difference between Transition Matrices Byl and Monte Carlo')
disp(Fnorm)
disp(Fnorm/Fnorm_byl*100)
disp('----------------------------------------')
SSE=sum(sum((T_sys-T_unsc).^2));
Fnorm=sqrt(SSE);
disp('Frobenius Norm of difference between Transition Matrices Byl and UF')
disp(Fnorm)
disp(Fnorm/Fnorm_byl*100)

figure()
e=eig(T_montecarlo);
plot(real(e),imag(e),'r*')
hold on
e=eig(T_sys);
plot(real(e),imag(e),'b*')
e=eig(T_unsc);
plot(real(e),imag(e),'k*')
e=eig(T_ext);
plot(real(e),imag(e),'g*')
grid on
legend('Monte Carlo','Byl','UT','Linearized')
title('Eigenvalue comparison','Interpreter','latex')




%%
figure()

subplot(4,1,1)
plot(height_sweep,T_unsc(:,2)),hold on
plot(height_sweep,T_sys(:,2),'r')
plot(height_sweep,T_montecarlo(:,2)),
plot(height_sweep,T_ext(:,2)),

grid on
%legend('UT','Byl','MC','Lin')
legend('Proposed','Systematic','MonteCarlo','Linearized','Location','southeast')

titlestr=strcat('Output PDF ($$x_{',string(2),'}$$='+string(height_sweep(2))+')');
title(titlestr,'Interpreter','latex')
axis tight

subplot(4,1,2)
plot(height_sweep,T_unsc(:,sample)),hold on
plot(height_sweep,T_sys(:,sample),'r')
plot(height_sweep,T_montecarlo(:,sample))
plot(height_sweep,T_ext(:,sample))
grid on
titlestr=strcat('Output PDF ($$x_{',string(sample),'}$$='+string(height_sweep(sample))+')');
title(titlestr,'Interpreter','latex')
axis tight

subplot(4,1,3)
plot(height_sweep,T_unsc(:,sample-20)),hold on
plot(height_sweep,T_sys(:,sample-20),'r')
plot(height_sweep,T_montecarlo(:,sample-20)),
plot(height_sweep,T_ext(:,sample-20)),
grid on
titlestr=strcat('Output PDF ($$x_{',string(sample-20),'}$$='+string(height_sweep(sample-20))+')');
title(titlestr,'Interpreter','latex')
axis tight
ylabel('Probability','Interpreter','latex')

subplot(4,1,4)
plot(height_sweep,T_unsc(:,end)),hold on
plot(height_sweep,T_sys(:,end),'r')
plot(height_sweep,T_montecarlo(:,end)),
plot(height_sweep,T_ext(:,end)),
grid on
titlestr=strcat('Output PDF ($$x_{',string(length(height_sweep)),'}$$='+string(height_sweep(end))+')');
title(titlestr,'Interpreter','latex')
axis tight
xlabel('States','Interpreter','latex')
sgtitle('Comparison of Output Distributions')

%%
sorted_exp=sortrows(experiment_mem);
mc_classify=[];
mc_means= [];

% figure()
% 
% subplot(2,1,1)
% for i=2:220
%     plot(height_sweep(i),COG(T_unsc(:,i))-COG(T_sys(:,i)),'ob'),hold on
%     plot(height_sweep(i),COG(T_sys(:,i))-COG(T_sys(:,i)),'or')
%     plot(height_sweep(i),COG(T_montecarlo(:,i))-COG(T_sys(:,i)),'ok'),grid on
% end
% subplot(2,1,2)
% 
% figure()
% for i=2:2:220
%     mc_classify=[];
% %     plot(height_sweep(i),mean_list(i-1),'.b'),hold on
% %     plot(height_sweep(i),mean_list(i)-sqrt(var_list(i-1)),'.b'),hold on
% %     plot(height_sweep(i),mean_list(i)+sqrt(var_list(i-1)),'.b'),hold on
%     errorbar(height_sweep(i),mean_list(i),sqrt(var_list(i-1)),'sb','MarkerSize',10),hold on
%     errorbar(height_sweep(i),ext_mean_list(i),sqrt(ext_var_list(i-1)),'sk','MarkerSize',10),hold on
% 
%     for j=1:length(sorted_exp(:,1))
%         if height_sweep(i)==sorted_exp(j,1)
%             mc_classify=[mc_classify;sorted_exp(j,4)];
%         end
%     end
%     mc_means=[mc_means mean(mc_classify)];
% %     plot(height_sweep(i),mean(mc_classify),'.r'),hold on
% %     plot(height_sweep(i),mean(mc_classify)+std(mc_classify),'.r'),hold on
% %     plot(height_sweep(i),mean(mc_classify)-std(mc_classify),'.r'),hold on    
%     errorbar(height_sweep(i),mean(mc_classify),std(mc_classify),'sr','MarkerSize',10)
% end
% %plot(height_sweep,mean(experiment_mem(:,4)),'ok'),grid on
% % legend('UT','Ex','MC')
% titlestr=strcat('Mean of output distributions');
% title(titlestr,'Interpreter','latex')
% axis tight
% axis square
%suptitle('Comparison of Output Distributions','Interpreter','latex')
%%
% figure()
% yyaxis left
% plot(height_sweep(2:end),mean_list)
% %sqrt(var_list),'sb','MarkerSize',10),
% hold on
% plot(height_sweep(2:end),ext_mean_list)%
% %,sqrt(ext_var_list),'sk','MarkerSize',10),hold on
% plot(height_sweep(2:end),mc_means)%,std(mc_classify),'sr','MarkerSize',10)
% yyaxis right
% plot(height_sweep(2:end),mean_list-mean_list)
% plot(height_sweep(2:end),ext_mean_list-mean_list)
% legend 

%%
sorted_exp=sortrows(experiment_mem);
mc_classify=[];
mc_means= [];
mc_std=[];
figure()
for i=2:220
    mc_classify=[];  
    for j=1:length(sorted_exp(:,1))
        if height_sweep(i)==sorted_exp(j,1)
            mc_classify=[mc_classify;sorted_exp(j,4)];
        end
    end
    mc_means=[mc_means mean(mc_classify)];
    mc_std=[mc_std std(mc_classify)];
end
sorted_exp_byl=sortrows(experiment_mem_byl);
byl_classify=[];
byl_classify_std=[];
byl_means= [];
byl_std=[];
for i=2:220
    byl_classify=[];  
    byl_classify_std=[];
    for j=1:length(sorted_exp_byl(:,1))
        if height_sweep(i)==sorted_exp_byl(j,1)
            byl_classify=[byl_classify;sorted_exp_byl(j,4)*sorted_exp_byl(j,5)];
            byl_classify_std=[byl_classify_std;sorted_exp_byl(j,4)*(1-sorted_exp_byl(j,5))];
        end
    end
    byl_means=[byl_means mean(byl_classify)];
    byl_std=[byl_std std(byl_classify_std)];
end
%%
figure()
plot(height_sweep(2:end),ext_mean_list,'LineWidth',2),hold on
plot(height_sweep(2:end),mean_list,'LineWidth',2),hold on
plot(height_sweep(2:end),mc_means,'LineWidth',2)
plot(height_sweep(2:end),byl_means*2e3,'LineWidth',2)
titlestr=strcat('Mean of output distributions');
title(titlestr,'Interpreter','latex')
axis tight
axis square
xlabel('States')
ylabel('Mean of output distribution')
legend('Linearized','Proposed','MonteCarlo','Systematic','Location','southeast')
figure()
plot(height_sweep(2:end),ext_mean_list-mean_list,'LineWidth',2),hold on
plot(height_sweep(2:end),abs(mc_means-mean_list),'LineWidth',2)
plot(height_sweep(2:end),byl_means*2e3-mean_list,'LineWidth',2)
legend('Linearized','MonteCarlo','Systematic')
titlestr=strcat('Mean Difference with Proposed Method');
title(titlestr,'Interpreter','latex')
axis tight
xlabel('States')
ylabel('Mean Difference of output distribution')
%%
figure()
plot(height_sweep(2:end),sqrt(ext_var_list),'LineWidth',2),hold on
plot(height_sweep(2:end),sqrt(var_list),'LineWidth',2)  
plot(height_sweep(2:end),mc_std,'LineWidth',2)
plot(height_sweep(2:end),byl_std*0.385783728798949,'LineWidth',2)
xlabel('States')
ylabel('Standard Deviation of output distribution')
titlestr=strcat('Standard Deviations of output distributions');
title(titlestr,'Interpreter','latex')
axis tight
axis square
legend('Linearized','Proposed','MonteCarlo','Systematic','Location','southeast')
figure()
plot(height_sweep(2:end),sqrt(ext_var_list)-sqrt(var_list),'LineWidth',2),hold on
plot(height_sweep(2:end),abs(mc_std-sqrt(var_list)),'LineWidth',2)
plot(height_sweep(2:end),abs(byl_std*0.385783728798949-sqrt(var_list)),'LineWidth',2)
xlabel('States')
ylabel('Standard Dev. Difference of output distribution')
titlestr=strcat('Standard Dev. Difference with Proposed Method');
title(titlestr,'Interpreter','latex')
axis tight
legend('Linearized','MonteCarlo','Systematic','Location','southeast')

disp(' Sum of Mean Errors:') 

sum(abs(ext_mean_list-mc_means))
sum(abs(mc_std-sqrt(ext_var_list)))
sum(abs(mean_list-mc_means))
sum(abs(mc_std-sqrt(var_list)))

% titlestr=strcat('Mean of output distributions');
% title(titlestr,'Interpreter','latex')
% axis tight
%suptitle('Comparison of Output Distributions','Interpreter','latex')

mean_error=mean_list-mc_means;
disp('Mean of mean errors')
mean(mean_error)
disp('Std of mean errors')
std(mean_error)
%%
KL_U=[];
KL_max=[];
figure()
for idx=1:length(height_sweep)
    kl_u=kldiv(height_sweep',T_unsc(:,idx),T_sys(:,idx));
    KL_U=[KL_U;kl_u];
    kl_m= kldiv(height_sweep',ones(length(T_sys(:,idx)),1)/length(T_sys(:,idx)),T_sys(:,idx));
    KL_max=[KL_max kl_m];
end
%errorbar(height_sweep(1:end),mean(MC_cell,1),std(MC_cell,1),'LineWidth',2)
x=height_sweep(1:end);
std_dev = std(MC_cell,1);
y=mean(MC_cell,1);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r','FaceAlpha',0.2,'LineStyle','none'),hold on

p1=plot(x, y, 'r-', 'LineWidth', 3);
p2=plot(height_sweep(1:end),KL_U,'b-','LineWidth',3);grid on
legend([p1,p2],'MC','UT','Location','southeast')
xlabel('States')
ylabel('Relative Entropy')
%ylim([-5e-3 5e-3])
axis tight
%p3=plot(height_sweep(1:end),KL_max,'-','LineWidth',3);grid on
title('Relative Entropies of Output Distributions')




% figure()
% yyaxis left
% p3=plot(height_sweep(1:end), KL_U'./KL_max,'b-','LineWidth',3);hold on
% yyaxis right
% p4=plot(height_sweep(1:end), y./KL_max,'r-','LineWidth',3);grid on
% title('Relative Tolerance / Max divergence ratio')
% legend([p3,p4],'UT','MC','Location','southeast')
% xlabel('States')

%%
figure()

[Vm,D]=eig(T_montecarlo);
plot(height_sweep(1:end),abs(Vm(:,2))/norm(Vm(:,2)),'r')
Md_m=Vm(2:end,2);
disp('----------------------------------------')
disp('Mean First Passage Time obtained from MC ')
mu_hat_mfpt=1/(1-max(real(D(2:end)))) ;
disp(mu_hat_mfpt)
hold on
[Vb,D]=eig(T_sys);
plot(height_sweep(1:end),abs(Vb(:,2))/norm(Vb(:,2)),'b')
Md_b=Vb(2:end,2);

disp('----------------------------------------')
disp('Mean First Passage Time obtained from the Truth matrix ')
mu_hat_mfpt=1/(1-max(real(D(2:end)))) ;
disp(mu_hat_mfpt)

[Vu,D]=eig(T_unsc);
plot(height_sweep(1:end),abs(Vu(:,2))/norm(Vu(:,2)),'k')
Md_u=Vu(2:end,2);

disp('----------------------------------------')
disp('Mean First Passage Time obtained from the UT ')
mu_hat_mfpt=1/(1-max(real(D(2:end)))) ;
disp(mu_hat_mfpt)


[Ve,D]=eig(T_ext);
plot(height_sweep(1:end),abs(Ve(:,2))/norm(Ve(:,2)),'k')
Md_e=Vu(2:end,2);

disp('----------------------------------------')
disp('Mean First Passage Time obtained from the Linearization ')
mu_hat_mfpt=1/(1-max(real(D(2:end)))) ;
disp(mu_hat_mfpt)


ylim([0 0.03])
grid on
legend('Monte Carlo','Byl','UT')
title('Eigenvector Comparison')
xlabel('States')
%%

[H1,H2]=meshgrid(height_sweep,height_sweep);

figure()
phi=abs(Vb(:,2))/norm(Vb(:,2));
phi(1)=0;
Metas=T_sys.*phi/sum(phi);
s=surfc(H1,H2,Metas);
s(1).EdgeColor = 'none';
view(135,45)
colormap(parula)
view(0,90)
s(2).EdgeColor = 'none';
s(2).ZLocation = 'zmax';
title('Metastable distribution by Exhaustive Method')
hold on
xlabel('$$h_k$$','Interpreter','latex')
ylabel('$$h_{k+1}$$','Interpreter','latex')
axis square
%%
[H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
figure()
phi=abs(Vb(:,2))/norm(Vb(:,2));
phi(1)=0;
s=surfc(H1,H2,T_sys(2:end,2:end));%.*phi/sum(phi));
s(1).EdgeColor = 'none';
s(1).FaceAlpha=1;
view(135,45)
colormap(parula)
view(0,90)
s(2).EdgeColor = 'none';
s(2).ZLocation = 'zmax';
title('Exhaustive Method')
hold on
plot3(0.614675,0.614675,1,'k*')
plot3(height_sweep,height_sweep,ones(1,length(height_sweep)),'--k')
axis square
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

increment=0.005;
height_sweep_tight=0.4:increment/100:1.5+increment/200;
Tdet=zeros(length(height_sweep_tight),1);
for i=1:length(height_sweep_tight)
    x_k=height_sweep_tight(i);
    h_next=f(x_k,0,0);
    Tdet(i)=h_next;
end
plot3(height_sweep_tight,Tdet,ones(1,length(height_sweep_tight)),'-k')

%set(gca,'xaxisLocation','top')
[M,c] =contour(H1,H2,Metas(2:end,2:end),0:1e-4:1e-3,'r');
c.LineWidth = 1;
c.ZLocation=1;
xlim([min(height_sweep(2:end)) max(height_sweep)])
ylim([min(height_sweep(2:end)) max(height_sweep)])
xlabel('$$h_k$$','Interpreter','latex')
ylabel('$$h_{k+1}$$','Interpreter','latex')
%%
[H1,H2]=meshgrid(height_sweep(2:end),height_sweep(2:end));
figure()
phi=abs(Vb(:,2))/norm(Vb(:,2));
phi(1)=0;
s=surf(H1,H2,T_sys(2:end,2:end));
s(1).EdgeColor = 'none';
s(1).FaceAlpha=1;
xlabel('$$h_k$$','Interpreter','latex')
ylabel('$$h_{k+1}$$','Interpreter','latex')
title('Exhaustive Method')
view(0,90)
xlim([min(height_sweep(2:end)) max(height_sweep)])
ylim([min(height_sweep(2:end)) max(height_sweep)])
axis square
%%
[H1,H2]=meshgrid(height_sweep,height_sweep);
figure()
phi=abs(Vm(:,2))/norm(Vm(:,2));
phi(1)=0;
s=surf(H1,H2,T_montecarlo);%.*phi/sum(phi));
s(1).EdgeColor = 'none';
%view(135,45)
view(0,90)
% s(2).EdgeColor = 'r';
% s(2).ZLocation = 'zmin';
title('Monte Carlo Method')
xlabel('$$h_k$$','Interpreter','latex')
ylabel('$$h_{k+1}$$','Interpreter','latex')
xlim([min(height_sweep) max(height_sweep)])
ylim([min(height_sweep) max(height_sweep)])
axis square

%%
figure()
phi=abs(Vu(:,2))/norm(Vu(:,2));
phi(1)=0;
s=surf(H1,H2,T_unsc(1:end,1:end));%.*phi/sum(phi));
s(1).EdgeColor = 'none';
%view(135,45)
view(0,90)
%s(2).EdgeColor = 'r';
%s(2).ZLocation = 'zmax';
title('Unscented Transform')
hold on
xlabel('$$h_k$$','Interpreter','latex')
ylabel('$$h_{k+1}$$','Interpreter','latex')
%plot3(0.614675,0.614675,max(max(T_unsc.*phi/sum(phi))),'bo')
%plot3(height_sweep,height_sweep,zeros(1,length(height_sweep))*max(max(T_unsc.*phi/sum(phi))),'--k')
xlim([min(height_sweep) max(height_sweep)])
ylim([min(height_sweep) max(height_sweep)])
axis tight
axis square

%%

plot3(height_sweep_tight,Tdet,zeros(1,length(height_sweep_tight))*max(max(T_unsc.*phi/sum(phi))),'-k')



%C=T_unsc.*phi/sum(phi);

%hold on
%contour3(C(5:10:150,:),'r')


%%
%% MFPT Vectors
figure()
plot(height_sweep(2:end),abs(inv(eye(219)-T_sys(2:end,2:end))*ones(219,1)),'b')
hold on
plot(height_sweep(2:end),abs(inv(eye(219)-T_montecarlo(2:end,2:end))*ones(219,1)),'r')
plot(height_sweep(2:end),abs(inv(eye(219)-T_unsc(2:end,2:end))*ones(219,1)),'k')
plot(height_sweep(2:end),abs(inv(eye(219)-T_ext(2:end,2:end))*ones(219,1)),'g')

legend('Byl','MC','UT','Lin')
title('State Dependent MFPT Vector')
%  set(gca, 'YScale', 'log')
phi=abs(Vb(:,2))/norm(Vb(:,2));
phi(1)=[];
M=abs(inv(eye(219)-T_sys(2:end,2:end))*ones(219,1))'*phi
phi=abs(Vu(:,2))/norm(Vu(:,2));
phi(1)=[];
M=abs(inv(eye(219)-T_unsc(2:end,2:end))*ones(219,1))'*phi
xlabel('States')
%%
disp('KL divergence of True Metastable distribution from uniform, (2nd eigenvector)')
disp(kldiv(height_sweep(2:end)',ones(length(Md_m),1)/length(Md_m),Md_b/sum(Md_b)))
disp('KL divergence of Metastable distribution of Monte Carlo, (2nd eigenvector)')
disp(kldiv(height_sweep(2:end)',Md_m/sum(Md_m),Md_b/sum(Md_b)))
disp('KL divergence of Metastable distribution of UT, (2nd eigenvector)')
disp(kldiv(height_sweep(2:end)',Md_u/sum(Md_u),Md_b/sum(Md_b)))
%%

