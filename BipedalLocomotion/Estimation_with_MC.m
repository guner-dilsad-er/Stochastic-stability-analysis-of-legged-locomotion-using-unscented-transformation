% 10 D noise experiments %
% Body angular velocity Markov Chain %
clear all
clc
%%
increment=0.004;
states=(-0.28:increment:0.18)-0.848597863783639;
all_states=-10:increment:10;
sweep=states+increment/2;
sweep(end)=[];
%%
noise_sweep=-0.3:5e-5:(0.3-5e-5);
noise_mean=0;
noise_var=1e-3;
discrete_noise_pdf=normpdf(noise_sweep,noise_mean,sqrt(noise_var));
discrete_noise_pdf=discrete_noise_pdf/sum(discrete_noise_pdf);
figure()
plot(noise_sweep,discrete_noise_pdf)

%%
addpath 'code/'
walk = load("periodicWalking_0.4m_0.5s_noPD.mat");
[x0m,a,tp,tm] = decode(walk.optimizedStack);
x0in= x0m';

f=@one_step_bipedal_for_noise_experiments;
exp_num=100000;%10000;
%%
% Choose an initial state and noise variance:
noise_var=1e-3;
noise_mean=0;

%% Deterministic Return Map

sweep_tight=(-0.28:increment:0.18)-0.848597863783639;

Tdet=zeros(length(sweep_tight),1);
for i=1:length(sweep_tight)
            x_k=sweep_tight(i);
            x0in(10)=x_k;
            h_next=f(x0in,0);
            Tdet(i)=h_next(10);
            %  end
end
%%
figure(15)
plot(sweep_tight,Tdet','LineWidth',2)
hold on
plot(sweep_tight,sweep_tight,'LineWidth',2)
title('Deterministic Return Map')
axis tight
grid on
plot(-0.848597863783639,-0.848597863783639,'bo')
legend('Return Map','Unity Slope','Fixed Point','Location','Southeast')
xlabel('$$\dot{q}^5_k$$','Interpreter','latex')
ylabel('$$\dot{q}^5_{k+1}$$','Interpreter','latex')


%% Find output distributions experimentally
% Monte Carlo Experiment
noiseMem=[];
NoisyStateTransition=zeros(exp_num,20);
for j=1:exp_num
    x0in(10)=rand(1)*(abs(states(1)-states(end)))+states(1);
    disturbance=[normrnd(noise_mean, sqrt(noise_var),[1,5])*0,normrnd(noise_mean, sqrt(noise_var),[1,5])];%normrnd(noise_mean, sqrt(noise_var),[1,10]);
    noiseMem(j,:)= disturbance;
    [xe]=one_step_bipedal_for_noise_experiments(x0in,disturbance);
    if isempty(xe)
      %  continue
             NoisyStateTransition(j,:)=[x0in,x0in*999];
    else
        % Store information in State Transition Map
        xe=xe(1,:);
        %idx=(k-1)*exp_num+j;
        NoisyStateTransition(j,:)=[x0in,xe];
    end
end
disp('100k exp done')
%%
MC_cell=[];
%for z=0:9
T=zeros(length(states)-1,length(states)-1);
for j=1:exp_num
     h=NoisyStateTransition(j,10);
    h_next=NoisyStateTransition(j,20);
    if h>max(states) || h<min(states) || isempty(h)
        out1=1;
        disp('Unexpected Error')
        disp(h)
    else
        out1=find(min(abs(sweep - h)) == abs(sweep - h),1);
    end
    
    if h_next>max(states) || h_next<min(states) || isempty(h_next)
        out2=1;
    else
        out2=find(min(abs(sweep - h_next)) == abs(sweep - h_next),1);
    end
    T(out2,out1)=T(out2,out1)+1;%*discrete_noise2_pdf(k);
end


T(1)=1;
Tmc=T./sum(T,1);

% MC_cell{z+1}=Tmc;
%end
%%
save('MC_100k_thesis_narrow.mat','NoisyStateTransition','T','Tmc','MC_cell')
disp('MC done')
%%

global noise_var cont_num

    for noise_var=1e-4*(8:1:11)
 for cont_num=1:5       
        tic
        Estimation_with_UT()
        disp(cont_num)
        disp(noise_var)
        toc
    end
end
disp('UT done')
noise_var=1e-3;
cont_num=1;
%%
tic

Estimation_with_Linearization()
disp('Lin done')
toc
%%
figure()
subplot(3,1,1)
plot(sweep,Tmc(:,35))
hold on
plot(sweep,Tunsc(:,35))
plot(sweep,T_ext(:,35))
xlim([min(sweep(2:end)) max(sweep(2:end))])
titlestr=strcat('Output PDF ($$x_{',string(35),'}$$='+string(sweep(35))+')');
title(titlestr,'Interpreter','latex')
subplot(4,1,2)
plot(sweep,Tmc(:,40))
hold on
plot(sweep,Tunsc(:,40))
plot(sweep,T_ext(:,40))
titlestr=strcat('Output PDF ($$x_{',string(40),'}$$='+string(sweep(40))+')');
title(titlestr,'Interpreter','latex')
ylabel('Probability')
xlim([min(sweep(2:end)) max(sweep(2:end))])
subplot(4,1,3)
plot(sweep,Tmc(:,75))
hold on
plot(sweep,Tunsc(:,75))
plot(sweep,T_ext(:,75))
titlestr=strcat('Output PDF ($$x_{',string(75),'}$$='+string(sweep(75))+')');
title(titlestr,'Interpreter','latex')
axis tight
subplot(4,1,4)
plot(sweep,Tmc(:,100))
hold on
plot(sweep,Tunsc(:,100))
plot(sweep,T_ext(:,100))
titlestr=strcat('Output PDF ($$x_{',string(100),'}$$='+string(sweep(100))+')');
title(titlestr,'Interpreter','latex')
xlabel('States')
xlim([min(sweep(2:end)) max(sweep(2:end))])
%ylim([0 0.2])
legend('MC','UT','Lin','Location','west')
suptitle('Comparison of output distributions of different states')

saveas(gcf,'Comparison_thesis_narrow.fig')
%%
figure()
subplot(2,2,1)
[X,Y]=meshgrid(sweep(2:end),sweep(2:end));
surf(X,Y,T_ext(2:end,2:end)-Tunsc(2:end,2:end),'EdgeColor','none')
view(0,90)
subplot(2,2,2)
surf(X,Y,Tmc(2:end,2:end)-Tunsc(2:end,2:end),'EdgeColor','none')
view(0,90)
subplot(2,2,3)
surf(X,Y,Tmc(2:end,2:end)-T_ext(2:end,2:end),'EdgeColor','none')
view(0,90)
subplot(2,2,4)
surf(X,Y,T_ext(2:end,2:end)-Tunsc(2:end,2:end),'EdgeColor','none')
hold on
plot3(sweep,Tdet(2:end),ones(1,length(sweep)),'-k','LineWidth',2)
xlim([min(sweep(2:end)) max(sweep(2:end))])
ylim([min(sweep(2:end)) max(sweep(2:end))])
view(0,90)
%%
figure()
[X,Y]=meshgrid(sweep(2:end),sweep(2:end));
surf(X,Y,Tmc(2:end,2:end),'EdgeColor','none')
xlim([min(sweep(2:end)) max(sweep(2:end))])
ylim([min(sweep(2:end)) max(sweep(2:end))])

view(0,90)
titlestr='Monte Carlo Method with '+string(exp_num)+' experiments';
title(titlestr,'interpreter','latex')
xlabel('$$\dot{q}^5_k$$','Interpreter','latex')
ylabel('$$\dot{q}^5_{k+1}$$','Interpreter','latex')
axis square

colorbar()
saveas(gcf,'MC_state_transition_thesis_narrow.fig')
%caxis([0 0.08])
axis square

%%
figure()
[X,Y]=meshgrid(sweep(2:end),sweep(2:end));
surf(X,Y,Tmc(2:end,2:end),'EdgeColor','none')
hold on
plot3(sweep,Tdet(2:end),ones(1,length(sweep)),'-k','LineWidth',2)
xlim([min(sweep(2:end)) max(sweep(2:end))])
ylim([min(sweep(2:end)) max(sweep(2:end))])

view(0,90)
titlestr='Monte Carlo Method with '+string(exp_num)+' experiments';
title(titlestr,'interpreter','latex')
xlabel('$$\dot{q}^5_k$$','Interpreter','latex')
ylabel('$$\dot{q}^5_{k+1}$$','Interpreter','latex')
axis square
colorbar()

%%
sorted_exp=sortrows(NoisyStateTransition,10);
%%
mc_classify=[];
mc_means= [];
mc_vars=[];
%figure()
for i=1:1:114
    mc_classify=[];
    for j=1:length(sorted_exp(:,1))
        if sweep(i)>(sorted_exp(j,10)-increment/2) &  sweep(i)<(sorted_exp(j,10)+increment/2)
            mc_classify=[mc_classify;sorted_exp(j,20)];
        end
    end
    mc_means=[mc_means mean(mc_classify)];
    mc_vars=[mc_vars (std(mc_classify))];
end
%%
figure()
subplot(2,1,1)
plot(sweep(2:end),unsc_mean_list,'LineWidth',2)
hold on
plot(sweep(2:end),ext_mean_list,'LineWidth',2)
plot(sweep(2:end),mc_means,'LineWidth',2)
xlim([min(sweep(2:end)) max(sweep(2:end))])
ylabel('Mean Value')
subplot(2,1,2)
plot(sweep(2:end),sqrt(unsc_var_list),'LineWidth',2)
hold on
plot(sweep(2:end),sqrt(ext_var_list),'LineWidth',2)
plot(sweep(2:end),mc_vars,'LineWidth',2)
xlim([min(sweep(2:end)) max(sweep(2:end))])
xlabel('States')
ylabel('Standard Deviation Value')

%legend('UT','Lin','MC')
legend('Proposed','Linearized','MontaCarlo','Location','southeast')
%% 
Tnorm=Tunsc;
figure(),
subtightplot(2,1,1,0.05)

plot3(sweep(2:end),Tdet(3:end),0*ones(1,length(sweep(2:end))),'-k','LineWidth',3)

hold on
surf(X,Y,Tnorm(2:end,2:end),'EdgeColor','none','FaceAlpha',0.7);%,view(0,90);
lins=round(linspace(2,length(sweep(2:end)),20));
mesh(X(lins,lins),Y(lins,lins),Tnorm(lins,lins),'FaceColor','none')
grid on
axis tight
ylabel('$$x_{k+1}$$','Interpreter','latex','Fontsize',15)
xlabel('$$x_k$$','Interpreter','latex','Fontsize',15)
zlabel('Transition probability','Interpreter','latex','Fontsize',15)
subtightplot(2,1,2,0.05)
plot3(sweep(2:end),Tdet(3:end),0*ones(1,length(sweep(2:end))),'-k','LineWidth',3)
hold on
surf(X,Y,Tnorm(2:end,2:end),'EdgeColor','none','FaceAlpha',0.9),view(0,90);
axis tight
ylabel('$$x_{k+1}$$','Interpreter','latex','Fontsize',15)
xlabel('$$x_k$$','Interpreter','latex','Fontsize',15)

set(gcf,'Renderer','Painter')

