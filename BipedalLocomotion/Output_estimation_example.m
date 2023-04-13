%% Look how good we are estimating the output disturbances:
clear all
clc
%%
increment=0.004;
states=(-0.28:increment:0.18)-0.848597863783639;
all_states=-10:increment:10;
sweep=states+increment/2;
sweep(end)=[];
walk = load("periodicWalking_0.4m_0.5s_noPD.mat");
[x0m,a,tp,tm] = decode(walk.optimizedStack);
x0in= x0m';
f=@one_step_bipedal_for_noise_experiments;
exp_num=10000;
%%
% Choose an initial state and noise variance:
k=70;
noise_var=1e-3;
noise_mean=0;

%% Find output distributions experimentally
% Monte Carlo Experiment
noiseMem=[];
NoisyStateTransition=zeros(exp_num,20);
x0in(10)=sweep(k);
parfor j=1:exp_num
    disturbance=normrnd(noise_mean, sqrt(noise_var),[1,10]);
    noiseMem(j,:)= disturbance;
    [xe]=one_step_bipedal_for_noise_experiments(x0in,disturbance);
    if isempty(xe)
        continue
    else
        % Store information in State Transition Map
        xe=xe(1,:);
        NoisyStateTransition(j,:)=[x0in,xe];
    end
    disp(string(k)+','+string(j))
end
%% Find output distributions with unscented transform
y_0=sweep(k);
P0=eps;
Rw=noise_var;
 Rwp=eps;
P=diag([P0,P0,P0,P0,P0,P0,P0,P0,P0,P0,Rwp,Rwp,Rwp,Rwp,Rwp,Rw,Rw,Rw,Rw,Rw]);
x0in(10)=y_0;
x_aug=[ x0in, zeros(1,10)];
%Sigma points around reference point
%Inputs:
%       x: reference point
%       P: covariance
%       c: coefficient
%Output:
%       X: Sigma points
n= length(x_aug);
W0=-1;
L=2*n+1;
% alpha=2e-3;                                 %default, tunable
% ki=0;                                       %default, tunable
% beta=1;                                     %default, tunable
% lambda=1;                                   %scaling factor
% c=n+lambda;                                 %scaling factor
% Wm=[lambda/c 0.5/c+zeros(1,2*n)];           %weights for means
% Wm=Wm/sum(Wm);
% Wc=Wm;
% Wc=Wc/sum(Wc);
% Wm=Wc;

 Wm=ones(1,41);
    Wm(1)=Wm(1)+2;
    Wm=Wm/sum(Wm);

    Wc=ones(1,41);
    Wc(1)=Wc(1)-0.5;
    Wc=Wc/sum(Wc);

A = sqrt(n/(1-W0))*chol(P)';
x_a1=x_aug';
X = [x_a1, x_a1+A ,x_a1-A];
y=zeros(n,1);
Y=zeros(n,L);
for k=1:L
    Y(:,k)=[f(X(1:10,k)',X(11:20,k)'),zeros(1,10)]';
    y=y+Wm(k)*Y(:,k);
end
Y1=Y-y(:,ones(1,L));
P=Y1*diag(Wc)*Y1';
P1=P(10,10);
 z1=y;  
DP=sqrt(diag(P(1:10,1:10))');


%% Lin
Q=diag([Rwp,Rwp,Rwp,Rwp,Rwp,Rw,Rw,Rw,Rw,Rw]);

 x0in(10)=y_0;
    fstate=@(x)(one_step_bipedal_for_noise_experiments(x0in,x));
    del_f=jcbn_finite_diff(fstate,zeros(1,10));
 x_k_f=one_step_bipedal_for_noise_experiments(x0in,0);

  P_k_f=del_f*Q*del_f';
  DP2=sqrt(diag(  P_k_f(1:10,1:10))');
  

%% Visualize the results
figure()
for ind=1:10    
    subplot(2,5,ind)
    h=histogram(NoisyStateTransition(:,10+ind),'BinEdges',[(z1(ind)-5*sqrt(DP(ind))):sqrt(DP(ind))/3:(z1(ind)+5*sqrt(DP(ind)))],'Normalization', 'pdf','Normalization', 'pdf','FaceColor','white');
    hold on  
   
    bins=h.BinEdges;
    discrete_pdf=normpdf(bins,z1(ind),sqrt(DP(ind)));
    discrete_pdf=discrete_pdf;
    p1=plot(bins,discrete_pdf,'LineWidth',2);
    p1.Color='#D95319';
    e1=errorbar(z1(ind),max( discrete_pdf)/2,sqrt(DP(ind)),'horizontal','LineWidth',3);
    e1.Color='#D95319';
    
    discrete_pdf=normpdf(bins,x_k_f(ind),sqrt(DP2(ind)));
    discrete_pdf=discrete_pdf;
    p2=plot(bins,discrete_pdf,'LineWidth',2);
    p2.Color='#77AC30';
    e2=errorbar(x_k_f(ind),max( discrete_pdf)/2,sqrt(DP2(ind)),'horizontal','LineWidth',3);
    e2.Color='#77AC30';
        e=errorbar(mean(NoisyStateTransition(:,10+ind)),0,std(NoisyStateTransition(:,10+ind)),'horizontal','LineWidth',3);
    e.Color='#0072BD';
 if ind<6
    ttlstr='$q_{'+string(ind)+'}$';
    else
    ttlstr='$\dot{q}_{'+string(ind-5)+'}$';
    end
    %title(ttlstr,'Interpreter','latex','FontSize', 15)
    xlabel(ttlstr)
end

subplot(2,5,8)
legend([h,p1,p2],'experimental','unscented','linearized','Location','southoutside')


subplot(2,5,1)
ylabel('PDF')
subplot(2,5,6)
ylabel('PDF')
