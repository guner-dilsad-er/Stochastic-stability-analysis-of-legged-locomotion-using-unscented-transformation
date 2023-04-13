close all
distq5 = -0.1:0.002:0.1;
distq10 = -1:0.02:1;
load('boaPD01_1_002_befIm.mat')
figure()
X = distq5+0.014721768962033;
Y = distq10-0.848597863783639;
[Xi, Yi] = meshgrid(X,Y);
surf(Xi,Yi,boa', 'EdgeColor', 'none');
view(0,90)
axis tight
xlabel('$$q_5$$','interpreter','latex')
ylabel('$$\dot{q}_5$$','interpreter','latex')
hold on
lins=linspace(2,length(boa(:,1))-1,15);
boa_modified=boa;
boa_modified(boa<10)=nan;
mesh(Xi(lins,lins),Yi(lins,lins),boa_modified(lins,lins)'+0.01,'FaceColor','none','EdgeColor','k')
colormap Parula;

boa_d=boa(lins,lins);
ten=(boa_d==10);
idx_list=[];
for cols=1:length(ten(:,1))
    ten_idx=find(ten(:,cols));
    if ~isempty(ten_idx)
        idx_list=[idx_list; ten_idx ,cols*ones(length(ten_idx),1)];
    end
end
distq5_d=distq5(2:7:100);
distq10_d=distq10(2:7:100);
states=[];
for i=1:length(idx_list(:,1))
    states=[states;distq5_d(idx_list(i,1)),distq10_d(idx_list(i,2))];
end
states = sortrows( states ,2 ) +[0.014721768962033,-0.848597863783639];
for i=1:length(states(:,1))
    plot3(states(i,1),states(i,2),10+0.01,'r*')
end
%%
shp = alphaShape(states(:,1),states(:,2));
k=convhull(states(:,1),states(:,2));
figure()
%subplot(1,2,2)
plot(shp)
hold on
plot(states(k,1),states(k,2))
xlabel('$$q_5$$','interpreter','latex')
ylabel('$$\dot{q}_5$$','interpreter','latex')
axis square
samplingnumber=100000;

xq = unifrnd(min(distq5)+0.014721768962033, max(distq5)+0.014721768962033,samplingnumber,1);
yq = unifrnd(min(distq10)-0.848597863783639, max(distq10)-0.848597863783639,samplingnumber,1);
xv=states(k,1);
yv=states(k,2);
in = inpolygon(xq,yq,xv,yv);
plot(xq(in),yq(in),'r+')
plot(xq(~in),yq(~in),'bo') % points outside
axis tight

disp('Successful Sampling: '+string(numel(xq(in)))+' out of '+string(samplingnumber))

initial_conditions=[xq(in),yq(in)];
%%

sweep=states;
addpath 'code/'
walk = load("periodicWalking_0.4m_0.5s_noPD.mat");
[x0m,a,tp,tm] = decode(walk.optimizedStack);
x0in= x0m';

f=@one_step_bipedal_for_noise_experiments;
% exp_num=100000;%10000;
%%
% Choose an initial state and noise variance:
noise_var=1e-3;
noise_mean=0;

%% Find output distributions experimentally
%Monte Carlo Experiment
noiseMem=[];
NoisyStateTransition=[];%zeros(exp_num,20);

for j=1:length(initial_conditions(:,1))
    x0in(10)=initial_conditions(j,2);
    x0in(5)=initial_conditions(j,1);
    disturbance=[normrnd(noise_mean, sqrt(noise_var),[1,5])*0,normrnd(noise_mean, sqrt(noise_var),[1,5])];%normrnd(noise_mean, sqrt(noise_var),[1,10]);
    noiseMem(j,:)= disturbance;
    [xe]=one_step_bipedal_for_noise_experiments(x0in,disturbance);
    if isempty(xe) %|| ~inpolygon(xe(5),xe(10),xv,yv)
        NoisyStateTransition=[NoisyStateTransition;x0in,x0in*nan];
        %         continue
    else
        % Store information in State Transition Map
        xe=xe(1,:);
    
        NoisyStateTransition=[NoisyStateTransition;x0in,xe];
    end
end
disp('30k exp done')
%%
save('TWOd_TWO_D_inputlimit.mat', 'NoisyStateTransition','noiseMem')
%load('TWOd_TWO_D_56653.mat')
%%
T=zeros(length(states)+1,length(states)+1);
for j=1:length(NoisyStateTransition(:,1))
    h=NoisyStateTransition(j,[5,10]);
    h_next=NoisyStateTransition(j,[15,20]);
    in = inpolygon(h(1),h(2),xv,yv);
    if ~in
        out1=1;
        disp('Unexpected Error')
    else
        out1=find(min(sqrt(sum((states-h).^2,2))) ==sqrt(sum((states-h).^2,2)),1)+1;
    end
    in = inpolygon(h_next(1),h_next(2),xv,yv);
    if ~in
        out2=1;
    else
        out2=find(min(sqrt(sum((states-h_next).^2,2))) ==sqrt(sum((states-h_next).^2,2)),1)+1;
    end
    T(out2,out1)=T(out2,out1)+1;%*discrete_noise2_pdf(k);
end

%%
T(1)=1;
Tmc=T./sum(T,1);
disp('MC done')
%%
%

distq5 = -0.1:0.002:0.1;
distq10 = -1:0.02:1;
load('boaPD01_1_002_befIm.mat')
dilsad_thesis_fig_defaults()
figure();
%subplot(1,2,1)
X = distq5+0.014721768962033;
Y = distq10-0.848597863783639;
[Xi, Yi] = meshgrid(X,Y);
surf(Xi,Yi,boa', 'EdgeColor', 'none','FaceAlpha',0.5);
view(0,90)
axis tight
hold on
for j=1:1000:length(NoisyStateTransition(:,1))
    h=[NoisyStateTransition(j,[5,10]),10.01];
    h_next=[NoisyStateTransition(j,[15,20]),10.01];
    vectarrow(h, h_next)
    hold on
end
xlabel('$$q_5$$','interpreter','latex')
ylabel('$$\dot{q}_5$$','interpreter','latex')
%%
figure()

surf(Tmc(2:end,2:end))
axis tight
view(0,90)
titlestr='Monte Carlo Method with '+string(length(initial_conditions(:,1)))+' experiments';
title(titlestr,'interpreter','latex')
xlabel('$$(q^5_k,\dot{q}^5_k)$$','Interpreter','latex')
ylabel('$$(q^5_{k+1},\dot{q}^5_{k+1})$$','Interpreter','latex')
colorbar()

for i=1:length(states(:,1))
    tick_labels{i}='('+string(states(i,1))+','+string(states(i,2))+')';
end
xticks(1:71)
xticklabels(tick_labels)
yticks(1:71)
yticklabels(tick_labels)
xtickangle(45)

%%
distq5l = (-1+0.006):0.014:1;
distq10l = (-10+0.06):0.14:10;
Xl = distq5l+0.014721768962033;
Yl = distq10l-0.848597863783639;
[Xll, Yll] = meshgrid(Xl,Yl);
all_states=[];
for j=1:length(Xll(:))
all_states=[all_states;Xll(j) ,Yll(j) ];
end
figure()
plot(all_states(:,1),all_states(:,2),'*'),hold on
plot(states(:,1),states(:,2),'b*')
%%
global cont_num
cont_num=1;
T=zeros(length(states)+1,length(states)+1);
f=@one_step_bipedal_for_noise_experiments;
idx=1;

for idx=1:length(states(:,1))
    
    P0=eps;
    Rw=noise_var;
    Rwp=eps;
    P=diag([P0,P0,P0,P0,P0,P0,P0,P0,P0,P0,Rwp,Rwp,Rwp,Rwp,Rwp,Rw,Rw,Rw,Rw,Rw]);
    x0in(5)=states(idx,1);
    x0in(10)=states(idx,2);
    x_aug=[ x0in, zeros(1,10)];
    
    %Sigma points around reference point
    %Inputs:
    %       x: reference point
    %       P: covariance
    %       c: coefficient
    %Output:
    %       X: Sigma points
    
    n= length(x_aug);
    W0=-2;
    L=2*n+1;%numel(x_aug);
    alpha=2e-3;                                 %default, tunable
    ki=0;                                       %default, tunable
    beta=1;                                     %default, tunable
    lambda=1;%alpha^2*(L+ki)-L;                    %scaling factor
    c=n+lambda;                                 %scaling factor
    Wm=[lambda/c 0.5/c+zeros(1,2*n)];           %weights for means
    
    Wm=Wm/sum(Wm);
    
    Wc=Wm;
    Wc=Wc/sum(Wc);
    Wm=Wc;
    A = sqrt(n/(1-W0))*chol(P)';
    x_a1=x_aug';
    X = [x_a1, x_a1+A ,x_a1-A];
    y=zeros(n,1);
    Y=zeros(n,L);
    
    for k=1:L
        if isempty(f(X(1:10,k)',X(11:20,k)'))
            continue
        end
        Y(:,k)=[f(X(1:10,k)',X(11:20,k)'),zeros(1,10)]';
        y=y+Wm(k)*Y(:,k);
    end
    %x=f(x1, w_k, v_k);
    Y1=Y-y(:,ones(1,L));
    P=Y1*diag(Wc)*Y1';
    mu= y([5,10])';
    Sigma=[P(5,5), P(10,10)];
    MVPDF = mvnpdf(states,mu,Sigma);
    all_mvpdf=mvnpdf(all_states,mu,Sigma);

    disp(sum(all_mvpdf))
    disp(sum(MVPDF))
    if(sum(MVPDF)>sum(all_mvpdf))
        disp('Normalization Error')
        disp(mu)
        disp(Sigma)
    end
    T(2:end,idx+1)=MVPDF'/sum(all_mvpdf);
    disp(idx+1)
end
%%
T(1)=1;
Tunsc=T;
%%
figure()
surf(Tunsc(2:end,2:end))
axis tight
view(0,90)
titlestr='Proposed Method with '+string(41*length(states(:,1)))+' experiments';
title(titlestr,'interpreter','latex')
xlabel('$$(q^5_k,\dot{q}^5_k)$$','Interpreter','latex')
ylabel('$$(q^5_{k+1},\dot{q}^5_{k+1})$$','Interpreter','latex')
colorbar()
caxis([0 1])
for i=1:length(states(:,1))
    tick_labels{i}='('+string(states(i,1))+','+string(states(i,2))+')';
end
xticks(1:71)
xticklabels(tick_labels)
yticks(1:71)
yticklabels(tick_labels)
xtickangle(45)
%%
[Vmc,Dmc]=eig(Tmc);
[Vut,Dut]=eig(Tunsc);