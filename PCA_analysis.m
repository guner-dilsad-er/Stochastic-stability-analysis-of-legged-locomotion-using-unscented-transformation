clear all
close all;
addpath 'bipedal_simulation/'
walk = load("periodicWalking_0.4m_0.5s_noPD.mat");
[x0m,a,tp,tm] = decode(walk.optimizedStack);
x0 = x0m';
stepCounter = 100;
global cont_num
cont_num=1;
p0 = [0,0];
tTotal = 0;
tInc = 0;
frameNo = 0;
noise_mean=0; noise_var=1e-3;
%%
Data=[];
DisturbanceData=[];
for k = 1:stepCounter
    disturbance=[normrnd(noise_mean, sqrt(noise_var),[1,5])*0,normrnd(noise_mean, sqrt(noise_var),[1,5])];
    if isempty(x0)
        break
    end
    x0=x0+disturbance;
    x0p = (impactModel(x0'))';
    options = odeset('RelTol',1e-8,'MaxStep',1e-2, 'Events', @(t,x)eventFunc(t,x));
    [t,x,te,xe,ie] = ode45(@(t,x)xDotFunc(t,x,a,tp,tm), 0:0.002:10, x0p, options);
    if isempty(xe)
        break
    else        
        p2=p2Func(xe');
        jp = jointPointsFunc(xe');
        if p2(1)<0.1 || p2(2)<-0.001 || jp(2,2)<0 || jp(5,2)<0.4 % step is shorter than 0.1m, toe inside ground
            break        
        else        
            x0 = xe;
            Data=[Data,xe'];
            DisturbanceData=[DisturbanceData; disturbance];
        end
    end
    tInc = tInc + te;
end
%%
Aug_Data=[x0m,Data(:,1:end-1)];
DevData=Data-Aug_Data;
DevData(:,1)=[];
[coeff,score,latent,tsquared,explained,mu] = pca(Data(1:10,:)');
Xcentered = score*coeff'+mu;
[RESIDUALS,RECONSTRUCTED] = pcares(Data,1);
%%
figure()
for i=1:10
    subplot(5,2,i)
    plot(DevData(i,:))
    hold on
    plot(DisturbanceData(:,i))
    grid on
    xlabel('Step Number')
    ylabel('Deviation')
end
sgtitle('Deviation Data')
legend('DevData','Location','south')
%%
vec=latent;
attrib='Scree Plot';
lev=0;
figure()
subplot(2,2,3)
screeplot(vec,attrib,lev)
subplot(2,2,[1 2])
biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10'});%
title('PCA biplot')
subplot(2,2,4)
plot(score(:,1),score(:,2),'x')
grid on, hold on
title('Score Plot')
sgtitle('Principal Component Analysis','interpreter','latex')
%%
figure()
for k=1:15
    subplot(5,3,k)
    i=randi(10,1,1);
    j=randi(10,1,1);
    if i==j
         i=randi(10,1,1);
         j=randi(10,1,1);
    end
    plot(DevData(i,:),DevData(j,:),'x')
    hold on
    grid on
    xl='q'+string(i);
    yl='q'+string(j);
    xlabel(xl)
    ylabel(yl)
    axis equal
end
sgtitle('Different Correlations')
