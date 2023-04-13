clear all
walk = load("periodicWalking_0.4m_0.5s_noPD.mat");
[x0m,a,tp,tm] = decode(walk.optimizedStack);
close all;
x0 = x0m';
stepCounter = 100;
uAll = [];
hAll = [];
tAll = [];
xAll = [];
FnAll = [];
FtAll = [];
vAll = [];
%figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
p0 = [0,0];
tTotal = 0;
tInc = 0;
frameNo = 0;
%clear frame
noise_mean=0;
noise_var=1e-4;
Data=[];

for k = 1:stepCounter
    %disturbance=normrnd(noise_mean, sqrt(noise_var),[1,10]);
    disturbance=[normrnd(noise_mean, sqrt(noise_var),[1,5])*0,normrnd(noise_mean, sqrt(noise_var),[1,5])];
    if isempty(x0)
        break
    end
    x0=x0+disturbance;
    x0p = (impactModel(x0'))';
    options = odeset('RelTol',1e-8,'MaxStep',1e-2, 'Events', @(t,x)eventFunc(t,x));
    [t,x,te,xe,ie] = ode45(@(t,x)xDotFunc(t,x,a,tp,tm), 0:0.002:10, x0p, options);
    x0 = xe;
%     disp(xe);
%     disp(te);
    Data=[Data,xe'];
%      for i = 1:length(x(:,1))
%         uAll(end+1,:) = (u96Func(x(i,:)',a,tp,tm))';
%         hAll(end+1,:) = (hFunc(x(i,:)', a, tp, tm))';
%         xAll(end+1,:) = x(i,:);
%         tTotal = tInc + t(i);
%         tAll(end+1,:) = tTotal;
%         F = conForFunc(x(i,:)',a,tp,tm);
%         FtAll(end+1) = F(1);
%         FnAll(end+1) = F(2);
%         
%         jp = jointPointsFunc(x(i,:)') + p0;
%         cmp = cmPointsFunc(x(i,:)') + p0;
%         clf;
%         plot(jp(:,1),jp(:,2), 'LineWidth', 2, 'Color', 'blue');
%         line([-100,100], [0,0], 'LineWidth', 2, 'Color', 'black');
%         xlim([-0.5,4.0])
%         ylim([-0.5,2])
      %  text(1.65,1.6,"time: " + num2str(tTotal))
        %drawnow();
%     end
    tInc = tInc + te;
end
%%
% Data=Data(:,10:end);
Aug_Data=[x0m,Data(:,1:end-1)];