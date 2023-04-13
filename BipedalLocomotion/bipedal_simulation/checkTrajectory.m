walk = load("periodicWalking_0.4m_0.5s_noPD.mat");
[x0m,a,tp,tm] = decode(walk.optimizedStack);
close all;
x0 = x0m';
%x0(5)=0.045+0.07
%test 
%test2
stepCounter = 1;
uAll = [];
hAll = [];
tAll = [];
xAll = [];
FnAll = [];
FtAll = [];
vAll = [];
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
p0 = [0,0];
tTotal = 0;
tInc = 0;
frameNo = 0;
clear frame

for k = 1:stepCounter
    x0p = (impactModel(x0'))';
   %x0p(5)=x0p(5)+0.05;
  %   x0p(5)=0.045+0.07
    options = odeset('RelTol',1e-8,'MaxStep',1e-2, 'Events', @(t,x)eventFunc(t,x));
    [t,x,te,xe,ie] = ode45(@(t,x)xDotFunc(t,x,a,tp,tm), 0:0.002:10, x0p, options);
    x0 = xe;
   
    disp(xe);
    disp(te);
    for i = 1:length(x(:,1))
        uAll(end+1,:) = (u96Func(x(i,:)',a,tp,tm))';
        hAll(end+1,:) = (hFunc(x(i,:)', a, tp, tm))';
        xAll(end+1,:) = x(i,:);
        tTotal = tInc + t(i);
        tAll(end+1,:) = tTotal;
        F = conForFunc(x(i,:)',a,tp,tm);
        FtAll(end+1) = F(1);
        FnAll(end+1) = F(2);
        
        jp = jointPointsFunc(x(i,:)') + p0;
        cmp = cmPointsFunc(x(i,:)') + p0;
        clf;
        plot(jp(:,1),jp(:,2), 'LineWidth', 2, 'Color', 'blue');
        line([-100,100], [0,0], 'LineWidth', 2, 'Color', 'black');
        xlim([-0.5,4.0])
        ylim([-0.5,2])
        text(1.65,1.6,"time: " + num2str(tTotal))
        drawnow();
    end
    p0 = jp(end,:);
    tInc = tInc + te;
end

figure();
subplot(4,1,1)
hold on;
plot(tAll, uAll(:,1), 'LineWidth', 2);
ylabel('u1 (Nm)');
xlabel('t');

subplot(4,1,2)
hold on;
plot(tAll, uAll(:,2), 'LineWidth', 2);
ylabel('u2 (Nm)');
xlabel('t');

subplot(4,1,3)
hold on;
plot(tAll, uAll(:,3), 'LineWidth', 2);
ylabel('u3 (Nm)');
xlabel('t');

subplot(4,1,4)
hold on;
plot(tAll, uAll(:,4), 'LineWidth', 2);
ylabel('u4 (Nm)');
xlabel('t');
sgtitle('Input Torques');

figure();
subplot(5,1,1)
hold on;
plot(tAll, xAll(:,1), 'LineWidth', 2);
ylabel('q1 (rad)');
xlabel('t');

subplot(5,1,2)
hold on;
plot(tAll, xAll(:,2), 'LineWidth', 2);
xlabel('t');
ylabel('q2 (rad)');

subplot(5,1,3)
hold on;
plot(tAll, xAll(:,3), 'LineWidth', 2);
xlabel('t');
ylabel('q3 (rad)');

subplot(5,1,4)
hold on;
plot(tAll, xAll(:,4), 'LineWidth', 2);
xlabel('t');
ylabel('q4 (rad)');

subplot(5,1,5)
hold on;
plot(tAll, xAll(:,5), 'LineWidth', 2);
xlabel('t');
ylabel('q5 (rad)');
sgtitle('Joint Positions');

figure();
subplot(5,1,1)
hold on;
plot(tAll, xAll(:,6), 'LineWidth', 2);
xlabel('t');
ylabel('dq1 (rad/s)');

subplot(5,1,2)
hold on;
plot(tAll, xAll(:,7), 'LineWidth', 2);
xlabel('t');
ylabel('dq2 (rad/s)');

subplot(5,1,3)
hold on;
plot(tAll, xAll(:,8), 'LineWidth', 2);
xlabel('t');
ylabel('dq3 (rad/s)');

subplot(5,1,4)
hold on;
plot(tAll, xAll(:,9), 'LineWidth', 2);
xlabel('t');
ylabel('dq4 (rad/s)');

subplot(5,1,5)
hold on;
plot(tAll, xAll(:,10), 'LineWidth', 2);
xlabel('t');
ylabel('dq5 (rad/s)');
sgtitle('Joint Velocities')

figure();
subplot(2,1,1)
hold on;
plot(tAll, cFunc()*xAll(:,1:5)', 'LineWidth', 2);
xlabel('t');
ylabel('theta (rad)');

subplot(2,1,2)
hold on;
plot(tAll, cFunc()*xAll(:,6:10)', 'LineWidth', 2);
xlabel('t');
ylabel('dtheta (rad/s)');
sgtitle('Internal Clock');

figure();
subplot(3,1,1)
hold on;
plot(tAll, FtAll, 'LineWidth', 2);
xlabel('t');
ylabel('Ft (N)');

subplot(3,1,2)
hold on;
plot(tAll, FnAll, 'LineWidth', 2);
xlabel('t');
ylabel('Fn (N)');

subplot(3,1,3)
hold on;
plot(tAll, abs(FtAll./FnAll), 'LineWidth', 2);
xlabel('t');
ylabel('|Ft/Fn|');
sgtitle('Contact Forces');
