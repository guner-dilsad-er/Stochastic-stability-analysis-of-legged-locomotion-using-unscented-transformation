%% Deterministic Returnmap
clear all
clc
close all
%%

%%
pos=0.1;%0.05
pos_div=0.025/2;%2e-2;
vel=0.5;
vel_div=0.125/2;%5e-2;
x1=-pos:pos_div:pos;
x2=x1;
x3=x1;
x4=x1;
x5=x1;
x6=-vel:vel_div:vel;
x7=x6;
x8=x6;
x9=x6;
x10=x6;
%[X1,X2,X3,X4,X5,X6,X7,X8,X9,X10] = ndgrid(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10);
%%
walk = load("periodicWalking_0.4m_0.5s_noPD.mat");
addpath('../code')
[x0m,a,tp,tm] = decode(walk.optimizedStack);
x0in= x0m';
%%
% 10 by 10 mesh
% tic
% Yi=zeros(length(x1)^5*length(x6)^5,10);
% parfor i=1:length(x1)^5*length(x6)^5
%     Yi(i,:)=x0m'+[X1(i),X2(i),X3(i),X4(i),X5(i),X6(i),X7(i),X8(i),X9(i),X10(i)];    
% end
% toc

%%
 tic
Yi=[];
for s=1:10
    if s<6
        id=length(x1);
    else
        id=length(x6);
    end
    for i=1:id
        add=zeros(1,10);
        if s<6
            add(s)=x1(i);
        else
            add(s)=x6(i);
        end
        Yi=[Yi;x0m'+add];
    end
end
toc
%%
disp('Initial conditions are ready!')
%f=@one_step_bipedal_for_noise_experiments;
disturbance=zeros(1,10);
%[xe]=one_step_bipedal_for_noise_experiments(x0in,disturbance);
%%
Yend=zeros(length(Yi(:,1)),10);
for  i=1:length(Yi(:,1))
    x0in=Yi(i,:);
     disp(i)
    %[xe]=one_step_bipedal_for_noise_experiments(x0in,disturbance);
    [xe]=one_step_bipedal_for_returnmap(x0in);
    if isempty(xe)
        Yend(i,:)=ones(1,10)*nan;
        % disp(i)
    else
        p2=p2Func(xe');
        jp = jointPointsFunc(xe');
        if p2(1)<0.1 || p2(2)<-0.001 || jp(2,2)<0 || jp(5,2)<0.4 % step is shorter than 0.1m, toe inside ground
            Yend(i,:)=ones(1,10)*nan;            
        else
            Yend(i,:)=xe;
        end
    end
end
disp('Returnmap is built.')
%%
% close all

x0in= x0m';
figure()
stlist=[1,3:5,7:9];
for st=stlist
    for s=1:10
        subplot(2,5,s)
        plot(Yi((st-1)*length(x1)+1:(st)*length(x1),s),Yend((st-1)*length(x1)+1:(st)*length(x1),s),'-x'),hold on

        hold on
       
        axis tight
    end
    NaNnum=sum(isnan(Yend((st-1)*length(x1)+1:(st)*length(x1),:)));
    NaNnum(1)
end
legend(string(stlist))
%%
