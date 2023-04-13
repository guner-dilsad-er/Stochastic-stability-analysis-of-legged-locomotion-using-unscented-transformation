% 10 D noise experiments %
%clear all
%clc
addpath('code/')
global noise_var cont_num 
%%
T=zeros(length(states)-1,length(states)-1);
idx=1;
unsc_var_list=[];
unsc_mean_list=[];

for y_0=sweep(2:end)
    idx=idx+1;
    % x_aug=[x, w_k1];
    P0=eps;
    %Rw=eps;
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
    W0=-1;%-0.3;
    L=2*n+1;%numel(x_aug);
    alpha=2e-3;%0.5;                                 %default, tunable
    ki=0;                                       %default, tunable
    beta=0;                                     %default, tunable
    lambda=1;%alpha^2*(L+ki)-L;                    %scaling factor
    c=n+lambda;                                 %scaling factor
    Wm=[lambda/c 0.5/c+zeros(1,2*n)];           %weights for means
        Wm=ones(1,41);
Wm(1)=Wm(1)+2;
    Wm=Wm/sum(Wm);
    
    Wc=Wm;
    Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance
    Wc=ones(1,41);
    Wc(1)=Wc(1)-0.5;%Wc(2:10)=0;
    Wc=Wc/sum(Wc);
   % Wm=Wc;
    %A = diag(sqrt(n/(1-W0))*chol(P)');
    A = sqrt(n/(1-W0))*chol(P)';
    x_a1=x_aug';
    %X = [x_a1, x_a1+A ,x_a1-A];
    X = [x_a1, x_a1+A ,x_a1-A];
    y=zeros(n,1);
    Y=zeros(n,L);
    
    for k=1:L
          if isempty(f(X(1:10,k)',X(11:20,k)'))
                   Y(:,k)=[x0in,zeros(1,10)]';
            continue
          end
        Y(:,k)=[f(X(1:10,k)',X(11:20,k)'),zeros(1,10)]';
        y=y+Wm(k)*Y(:,k);
    end
    %x=f(x1, w_k, v_k);
    Y1=Y-y(:,ones(1,L));
    P=Y1*diag(Wc)*Y1';
    P1=P(10,10);
    z1=y(10);  
    unsc_var_list=[unsc_var_list P1];
    unsc_mean_list=[unsc_mean_list z1];
    height_pdf=normpdf(sweep,z1,sqrt(P1));
    all_height_pdf=normpdf(all_states,z1,sqrt(P1));
    height_pdf=height_pdf/sum(all_height_pdf);
    height_pdf(1)=height_pdf(1)+normcdf(sweep(1),z1,sqrt(P1))+normcdf(all_states(end),z1,sqrt(P1))-normcdf(sweep(end)+increment,z1,sqrt(P1));
    T(:,idx)=(height_pdf'/sum(height_pdf));
    %    UTpts=[UTpts;X(1,:)];
    %disp(idx)

end
%%
savename='Controller_Comparison_mat_files/Tunsc_'+string(noise_var)+'_C'+string(cont_num)+'.mat';
T(1)=1;
Tunsc=T./sum(T,1);
save(savename,'T','Tunsc','states','noise_var','unsc_var_list','unsc_mean_list')
% end
%%
