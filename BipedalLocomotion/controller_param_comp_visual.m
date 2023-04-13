%% Look how we are comparing different controllers for different levels of noise
clear all
clc

dilsad_thesis_fig_defaults()
close all
marker_list={'ks-','ko-','kv-','k*-','k--'};
figure()   
noise_var_sweep=[1e-4*(0.5:0.1:1.3) 1e-4*(1.4:0.2:2.2) 1e-4*(2.5:0.5:7) 1e-4*(8:1:11)];
marker_list={'s-','o-','v-','*-','d-'};

for cont_num=1:5
    mfpt=[];
    %     switch cont
    %         case 1
    %         case 2
    %             noise_var_sweep=1e-4*(0.5:0.1:1.3);
    %         case 3
    %             noise_var_sweep=1e-4*(0.5:0.1:1.3);
    %         case 4
    %             noise_var_sweep=1e-4*(0.5:0.1:1.3);
    %         case 5
    %             noise_var_sweep=1e-4*(0.5:0.1:1);
    %     end
    
    
    for noise_var=noise_var_sweep
        
        disp(cont_num)
        disp(noise_var)
        
        savename='Controller_Comparison_mat_files/Tunsc_'+string(noise_var)+'_C'+string(cont_num)+'.mat';
        marker=marker_list{cont_num};
        
        
        load(savename)
        [VVV,DDD]=eig(Tunsc(2:end,2:end));
        %DDD(DDD>=1)=0;
        % disp('Calculated MFPT')
        %         if cont ==5
        %             if abs(max(max(DDD)))>1
        %                 mfpt=[mfpt 1/(1-(1-(max(max(DDD))-1)))];
        %             else
        %                         mfpt=[mfpt 1/(1-max(max(DDD)))];
        %
        %             end
        %         else
        
        max(max(DDD))
        max_lamb=max(max(DDD));
        mfpt=[mfpt 1/(1- max_lamb)]; 
        
        % end
    end
    A=find(isinf(mfpt))-1;
    if ~isempty(A)
    if A(1)== 0
        mfpt(1)=1e15;
    end
    end
    mfpt(find(isinf(mfpt)))=   (mfpt(find(isinf(mfpt))-1) + mfpt(find(isinf(mfpt))+1))/2
    %abs(mfpt)
    plot(sqrt(noise_var_sweep),abs(mfpt),marker,'Linewidth',2)
    hold on
    set(gca, 'YScale', 'log')
    
    
    
    
    
end
grid on
legend('$$C_1$$','$$C_2$$','$$C_3$$','$$C_4$$','$$C_5$$','Interpreter','latex','Location','southwest')
xlabel('Standard Deviation of Noise','Interpreter','latex','Fontsize',12)
ylabel('MFPT of 5-link bipedal','Interpreter','latex','Fontsize',12)
title('Comparison of different controllers','Interpreter','latex','Fontsize',12)
xlim([sqrt(min(noise_var_sweep)) sqrt(max(noise_var_sweep))])

%%
noise1_sweep=(-0.1:1e-4:(0.1-1e-4))/2;
noise1_mean=0;
noise_vars=1e-4*(0.5:0.1:1.3);
figure()
for noise1_var=noise_vars
    discrete_noise1_pdf=normpdf(noise1_sweep,noise1_mean,sqrt(noise1_var));
    discrete_noise1_pdf=discrete_noise1_pdf/sum(discrete_noise1_pdf);
    plot(noise1_sweep,discrete_noise1_pdf)
    hold on
end
grid on
title('Body Angle Noise')
%%
close all
marker_list={'s-','o-','v-','*-','--'};

for cont_num=1:5
mfpt=[];
marker=marker_list{cont_num};
figure(cont_num)
for i=1:9
    
    noise_var=noise_var_sweep(i);
    savename='Controller_Comparison_mat_files/Tunsc_'+string(noise_var)+'_C'+string(cont_num)+'.mat';
        
        
    load(savename)
    subplot(3,3,i)
   % surf(Tunsc(2:end,2:end),'EdgeColor','none')
    D=sort(eig(Tunsc),'descend'); 
    view(0,90)

    1-D(2)
    1-D(3)
    
    mfpt=[mfpt 1/(1- D(2))]; 

end
    figure(6)
    mfpt(isinf(mfpt))=1e16;
    plot(noise_var_sweep,abs(mfpt),marker)
    hold on
    set(gca, 'YScale', 'log')
end
legend