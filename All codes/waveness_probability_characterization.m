N_sh=500;
n_sessions=size(prob_data,1);

for i=1:n_sessions
    waveness_prob(i)=sum(prob_data(i,3:end));
    
    for sh=1:N_sh
        waveness_prob_sh(i,sh)= sum(prob_data_sh_allsessions{1,i}(sh,3:end));
    end
end

thr=prctile(waveness_prob_sh',99);

sig_waveness=find(waveness_prob>thr);

waveness_prob_PaS.waveness=waveness_prob;
waveness_prob_PaS.waveness_sh=waveness_prob_sh;
waveness_prob_PaS.sig_waveness=sig_waveness;

%% Waveness for MEC sessions without waves

clear all
close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

mice_number=12;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];

ml=[2.64,3.7,2.7,2.64,3.2,3.3,2.52,3.3,3.5,-10,-10,-10]; %ML coordinates

count=0;
for m=1:mice_number
    
    if mice(m,1)~= 'L'
        mouse=mice(m,:);
    else
        if mice(m,2)=='0'
            mouse=[mice(m,1),mice(m,3:5)];
            mouse_l=mouse(2);
            mouse_a=mouse(4);
        else
            mouse=mice(m,:);
            mouse_l=mouse(2:3);
            mouse_a=mouse(5);
        end
    end
    
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    load([save_data_path ['WS_Osc_14_',mice(m,:),'.mat']]);
    ws_ent=load([save_data_path ['WS_Entropy_',mice(m,:),'.mat']]);
    
    for day=1:dates.daysnum
        for s=1:dates.sesnum(day)
            disp(s)
            munit=dates.ses{day}(s);
            
            file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
            
            if (exist(file_name_spk) == 2)
                disp(day)
                spk=load(file_name_spk);
                
                if s<100
                    count=count+1;
                    
                    if isfield(dates,'actual_day')  == 1
                        big_table(count,3)=dates.actual_day(day);   %Day # on the wheel
                    else
                        big_table(count,3)=day; %Day # on the wheel
                    end
                    
                    if mouse(1) == '9' || mouse(1) == '6'
                        big_table(count,1)=1;%str2num(mice(m,2:3));
                        big_table(count,10)=-10; %ML position
                        if mouse=='92227'
                            big_table(count,2)=1;
                        elseif mouse == '92229'
                            big_table(count,2)=2;
                        elseif mouse == '60961'
                            big_table(count,2)=3;
                        end
                    else
                        big_table(count,1)=str2num(mice(m,2:3)); %Litter
                        big_table(count,2)=str2num(mice(m,end)); %Mouse number
                        big_table(count,10)=ml(m); %ML position
                    end
                    
                    big_table(count,4)=s; %Session
                    big_table(count,5)=munit; %munit
                    big_table(count,6)=WS_stat.WS(day,s); %WS - Osc
                    big_table(count,11)=ws_ent.WS_stat.wave_score_ent(day,s); %WS - Ent
                    big_table(count,12)=size(spk.spikes_d_s,1);
                    
                    threshold_kl=1;
                    if ( big_table(count,6)>=threshold_kl )
                        big_table(count,7)=1;
                    else
                        big_table(count,7)=0;
                    end
                    
                    if (isnan(big_table(count,6)))
                        big_table(count,7)=NaN;
                        big_table(count,8)=NaN; %dt
                        big_table(count,9)=NaN; %fft
                    else
                        big_table(count,8)=WS_stat.dt(day,s); %dt
                        big_table(count,9)=WS_stat.FFT(day,s); %fft
                    end
                    
                end
                clear spk
            end
        end
        
    end
    clear WS_stat
end

% No wave sessions

MEC_ses=find(big_table(:,10)>3); %Sessions inthe in MEC
adults=find(big_table(:,3)>15); %Sessions inthe of adults
sessions=intersect(MEC_ses,adults);
count=0;

ws=big_table(sessions,6);
index_nan=isnan(ws);
sessions(index_nan)=[];
session_number=big_table(sessions,4);
days=big_table(sessions,3);
repeated_days=find(diff(days)==0);
clear index_to_delete

count=0;
for w=1:length(repeated_days)
    ws1=big_table(sessions(repeated_days(w)),6);
    ws2=big_table(sessions(repeated_days(w)+1),6);
    
    if (ws1==1 && ws2==0)
        count=count+1;
        index_to_delete(count)=repeated_days(w)+1;
    elseif  (ws1==0 && ws2==1)
        count=count+1;
        index_to_delete(count)=repeated_days(w);
    elseif (ws1==0 && ws2==0)
        count=count+1;
        index_to_delete(count)=repeated_days(w);
    else
        print('Other case');
    end
end

sessions(index_to_delete)=[];

big_table_2=big_table(sessions,:);

mec_sessions=sessions;
no_waves=mec_sessions(find(big_table(mec_sessions,6)==0));

% Calculates the transition matrices
N_sh=500;
clus=10;
count=0;

for w=1:length(no_waves)
    row_w=no_waves(w);
    disp(w);
    
    count=count+1;
    mouse=['L',num2str(big_table(no_waves(w),1)),'M',num2str(big_table(no_waves(w),2))];
    day=big_table(no_waves(w),3);
    s=big_table(no_waves(w),4);
    munit=big_table(no_waves(w),5);
    
    %     dt=floor(big_table(no_waves(w),8));
    %     if isinteger(dt)
    %     else
    %         dt=floor(dt);
    %     end
    
   
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    if isfield(dates,'actual_day')==1
        day=find(dates.actual_day==day);
    end
    
    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_spk,'-mat');
    spikes=full(spikes_d_s);
    
    spikes_w=spikes;
    
    %Calculates the probability of wave length
    
    dt=66;
    
    
    [prob_data(w,:)] = comp_wave_prob_f(spikes_w,floor(dt),clus);
    
    for i=1:N_sh
        mat_sh=shuffle(spikes_w')';
        [prob_data_sh(i,:)]=comp_wave_prob_f(mat_sh,floor(dt),clus);
        clear mat_sh
    end
    prob_data_sh_allsessions{count}=prob_data_sh;
    
    %TM
    %     [TM(:,:,count),adj_sh] = comp_TM_quick(spikes_w,floor(dt),clus,N_sh);
    %     TM_sh_allsessions{count}=adj_sh;
    
    
    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2
    
end

%% comparison

fraction_mec_w=length(waveness_prob_MEC_waves.sig_waveness)/length(waveness_prob_MEC_waves.waveness);
fraction_mec_nw=length(waveness_prob_MEC_no_waves.sig_waveness)/length(waveness_prob_MEC_no_waves.waveness);
fraction_PaS=length(waveness_prob_PaS.sig_waveness)/length(waveness_prob_PaS.waveness);
fraction_V1=length(waveness_prob_V1.sig_waveness)/length(waveness_prob_V1.waveness);

mean_waveness_mec_w=mean(waveness_prob_MEC_waves.waveness(waveness_prob_MEC_waves.sig_waveness));
mean_waveness_mec_nw=mean(waveness_prob_MEC_no_waves.waveness(waveness_prob_MEC_no_waves.sig_waveness));
mean_waveness_pas=mean(waveness_prob_PaS.waveness(waveness_prob_PaS.sig_waveness));
mean_waveness_v1=mean(waveness_prob_V1.waveness(waveness_prob_V1.sig_waveness));


sem_waveness_mec_w=std(waveness_prob_MEC_waves.waveness(waveness_prob_MEC_waves.sig_waveness))/sqrt(length(waveness_prob_MEC_waves.sig_waveness));
sem_waveness_mec_nw=std(waveness_prob_MEC_no_waves.waveness(waveness_prob_MEC_no_waves.sig_waveness))/sqrt(length(waveness_prob_MEC_no_waves.sig_waveness));
sem_waveness_pas=std(waveness_prob_PaS.waveness(waveness_prob_PaS.sig_waveness))/sqrt(length(waveness_prob_PaS.sig_waveness));
sem_waveness_v1=std(waveness_prob_V1.waveness(waveness_prob_V1.sig_waveness))/sqrt(length(waveness_prob_V1.sig_waveness));

% 
% figure
% bar(100*[fraction_mec_w,fraction_mec_nw,fraction_PaS,fraction_V1])
% ylabel({'Sessions with'; 'significant waveness %'})
% axis([0.5 4.5 0 110])
% box off
% set(gca,'fontsize',18)
% xticklabels({'MEC/W','MEC/NW' , 'PaS', 'V1'})
% 
% figure
% bar(100*[fraction_mec_w,fraction_PaS,fraction_V1])
% ylabel({'Sessions with'; 'significant waveness %'})
% axis([0.5 3.5 0 110])
% box off
% set(gca,'fontsize',18)
% xticklabels({'MEC/W','PaS', 'V1'})

figure
bar(100*[fraction_PaS,fraction_V1,fraction_mec_w])
ylabel({'Sessions with'; 'significant waveness %'})
axis([0.5 3.5 0 110])
box off
set(gca,'fontsize',18)
xticklabels({'PaS', 'V1','MEC/W'})



figure
hold on
bar([mean_waveness_pas,mean_waveness_v1,mean_waveness_mec_w])
errorbar([1,2,3],[mean_waveness_pas,mean_waveness_v1,mean_waveness_mec_w],[sem_waveness_pas,sem_waveness_v1,sem_waveness_mec_w],'k*','linewidth',2)
ylabel('Waveness')
axis([0.5 3.5 0 1])
box off
set(gca,'fontsize',18)
xticks([1 2 3])
xticklabels({'PaS', 'V1','MEC/W'})


[p,h_,stats]=ranksum(waveness_prob_MEC_no_waves.waveness(waveness_prob_MEC_no_waves.sig_waveness),waveness_prob_PaS.waveness(waveness_prob_PaS.sig_waveness));


figure
bar(100*[fraction_mec_w,fraction_mec_nw])
ylabel({'Sessions with'; 'significant waveness %'})
axis([0.5 2.5 0 110])
box off
set(gca,'fontsize',18)
xticklabels({'MEC/W','MEC/NW'})

figure
hold on
bar([mean_waveness_mec_w,mean_waveness_mec_nw])
errorbar([1,2],[mean_waveness_mec_w,mean_waveness_mec_nw],[sem_waveness_mec_w,sem_waveness_mec_nw],'k*','linewidth',2)
ylabel('Waveness')
axis([0.5 2.5 0 1])
box off
set(gca,'fontsize',18)
xticks([1 2 3])
xticklabels({'MEC/W','MEC/NW'})

