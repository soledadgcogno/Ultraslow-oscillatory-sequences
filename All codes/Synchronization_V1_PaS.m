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

%% PaS sessions

PaS_ses1=find(big_table(:,10)<=3); %Sessions inthe in MEC
PaS_ses2=find(big_table(:,10)>0); %Sessions inthe in MEC
PaS_ses=intersect(PaS_ses1,PaS_ses2);
adults=find(big_table(:,3)>15); %Sessions inthe of adults
sessions=intersect(PaS_ses,adults);
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
    N1=big_table(sessions(repeated_days(w)),12);
    N2=big_table(sessions(repeated_days(w)+1),12);
    
    if (N1>N2)
        count=count+1;
        index_to_delete(count)=repeated_days(w)+1;
    elseif  (N1<N2)
        count=count+1;
        index_to_delete(count)=repeated_days(w);
    else
        print('Other case');
    end
end

sessions(index_to_delete)=[];
big_table_2=big_table(sessions,:);
PaS_sessions=sessions;

%% PaS
N_sh=500;
clus=10;
count=0;

for w=1:length(PaS_sessions)
    row_w=PaS_sessions(w);
    disp(w)
    
    mouse=['L',num2str(big_table(PaS_sessions(w),1)),'M',num2str(big_table(PaS_sessions(w),2))];
    day=big_table(PaS_sessions(w),3);
    s=big_table(PaS_sessions(w),4);
    munit=big_table(PaS_sessions(w),5);
    
    dt=66;
    
    if isnan(dt)
        disp('Error')
        continue
    else
        count=count+1;
        num_clus_discr=10;
        downsample_factor=1;
        load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
        if isfield(dates,'actual_day')==1
            day=find(dates.actual_day==day);
        end
        file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        load(file_name_spk,'-mat');
        spikes=full(spikes_d_s);       
        
        spikes_w=spikes;
        
        coac=sum(spikes_w)/size(spikes_w,1);
        
        C=corr(spikes_w');
        C_vec=C(:);
        C_vec(find(C_vec==1))=[];
               
        prob_coac_pas(w,:)=histcounts(coac,0:0.01:0.5,'Normalization','Probability');        
        prob_C_pas(w,:)=histcounts(abs(C_vec),0:0.05:1,'Normalization','Probability');
        
         N=size(spikes_w,1);
        [~,sorting_w,~]=get_sorting(spikes_w);
        
        cells_per_ens=floor(N/10);
        for e=1:10
            ens(e,:)=sorting_w((e-1)*cells_per_ens+1 : e*cells_per_ens);
            ense_n(e,:)=e*ones(1,cells_per_ens);
        end
        Ens(:,1)=reshape(ens,cells_per_ens*10,1);
        Ens(:,2)=reshape(ense_n,cells_per_ens*10,1);
        if (N-10*cells_per_ens)>0
            Ens(10*cells_per_ens+1:N,1)=sorting_w(10*cells_per_ens+1:N);
            Ens(10*cells_per_ens+1:N,2)=ones(length(sorting_w(10*cells_per_ens+1:N)),1)*10;
        end
        
        C_ens_vec=[];
        for i=1:10
            cells_ens=find(Ens(:,2)==i);
            C_ens=corr(spikes_w(Ens(cells_ens,1),:)');
            C_ens_vec=[C_ens_vec;C_ens(:)];
        end
        C_ens_vec(find(C_ens_vec==1))=[];

        prob_C_end_pas(w,:)=histcounts(abs(C_ens_vec),0:0.05:1,'Normalization','Probability');
        
        

    end
    clear spikes spikes_w dates C_vec C C_ens_vec Ens ense_n ens cells_per_ens sorting_w coac C_ens cells_d

end



%% V1 sessions


V1_ses=find(big_table(:,10)<0); %Sessions inthe in MEC
adults=find(big_table(:,3)>15); %Sessions inthe of adults
sessions=intersect(V1_ses,adults);
count=0;

ws=big_table(sessions,6);
index_nan=isnan(ws);
sessions(index_nan)=[];
session_number=big_table(sessions,4);
days=big_table(sessions,3);
repeated_days=find(diff(days)==0);
clear index_to_delete

index_to_delete=[];
count=0;
for w=1:length(repeated_days)
    N1=big_table(sessions(repeated_days(w)),12);
    N2=big_table(sessions(repeated_days(w)+1),12);
    
    if (N1>N2)
        count=count+1;
        index_to_delete(count)=repeated_days(w)+1;
    elseif  (N1<N2)
        count=count+1;
        index_to_delete(count)=repeated_days(w);    
    else
        print('Other case');
    end   
end

sessions(index_to_delete)=[];
big_table_2=big_table(sessions,:);
V1_sessions=sessions;

%% coac v1


for w=1:length(V1_sessions)
    row_w=V1_sessions(w);
    disp(w)
    
    if num2str(big_table(V1_sessions(w),2))=='1'
        mouse=92227;
    elseif num2str(big_table(V1_sessions(w),2))=='2'
        mouse=92229;
    elseif num2str(big_table(V1_sessions(w),2))=='3'
        mouse=60961;
    end
    
    day=big_table(V1_sessions(w),3);
    s=big_table(V1_sessions(w),4);
    munit=big_table(V1_sessions(w),5);
    
    dt=floor(big_table(V1_sessions(w),8));
    if isinteger(dt)
    else
        dt=floor(dt);
    end
    
    if isnan(dt)
        disp('Error')
        continue
    else
        count=count+1;
        num_clus_discr=10;
        downsample_factor=1;
        load([rec_data_path,strcat('recording_dates_',num2str(mouse),'.mat')]);
        if isfield(dates,'actual_day')==1
            day=find(dates.actual_day==day);
        end
        file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',num2str(mouse),'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        load(file_name_spk,'-mat');
        spikes=full(spikes_d_s);
        
        spikes_w=spikes;
        
        coac=sum(spikes_w)/size(spikes_w,1);
        
        prob_coac_v1(w,:)=histcounts(coac,0:0.01:0.5,'Normalization','Probability');
        
          C=corr(spikes_w');
        C_vec=C(:);
        C_vec(find(C_vec==1))=[];
        
        N=size(spikes_w,1);
        [~,sorting_w,~]=get_sorting(spikes_w);
        
        cells_per_ens=floor(N/10);
        for e=1:10
            ens(e,:)=sorting_w((e-1)*cells_per_ens+1 : e*cells_per_ens);
            ense_n(e,:)=e*ones(1,cells_per_ens);
        end
        Ens(:,1)=reshape(ens,cells_per_ens*10,1);
        Ens(:,2)=reshape(ense_n,cells_per_ens*10,1);
        if (N-10*cells_per_ens)>0
            Ens(10*cells_per_ens+1:N,1)=sorting_w(10*cells_per_ens+1:N);
            Ens(10*cells_per_ens+1:N,2)=ones(length(sorting_w(10*cells_per_ens+1:N)),1)*10;
        end
        
        C_ens_vec=[];
        for i=1:10
            cells_ens=find(Ens(:,2)==i);
            C_ens=corr(spikes_w(Ens(cells_ens,1),:)');
            C_ens_vec=[C_ens_vec;C_ens(:)];
        end
        C_ens_vec(find(C_ens_vec==1))=[];

        prob_C_v1(w,:)=histcounts(abs(C_vec),0:0.05:1,'Normalization','Probability');
        prob_C_end_v1(w,:)=histcounts(abs(C_ens_vec),0:0.05:1,'Normalization','Probability');

    end
    
    clear spikes spikes_w dates C_vec C C_ens_vec Ens ense_n ens cells_per_ens sorting_w coac C_ens cells_d
end

%% MEC wave sessions


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
waves=mec_sessions(find(big_table(mec_sessions,6)==1));
% coac mec


for w=1:length(waves)
    row_w=waves(w);
    disp(w)
        
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);
    munit=big_table(waves(w),5);
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];

    dt=floor(big_table(waves(w),8));
    if isinteger(dt)
    else
        dt=floor(dt);
    end
    
    if isnan(dt)
        disp('Error')
        continue
    else
        count=count+1;
        num_clus_discr=10;
        downsample_factor=1;
        load([rec_data_path,strcat('recording_dates_',num2str(mouse),'.mat')]);
        if isfield(dates,'actual_day')==1
            day=find(dates.actual_day==day);
        end
        file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',num2str(mouse),'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        load(file_name_spk,'-mat');
        spikes=full(spikes_d_s);
        
        spikes_w=spikes;
        
        coac=sum(spikes_w)/size(spikes_w,1);
        
        prob_coac_mec(w,:)=histcounts(coac,0:0.01:0.5,'Normalization','Probability');
        
          C=corr(spikes_w');
        C_vec=C(:);
        C_vec(find(C_vec==1))=[];
        
         N=size(spikes_w,1);
        [~,sorting_w,~]=get_sorting(spikes_w);
        
        cells_per_ens=floor(N/10);
        for e=1:10
            ens(e,:)=sorting_w((e-1)*cells_per_ens+1 : e*cells_per_ens);
            ense_n(e,:)=e*ones(1,cells_per_ens);
        end
        Ens(:,1)=reshape(ens,cells_per_ens*10,1);
        Ens(:,2)=reshape(ense_n,cells_per_ens*10,1);
        if (N-10*cells_per_ens)>0
            Ens(10*cells_per_ens+1:N,1)=sorting_w(10*cells_per_ens+1:N);
            Ens(10*cells_per_ens+1:N,2)=ones(length(sorting_w(10*cells_per_ens+1:N)),1)*10;
        end
        
        C_ens_vec=[];
        for i=1:10
            cells_ens=find(Ens(:,2)==i);
            C_ens=corr(spikes_w(Ens(cells_ens,1),:)');
            C_ens_vec=[C_ens_vec;C_ens(:)];
        end
        C_ens_vec(find(C_ens_vec==1))=[];

        prob_C_end_mec(w,:)=histcounts(abs(C_ens_vec),0:0.05:1,'Normalization','Probability');
        
       
        prob_C_mec(w,:)=histcounts(abs(C_vec),0:0.05:1,'Normalization','Probability');
    end
    
    clear spikes spikes_w dates C_vec C C_ens_vec Ens ense_n ens cells_per_ens sorting_w coac C_ens cells_d
end



%%

mean_c_v1=mean(prob_C_v1);
mean_c_pas=mean(prob_C_pas);
sem_c_v1=std(prob_C_v1)/sqrt(size(prob_C_v1,1));
sem_c_pas=std(prob_C_pas)/sqrt(size(prob_C_pas,1));


ranksum();
% mean_c_mec=mean(prob_C_mec);
% sem_c_mec=std(prob_C_mec)/sqrt(size(prob_C_mec,1));

% figure
% errorbar([0:0.05:1-0.05],mean_c_v1',sem_c_v1,'linewidth',2);
% hold on
% errorbar([0:0.05:1-0.05],mean_c_pas',sem_c_pas,'linewidth',2);
% errorbar([0:0.05:1-0.05],mean_c_mec',sem_c_mec,'linewidth',2);
% legend('V1','PaS','MEC');
% box off
% set(gca,'Yscale','log')
% axis([0 0.5 0.000001 1])
% ylabel('Probability');
% xlabel('|Correlation|');
% set(gca,'fontsize',16)
% legend boxoff


figure
errorbar([0:0.05:1-0.05],mean_c_v1',sem_c_v1,'linewidth',2);
hold on
errorbar([0:0.05:1-0.05],mean_c_pas',sem_c_pas,'linewidth',2);
legend('V1','PaS');
box off
set(gca,'Yscale','log')
axis([0 0.7 0.0000001 1])
ylabel('Probability');
xlabel('|Correlation|');
set(gca,'fontsize',16)
xticks([0 0.3 0.6])
legend boxoff


mean_prob_coac_v1=mean(prob_coac_v1);
mean_prob_coac_pas=mean(prob_coac_pas);
sem_prob_coac_v1=std(prob_coac_v1)/sqrt(size(prob_coac_v1,1));
sem_prob_coac_pas=std(prob_coac_pas)/sqrt(size(prob_coac_pas,1));
mean_prob_coac_mec=mean(prob_coac_mec);
sem_prob_coac_mec=std(prob_coac_mec)/sqrt(size(prob_coac_mec,1));
% 
% figure
% errorbar([0:0.01:0.5-0.01],mean_prob_coac_v1',sem_prob_coac_v1,'linewidth',2);
% hold on
% errorbar([0:0.01:0.5-0.01],mean_prob_coac_pas',sem_prob_coac_pas,'linewidth',2);
% errorbar([0:0.01:0.5-0.01],mean_prob_coac_mec',sem_prob_coac_mec,'linewidth',2);
% legend('V1','PaS','MEC');
% box off
% set(gca,'Yscale','log')
% axis([0 0.2 0.000001 1])
% ylabel('Probability');
% xlabel('Fraction of coactive cells');
% set(gca,'fontsize',16)
% title('Coactivity')
% legend boxoff

figure
errorbar([0:0.01:0.5-0.01],mean_prob_coac_v1',sem_prob_coac_v1,'linewidth',2);
hold on
errorbar([0:0.01:0.5-0.01],mean_prob_coac_pas',sem_prob_coac_pas,'linewidth',2);
legend('V1','PaS');
box off
set(gca,'Yscale','log')
axis([0 0.3 0.000001 1])
ylabel('Probability');
xlabel('Fraction of coactive cells');
set(gca,'fontsize',16)
title('Coactivity')
legend boxoff



% 
% mean_c_end_v1=mean(prob_C_end_v1);
% mean_c_end_pas=mean(prob_C_end_pas);
% sem_c_end_v1=std(prob_C_end_v1)/sqrt(size(prob_C_end_v1,1));
% sem_c_end_pas=std(prob_C_end_pas)/sqrt(size(prob_C_end_pas,1));
% mean_c_end_mec=mean(prob_C_end_mec);
% sem_c_end_mec=std(prob_C_end_mec)/sqrt(size(prob_C_end_mec,1));
% 
% figure
% errorbar([0:0.05:1-0.05],mean_c_end_v1',sem_c_end_v1,'linewidth',2);
% hold on
% errorbar([0:0.05:1-0.05],mean_c_end_pas',sem_c_end_pas,'linewidth',2);
% errorbar([0:0.05:1-0.05],mean_c_end_mec',sem_c_end_mec,'linewidth',2);
% legend('V1','PaS','mec');
% box off
% set(gca,'Yscale','log')
% axis([0 0.5 0.000001 1])
% ylabel('Probability');
% xlabel('|Correlation| in ensembles');
% set(gca,'fontsize',16)
% legend boxoff