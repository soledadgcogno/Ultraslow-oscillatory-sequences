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

%% Wave sessions


figpath='C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Raster Plots all MEC sessions\';
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


%% IWI distribution
clus=10;
IWI_ses_pooled=[];
sf=7.73; 

for w=1:length(waves)
    row_w=waves(w);
    disp(w);
    
    count=count+1;
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);
    munit=big_table(waves(w),5);
    
    dt=floor(big_table(waves(w),8));
    if isinteger(dt)
    else
        dt=floor(dt);
    end
    
%     dt_Ent=floor(big_table(waves(w),11));

    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_spk,'-mat');
    spikes=full(spikes_d_s);
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);

    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);        
    
    clear IWI_ses
    for i=2:size(table_u,1)
        IWI_ses(i-1)=(table_u(i,1)-table_u(i-1,2))./sf;
    end
    IWI_ses_pooled=[IWI_ses_pooled,IWI_ses];  

    IWI_per_session{1,w}=IWI_ses;
    
%     figure
%     histogram(IWI_ses,0:2:100)
    
    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2 IWI_ses
    
end

figure
h=histogram(IWI_ses_pooled,[0:20:500]);
figure
bar(h.BinEdges(1:end-1),h.Values);
ylabel("Counts - ICI");
xlabel("ICI (s)");
set(gca,'fontsize',16)
% xticks([0,100,200,300])
axis([-20 500 -0.1 320])
box off
set(gca,'Yscale','log')

min_IWI=min(IWI_ses_pooled);
max_IWI=max(IWI_ses_pooled);

for i=1:15
    delta_IWI{1,i}=diff(IWI_per_session{1,i});
end

length(find(IWI_ses_pooled<0.130));