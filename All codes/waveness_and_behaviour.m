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
    load([save_data_path ['WS_Osc_15_sf7p73II_',mice(m,:),'.mat']]);
%     ws_ent=load([save_data_path ['WS_Entropy_',mice(m,:),'.mat']]);
    ws_ent=load([save_data_path ['WS_Entropy_dt66_',mice(m,:),'.mat']]);
    ws_prob=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_',mice(m,:),'.mat']]);
    ws_prob_sig=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_',mice(m,:),'with_significance.mat']]);

    
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
                    big_table(count,12)=size(spk.spikes_d_s,1); %N
                    big_table(count,13)=WS_stat.fraction(day,s); %fraction - Osc
                    big_table(count,14)=ws_prob.WS_stat.wave_score_prob(day,s); %WS - Prob
                    big_table(count,15)=ws_prob.WS_stat.wave_score_wave(day,s); %WS - Prob - check the wave sessions
                    big_table(count,16)=ws_prob_sig.WS_stat.wave_score_prob_sig(day,s); %WS - Prob



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
    clear WS_stat ws_ent ws_prob
end

 
%% V1

sig_waveness_v1=14;

sig_v1=zeros(1,19);
sig_v1(14)=1;

% V1 sessions


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

% Computes behavioral variables

clus=10;
count=0;

speed_smoothing_kernel=15;
threshold_still=2;
gap_for_motion=0;
gap_for_still=0;

R=8.5145;
ind_still=0;
ind_still_sh=0;
ind_we=0;
cont_w=0;

full_position=[]; %Concatenates position on the wheel across all sessions with waves
full_tot_position=[]; %Concatenates total position across all sessions with waves
full_speed=[]; %Concatenates speed across all sessions with waves
full_ac=[]; %Concatenates ac across all sessions with waves
full_wave_vec=[];   %Concatenates wave_vec across all sessions with waves
number_of_laps_per_wave=[];
duration_waves_pooled=[];
epochs_motion_pooled=[];
epochs_still_pooled=[];

cont=0;
count_iwi=0;
for w=1:length(V1_sessions) %We exclude litter 5 in this analysis
    row_w=V1_sessions(w);
    disp(w)
    
    %Parameters
    count=count+1;
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
    
    load([rec_data_path,strcat('recording_dates_',num2str(mouse),'.mat')]);    
    if isfield(dates,'actual_day')==1
        day=find(dates.actual_day==day);
    end
    
    %Load spike data
    file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',num2str(mouse),'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name,'-mat');
    spikes=full(spikes_d_s);
    [~,T]=size(spikes);
    file_name_osf=[dpath ['spikes_30ms_Do_THR1p5_SNRH_FoV1','_',num2str(mouse),'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_osf,'-mat');
    spikes_d_osf=full(spikes_d_s);
    [N,T_osf]=size(spikes_d_osf);
    
    % Load behavioral data
    if length(num2str(dates.days(day)))==7
        st=num2str(dates.days(day));
        st_new=['0',st];
        file_name_beh=[dbeh_path,num2str(mouse),'\Flavio_2P_V1_Spring 2020_',num2str(mouse),'_',st_new,'_MUnit_',num2str(munit),'_TRACKING.csv'];
    else
        file_name_beh=[dbeh_path,num2str(mouse),'\Flavio_2P_V1_Spring 2020_',num2str(mouse),'_',num2str(dates.days(day)),'_MUnit_',num2str(munit),'_TRACKING.csv'];
    end
    
    alfa= table2array(readtable(file_name_beh)); % table with behavioral information
    timestam=alfa(:,1); %Camera time points
    
    if w>0
        %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
        tot_position=alfa(:,3)/10; %in cm
        lap_position=alfa(:,4)/10; %in cm
        lap_index=alfa(:,5);
        motor=alfa(:,6);
        
        %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
        R=max(lap_position)/(2*pi); %in mm
        speed=diff(tot_position)./diff(timestam); %cm/s % No smoothing
%         ac=diff(speed)./diff(timestam(1:end-1));
        angle= lap_position/R;
        angular_speed=angdiff(angle)./diff(timestam); 
        
        %%%%%%%%%%%%%%%%%%%%%% Imaging time steps and camera time steps
        sampling_rate_tracking=1/(alfa(3,1)-alfa(2,1));
        dt_rate_imaging=1/30.95;
        sampling_rate_imaging_d=1/(dt_rate_imaging*4);
        
        times_s=0.0323:0.0323:(T_osf*0.0323);
        times_s_d=downsample(times_s,4);
        times_tracking=times_s_d(1:T); %Imaging time points
        
        %%%%%%%%%%%%%%%%%%%%%% Interpolated and downsampled quantities
%         angle_d=interp1(timestam,angle,times_tracking);
        tot_position_d=interp1(timestam,tot_position,times_tracking);
        lap_index_d=interp1(timestam,lap_index,times_tracking);
        lap_position_d=interp1(timestam,lap_position,times_tracking);
        motor_d=interp1(timestam,motor,times_tracking);
        
        speed_d=interp1(timestam(1:end-1),speed,times_tracking(1:end-1));
        speed_d=smooth(speed_d,speed_smoothing_kernel);
        speed_d2=speed_d;
        speed_d2(speed_d2<threshold_still)=0;
        speed_d2(end+1)=speed_d2(end);
        ac_d=diff(speed_d2)./diff(times_tracking(1:end));
        ac_d(end+1)=ac_d(end);
        %         ac_d(end+2)=ac_d(end);
        
%     elseif w<=3
%         
%         tot_position=alfa(:,2)/10; %in cm
%         
%         %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
%         
%         speed=diff(tot_position)./diff(timestam); %cm/s % Smoothing kernel of 2.3 seconds
%         
%         %%%%%%%%%%%%%%%%%%%%%% Imaging time steps and camera time steps
%         
%         times_s=0.0323:0.0323:(T_osf*0.0323);
%         times_s_d=downsample(times_s,4);
%         times_tracking=times_s_d(1:T); %Imaging time points
%         
%         %%%%%%%%%%%%%%%%%%%%%% Interpolated and downsampled quantities
%         tot_position_d=interp1(timestam(1:end),tot_position,times_tracking(1:end)); %Interpolated to imagining time points
%         lap_index_d = floor(abs(tot_position_d./(2*pi*R)));
%         lap_index_d(find(lap_index_d==0))=0;
%         lap_position_d=tot_position_d-(2*pi*R.*(lap_index_d));
% %         angle_d=lap_position_d./R;
%         
%         speed_d=interp1(timestam(1:end-1),speed,times_tracking(1:end-1)); %Interpolated to imagining time points
%         speed_d=smooth(speed_d,speed_smoothing_kernel);
%         speed_d2=speed_d;
%         speed_d2(speed_d2<threshold_still)=0;
%         speed_d2(end+1)=speed_d2(end);
%         ac_d=diff(speed_d2)./diff(times_tracking(1:end));
%         ac_d(end+1)=ac_d(end);

    end
    
    full_tot_position=[full_tot_position;tot_position_d];
    full_position=[full_position;lap_position_d];
    full_speed=[full_speed;speed_d2];   %Concatenates speed across all sessions with waves
    full_ac=[full_ac;ac_d];   %Concatenates ac across all sessions with waves
   
    % Idenfity running and stillness epochs
    
    [table_motion,table_still,motion]=get_table_motion_updated_2(speed_d,threshold_still,gap_for_motion,gap_for_still);
    
    epochs_motion=(table_motion(:,2)-table_motion(:,1))/8;
    epochs_still=(table_still(:,2)-table_still(:,1))/8;
    epochs_motion_pooled=[epochs_motion_pooled;epochs_motion];
    epochs_still_pooled=[epochs_still_pooled;epochs_still];

   
    fraction_running(w)= length(find(motion==1))/length(motion);
    fraction_immobility(w)= length(find(motion==0))/length(motion);
%     fraction_wave(w)= length(find(wave_vec==1))/length(wave_vec);
%     fraction_no_wave(w)= length(find(wave_vec==0))/length(wave_vec);
    mean_speed(w)=mean(speed_d2);
    median_speed(w)=median(speed_d2);
    run_distance(w)=tot_position_d(end);
    min_speed(w)=min(speed_d2);
    max_speed(w)=max(speed_d2);
    mean_ac(w)=mean(ac_d);
    median_ac(w)=median(ac_d);
    min_ac(w)=min(ac_d);
    max_ac(w)=max(ac_d);
    
    aux=find(speed_d2>0);
    mean_speed_lzero(w)=mean(speed_d2(aux));
    median_speed_lzero(w)=median(speed_d2(aux));

    % Information about total number of laps in the session
    number_of_laps(w)=max(lap_index_d);
     
   
        
    clear spikes sorting sorting spikes_d_s table_u table_motion table_still motion
    clear alfa angle_d cells_d lap_index_d lap_position_d speed speed_d speed_d2 spikes spikes_d_osf spikes_d_s subset_waves times_s times_s_d times_tracking ...
        timestam tot_position tot_position_d ac_d motion subset_nowaves subset_waves subset_norunning subset_running wave_vec wave_epoch_speed_per_ses ... 
        wave_epoch_ac_per_ses T T_h T_osf table_u_copy angle angular_speed  lap_index lap_position motor motor_d fraction_still_waves epochs_motion epochs_still
    clear idx_iwi iwi aux
end

v1_data(:,1)=big_table(V1_sessions,14); %waveness probability
v1_data(:,2)=big_table(V1_sessions,6); %wavescore
v1_data(:,3)=mean_speed;
v1_data(:,4)=median_speed;
v1_data(:,5)=mean_ac;
v1_data(:,6)=median_ac;
v1_data(:,7)=min_speed;
v1_data(:,8)=max_speed;
v1_data(:,9)=min_ac;
v1_data(:,10)=max_ac;
v1_data(:,11)=run_distance;
v1_data(:,12)=fraction_running;
v1_data(:,13)=mean_speed_lzero;
v1_data(:,14)=median_speed_lzero;

clear mean_speed median_speed mean_ac median_ac min_speed max_speed min_ac max_ac run_distance fraction_running  mean_speed_lzero median_speed_lzero

%% PaS

sig_waveness_pas=[5,6,7,11,12,13,15];
sig_pas=zeros(1,18);
sig_pas(sig_waveness_pas)=1;

%PaS sessions
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

% Behavioural data
clus=10;
count=0;

speed_smoothing_kernel=15;
threshold_still=2;
gap_for_motion=0;
gap_for_still=0;

R=8.5145;
ind_still=0;
ind_still_sh=0;
ind_we=0;
cont_w=0;

full_position=[]; %Concatenates position on the wheel across all sessions with waves
full_tot_position=[]; %Concatenates total position across all sessions with waves
full_speed=[]; %Concatenates speed across all sessions with waves
full_ac=[]; %Concatenates ac across all sessions with waves
full_wave_vec=[];   %Concatenates wave_vec across all sessions with waves
number_of_laps_per_wave=[];
duration_waves_pooled=[];
epochs_motion_pooled=[];
epochs_still_pooled=[];

cont=0;
count_iwi=0;
for w=1:18%length(PaS_sessions) %We exclude litter 5 in this analysis
    row_w=PaS_sessions(w);
    disp(w)
    
    %Parameters
    count=count+1;
    mouse=['L',num2str(big_table(PaS_sessions(w),1)),'M',num2str(big_table(PaS_sessions(w),2))];

    day=big_table(PaS_sessions(w),3);
    s=big_table(PaS_sessions(w),4);
    munit=big_table(PaS_sessions(w),5);
    

    dt=66;
    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    
    %Load spike data
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    
    if isfield(dates,'actual_day')==1
        day=find(dates.actual_day==day);
    end
    
         
    file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name,'-mat');
    spikes=full(spikes_d_s);
    [~,T]=size(spikes);
    file_name_osf=[dpath ['spikes_30ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_osf,'-mat');
    spikes_d_osf=full(spikes_d_s);
    [N,T_osf]=size(spikes_d_osf);
    
  
    % Load behavioral data
    if length(num2str(dates.days(day)))==7
        st=num2str(dates.days(day));
        st_new=['0',st];
        file_name_beh=[dbeh_path,num2str(mouse),'\Flavio_2P_',mouse,'_',st_new,'_MUnit_',num2str(munit),'_TRACKING.csv'];
    else
        file_name_beh=[dbeh_path,num2str(mouse),'\Flavio_2P_',mouse,'_',num2str(dates.days(day)),'_MUnit_',num2str(munit),'_TRACKING.csv'];
    end
    
    alfa= table2array(readtable(file_name_beh)); % table with behavioral information
    timestam=alfa(:,1); %Camera time points
    
    if w>0
% % % %         %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
% % % %         tot_position=alfa(:,3)/10; %in cm
% % % %         lap_position=alfa(:,4)/10; %in cm
% % % %         lap_index=alfa(:,5);
% % % %         motor=alfa(:,6);
% % % %         
% % % %         %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
% % % %         R=max(lap_position)/(2*pi); %in mm
% % % %         speed=diff(tot_position)./diff(timestam); %cm/s % No smoothing
% % % % %         ac=diff(speed)./diff(timestam(1:end-1));
% % % %         angle= lap_position/R;
% % % %         angular_speed=angdiff(angle)./diff(timestam); 
% % % %         
% % % %         %%%%%%%%%%%%%%%%%%%%%% Imaging time steps and camera time steps
% % % %         sampling_rate_tracking=1/(alfa(3,1)-alfa(2,1));
% % % %         dt_rate_imaging=1/30.95;
% % % %         sampling_rate_imaging_d=1/(dt_rate_imaging*4);
% % % %         
% % % %         times_s=0.0323:0.0323:(T_osf*0.0323);
% % % %         times_s_d=downsample(times_s,4);
% % % %         times_tracking=times_s_d(1:T); %Imaging time points
% % % %         
% % % %         %%%%%%%%%%%%%%%%%%%%%% Interpolated and downsampled quantities
% % % % %         angle_d=interp1(timestam,angle,times_tracking);
% % % %         tot_position_d=interp1(timestam,tot_position,times_tracking);
% % % %         lap_index_d=interp1(timestam,lap_index,times_tracking);
% % % %         lap_position_d=interp1(timestam,lap_position,times_tracking);
% % % %         motor_d=interp1(timestam,motor,times_tracking);
% % % %         
% % % %         speed_d=interp1(timestam(1:end-1),speed,times_tracking(1:end-1));
% % % %         speed_d=smooth(speed_d,speed_smoothing_kernel);
% % % %         speed_d2=speed_d;
% % % %         speed_d2(speed_d2<threshold_still)=0;
% % % %         speed_d2(end+1)=speed_d2(end);
% % % %         ac_d=diff(speed_d2)./diff(times_tracking(1:end));
% % % %         ac_d(end+1)=ac_d(end);
% % % %         %         ac_d(end+2)=ac_d(end);
% % % %         
% % % % %     elseif w<=3
% % % % %         
        tot_position=alfa(:,2)/10; %in cm
        
        %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
        
        speed=diff(tot_position)./diff(timestam); %cm/s % Smoothing kernel of 2.3 seconds
        
        %%%%%%%%%%%%%%%%%%%%%% Imaging time steps and camera time steps
        
        times_s=0.0323:0.0323:(T_osf*0.0323);
        times_s_d=downsample(times_s,4);
        times_tracking=times_s_d(1:T); %Imaging time points
        
        %%%%%%%%%%%%%%%%%%%%%% Interpolated and downsampled quantities
        tot_position_d=interp1(timestam(1:end),tot_position,times_tracking(1:end)); %Interpolated to imagining time points
        lap_index_d = floor(abs(tot_position_d./(2*pi*R)));
        lap_index_d(find(lap_index_d==0))=0;
        lap_position_d=tot_position_d-(2*pi*R.*(lap_index_d));
%         angle_d=lap_position_d./R;
        
        speed_d=interp1(timestam(1:end-1),speed,times_tracking(1:end-1)); %Interpolated to imagining time points
        speed_d=smooth(speed_d,speed_smoothing_kernel);
        speed_d2=speed_d;
        speed_d2(speed_d2<threshold_still)=0;
        speed_d2(end+1)=speed_d2(end);
        ac_d=diff(speed_d2)./diff(times_tracking(1:end));
        ac_d(end+1)=ac_d(end);

    end
    
    full_tot_position=[full_tot_position;tot_position_d];
    full_position=[full_position;lap_position_d];
    full_speed=[full_speed;speed_d2];   %Concatenates speed across all sessions with waves
    full_ac=[full_ac;ac_d];   %Concatenates ac across all sessions with waves

    % Idenfity running and stillness epochs
    
    [table_motion,table_still,motion]=get_table_motion_updated_2(speed_d,threshold_still,gap_for_motion,gap_for_still);
    
    epochs_motion=(table_motion(:,2)-table_motion(:,1))/8;
    epochs_still=(table_still(:,2)-table_still(:,1))/8;
    epochs_motion_pooled=[epochs_motion_pooled;epochs_motion];
    epochs_still_pooled=[epochs_still_pooled;epochs_still];

   
    fraction_running(w)= length(find(motion==1))/length(motion);
    fraction_immobility(w)= length(find(motion==0))/length(motion);
%     fraction_wave(w)= length(find(wave_vec==1))/length(wave_vec);
%     fraction_no_wave(w)= length(find(wave_vec==0))/length(wave_vec);
    mean_speed(w)=mean(speed_d2);
    median_speed(w)=median(speed_d2);
    run_distance(w)=tot_position_d(end);
    min_speed(w)=min(speed_d2);
    max_speed(w)=max(speed_d2);
    mean_ac(w)=mean(ac_d);
    median_ac(w)=median(ac_d);
    min_ac(w)=min(ac_d);
    max_ac(w)=max(ac_d);
    
     aux=find(speed_d2>0);
    mean_speed_lzero(w)=mean(speed_d2(aux));
    median_speed_lzero(w)=median(speed_d2(aux));
    
    % Information about total number of laps in the session
    number_of_laps(w)=max(lap_index_d);
   
        
    clear spikes sorting sorting spikes_d_s table_u table_motion table_still motion
    clear alfa angle_d cells_d lap_index_d lap_position_d speed speed_d speed_d2 spikes spikes_d_osf spikes_d_s subset_waves times_s times_s_d times_tracking ...
        timestam tot_position tot_position_d ac_d motion subset_nowaves subset_waves subset_norunning subset_running wave_vec wave_epoch_speed_per_ses ... 
        wave_epoch_ac_per_ses T T_h T_osf table_u_copy angle angular_speed  lap_index lap_position motor motor_d fraction_still_waves epochs_motion epochs_still
    clear idx_iwi iwi aux
end


PaS_data(:,1)=big_table(PaS_sessions(1:18),14); %waveness probability
PaS_data(:,2)=big_table(PaS_sessions(1:18),6); %wavescore
PaS_data(:,3)=mean_speed;
PaS_data(:,4)=median_speed;
PaS_data(:,5)=mean_ac;
PaS_data(:,6)=median_ac;
PaS_data(:,7)=min_speed;
PaS_data(:,8)=max_speed;
PaS_data(:,9)=min_ac;
PaS_data(:,10)=max_ac;
PaS_data(:,11)=run_distance;
PaS_data(:,12)=fraction_running;
PaS_data(:,13)=mean_speed_lzero;
PaS_data(:,14)=median_speed_lzero;


clear mean_speed median_speed mean_ac median_ac min_speed max_speed min_ac max_ac run_distance fraction_running median_speed_lzero mean_speed_lzero

%% MEC

sig_waveness_mec=[3,4,5,6,7,8,9,10,11,12,13,14,15,22,23,24,25,26,27];
sig_mec_aux=[3,4,5,6,7,8,9,10,11,12,13,14]; % Here we discard those sessions for which there are no behavioural data
sig_mec=zeros(1,14);
sig_mec(sig_mec_aux)=1;

%------------------------- MEC


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

% Behavoiural data

clus=10;
count=0;

speed_smoothing_kernel=15;
threshold_still=2;
gap_for_motion=0;
gap_for_still=0;

R=8.5145;
ind_still=0;
ind_still_sh=0;
ind_we=0;
cont_w=0;

full_position=[]; %Concatenates position on the wheel across all sessions with waves
full_tot_position=[]; %Concatenates total position across all sessions with waves
full_speed=[]; %Concatenates speed across all sessions with waves
full_ac=[]; %Concatenates ac across all sessions with waves
full_wave_vec=[];   %Concatenates wave_vec across all sessions with waves
number_of_laps_per_wave=[];
subset_waves_pooled=[];
subset_nowaves_pooled=[];
subset_running_pooled=[];
subset_norunning_pooled=[];
duration_waves_pooled=[];
epochs_motion_pooled=[];
epochs_still_pooled=[];
phase_pooled=[];
session_pooled=[];
cont=0;
count_iwi=0;

for w=1:14 %We exclude litter 5 in this analysis
    row_w=mec_sessions(w);
    disp(w)
    
    %Parameters
    count=count+1;
    mouse=['L',num2str(big_table(mec_sessions(w),1)),'M',num2str(big_table(mec_sessions(w),2))];
    day=big_table(mec_sessions(w),3);
    s=big_table(mec_sessions(w),4);    
    munit=big_table(mec_sessions(w),5); 
    dt=floor(big_table(mec_sessions(w),8));
    if isinteger(dt)
    else
        dt=floor(dt);
    end
    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    
    %Load spike data
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name,'-mat');
    spikes=full(spikes_d_s);
    [~,T]=size(spikes);
    file_name_osf=[dpath ['spikes_30ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_osf,'-mat');
    spikes_d_osf=full(spikes_d_s);
    [N,T_osf]=size(spikes_d_osf);
    
 
        
    % Load behavioral data
    if length(num2str(dates.days(day)))==7
        st=num2str(dates.days(day));
        st_new=['0',st];
        file_name_beh=[dbeh_path,mouse,'\Flavio_2P_',mouse,'_',st_new,'_MUnit_',num2str(munit),'_TRACKING.csv'];
    else
        file_name_beh=[dbeh_path,mouse,'\Flavio_2P_',mouse,'_',num2str(dates.days(day)),'_MUnit_',num2str(munit),'_TRACKING.csv'];
    end
    
    alfa= table2array(readtable(file_name_beh)); % table with behavioral information
    timestam=alfa(:,1); %Camera time points
    
    if w>6 
        %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
        tot_position=alfa(:,3)/10; %in cm
        lap_position=alfa(:,4)/10; %in cm
        lap_index=alfa(:,5);
        motor=alfa(:,6);
        
        %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
        R=max(lap_position)/(2*pi); %in mm
        speed=diff(tot_position)./diff(timestam); %cm/s % No smoothing
%         ac=diff(speed)./diff(timestam(1:end-1));
        angle= lap_position/R;
        angular_speed=angdiff(angle)./diff(timestam); 
        
        %%%%%%%%%%%%%%%%%%%%%% Imaging time steps and camera time steps
        sampling_rate_tracking=1/(alfa(3,1)-alfa(2,1));
        dt_rate_imaging=1/30.95;
        sampling_rate_imaging_d=1/(dt_rate_imaging*4);
        
        times_s=0.0323:0.0323:(T_osf*0.0323);
        times_s_d=downsample(times_s,4);
        times_tracking=times_s_d(1:T); %Imaging time points
        
        %%%%%%%%%%%%%%%%%%%%%% Interpolated and downsampled quantities
%         angle_d=interp1(timestam,angle,times_tracking);
        tot_position_d=interp1(timestam,tot_position,times_tracking);
        lap_index_d=interp1(timestam,lap_index,times_tracking);
        lap_position_d=interp1(timestam,lap_position,times_tracking);
        motor_d=interp1(timestam,motor,times_tracking);
        
        speed_d=interp1(timestam(1:end-1),speed,times_tracking(1:end-1));
        speed_d=smooth(speed_d,speed_smoothing_kernel);
        speed_d2=speed_d;
        speed_d2(speed_d2<threshold_still)=0;
        speed_d2(end+1)=speed_d2(end);
        ac_d=diff(speed_d2)./diff(times_tracking(1:end));
        ac_d(end+1)=ac_d(end);
        %         ac_d(end+2)=ac_d(end);
        
    elseif w<=6
        
        tot_position=alfa(:,2)/10; %in cm
        
        %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
        
        speed=diff(tot_position)./diff(timestam); %cm/s % Smoothing kernel of 2.3 seconds
        
        %%%%%%%%%%%%%%%%%%%%%% Imaging time steps and camera time steps
        
        times_s=0.0323:0.0323:(T_osf*0.0323);
        times_s_d=downsample(times_s,4);
        times_tracking=times_s_d(1:T); %Imaging time points
        
        %%%%%%%%%%%%%%%%%%%%%% Interpolated and downsampled quantities
        tot_position_d=interp1(timestam(1:end),tot_position,times_tracking(1:end)); %Interpolated to imagining time points
        lap_index_d = floor(abs(tot_position_d./(2*pi*R)));
        lap_index_d(find(lap_index_d==0))=0;
        lap_position_d=tot_position_d-(2*pi*R.*(lap_index_d));
%         angle_d=lap_position_d./R;
        
        speed_d=interp1(timestam(1:end-1),speed,times_tracking(1:end-1)); %Interpolated to imagining time points
        speed_d=smooth(speed_d,speed_smoothing_kernel);
        speed_d2=speed_d;
        speed_d2(speed_d2<threshold_still)=0;
        speed_d2(end+1)=speed_d2(end);
        ac_d=diff(speed_d2)./diff(times_tracking(1:end));
        ac_d(end+1)=ac_d(end);

    end
    
    

    % Idenfity running and stillness epochs
    
    [table_motion,table_still,motion]=get_table_motion_updated_2(speed_d,threshold_still,gap_for_motion,gap_for_still);
    
    epochs_motion=(table_motion(:,2)-table_motion(:,1))/8;
    epochs_still=(table_still(:,2)-table_still(:,1))/8;
    epochs_motion_pooled=[epochs_motion_pooled;epochs_motion];
    epochs_still_pooled=[epochs_still_pooled;epochs_still];

    
    cont=length(session_pooled);   
    %Random version of subset_waves

  fraction_running(w)= length(find(motion==1))/length(motion);
    fraction_immobility(w)= length(find(motion==0))/length(motion);
%     fraction_wave(w)= length(find(wave_vec==1))/length(wave_vec);
%     fraction_no_wave(w)= length(find(wave_vec==0))/length(wave_vec);
    mean_speed(w)=mean(speed_d2);
    median_speed(w)=median(speed_d2);
    run_distance(w)=tot_position_d(end);
    min_speed(w)=min(speed_d2);
    max_speed(w)=max(speed_d2);
    mean_ac(w)=mean(ac_d);
    median_ac(w)=median(ac_d);
    min_ac(w)=min(ac_d);
    max_ac(w)=max(ac_d);
    
    
     aux=find(speed_d2>0);
    mean_speed_lzero(w)=mean(speed_d2(aux));
    median_speed_lzero(w)=median(speed_d2(aux));
    
    % Information about total number of laps in the session
    number_of_laps(w)=max(lap_index_d);
    
    clear spikes sorting sorting spikes_d_s table_u table_motion table_still motion
    clear alfa angle_d cells_d lap_index_d lap_position_d speed speed_d speed_d2 spikes spikes_d_osf spikes_d_s subset_waves times_s times_s_d times_tracking ...
        timestam tot_position tot_position_d ac_d motion subset_nowaves subset_waves subset_norunning subset_running wave_vec wave_epoch_speed_per_ses ... 
        wave_epoch_ac_per_ses T T_h T_osf table_u_copy angle angular_speed  lap_index lap_position motor motor_d fraction_still_waves epochs_motion epochs_still
    clear idx_iwi iwi aux
end

mec_data(:,1)=big_table(mec_sessions(1:14),14); %waveness probability
mec_data(:,2)=big_table(mec_sessions(1:14),6); %wavescore
mec_data(:,3)=mean_speed;
mec_data(:,4)=median_speed;
mec_data(:,5)=mean_ac;
mec_data(:,6)=median_ac;
mec_data(:,7)=min_speed;
mec_data(:,8)=max_speed;
mec_data(:,9)=min_ac;
mec_data(:,10)=max_ac;
mec_data(:,11)=run_distance;
mec_data(:,12)=fraction_running;
mec_data(:,13)=mean_speed_lzero;
mec_data(:,14)=median_speed_lzero;

clear mean_speed median_speed mean_ac median_ac min_speed max_speed min_ac max_ac run_distance fraction_running median_speed_lzero mean_speed_lzero

%% Figures
cc=plasma(5);
cc_v1=cc(1,:);
cc_pas=cc(2,:);
cc_mec=cc(3,:);

figure
scatter(mec_data(find(sig_mec),1),mec_data(find(sig_mec),3),60,cc_mec,'filled');
hold on
scatter(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),3),60,cc_pas,'filled');
scatter(v1_data(find(sig_v1),1),v1_data(find(sig_v1),3),60,cc_v1,'filled');
alpha 0.8
xlabel('Significant waveness');
ylabel('Mean speed (cm/s)');
set(gca,'fontsize',16)
legend('MEC','PaS','V1');
legend boxoff

[rho_mec,pval_mec]=corr(mec_data(find(sig_mec),1),mec_data(find(sig_mec),3),'type','spearman');
[rho_pas,pval_pas]=corr(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),3),'type','spearman');
[rho_v1,pval_v1]=corr(v1_data(find(sig_v1),1),v1_data(find(sig_v1),3),'type','spearman');
pool_1=[mec_data(find(sig_mec),1);PaS_data(find(sig_pas),1);v1_data(find(sig_v1),1)];
pool_2=[mec_data(find(sig_mec),3);PaS_data(find(sig_pas),3);v1_data(find(sig_v1),3)];

[rho_tot,pval_tot]=corr(pool_1,pool_2,'type','spearman');



figure
scatter(mec_data(find(sig_mec),1),mec_data(find(sig_mec),4),'filled');
hold on
scatter(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),4),'filled');
scatter(v1_data(find(sig_v1),1),v1_data(find(sig_v1),4),'filled');
alpha 0.8
xlabel('Waveness');
ylabel('Median speed (cm/s)');
set(gca,'fontsize',16)
[rho_mec,pval_mec]=corr(mec_data(find(sig_mec),1),mec_data(find(sig_mec),4),'type','spearman');
[rho_pas,pval_pas]=corr(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),4),'type','spearman');
[rho_v1,pval_v1]=corr(v1_data(find(sig_v1),1),v1_data(find(sig_v1),4),'type','spearman');
pool_1=[mec_data(find(sig_mec),1);PaS_data(find(sig_pas),1);v1_data(find(sig_v1),1)];
pool_2=[mec_data(find(sig_mec),4);PaS_data(find(sig_pas),4);v1_data(find(sig_v1),4)];
[rho_tot,pval_tot]=corr(pool_1,pool_2,'type','spearman');

figure
scatter(mec_data(find(sig_mec),1),mec_data(find(sig_mec),5),'filled');
hold on
scatter(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),5),'filled');
scatter(v1_data(find(sig_v1),1),v1_data(find(sig_v1),5),'filled');
alpha 0.8
xlabel('Waveness');
ylabel('Mean acceleration (cm/s2)');
set(gca,'fontsize',16)
[rho_mec,pval_mec]=corr(mec_data(find(sig_mec),1),mec_data(find(sig_mec),5),'type','spearman');
[rho_pas,pval_pas]=corr(PaS_data(:,1),PaS_data(:,5),'type','spearman');
[rho_v1,pval_v1]=corr(v1_data(find(sig_v1),1),v1_data(find(sig_v1),5),'type','spearman');
pool_1=[mec_data(find(sig_mec),1);PaS_data(find(sig_pas),1);v1_data(find(sig_v1),1)];
pool_2=[mec_data(find(sig_mec),5);PaS_data(find(sig_pas),5);v1_data(find(sig_v1),5)];
[rho_tot,pval_tot]=corr(pool_1,pool_2,'type','spearman');

figure
scatter(mec_data(find(sig_mec),6),mec_data(find(sig_mec),6),'filled');
hold on
scatter(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),6),'filled');
scatter(v1_data(find(sig_v1),1),v1_data(find(sig_v1),6),'filled');
alpha 0.8
ylabel('Waveness');
xlabel('Median acceleration (cm/s2)');
set(gca,'fontsize',16)
[rho_mec,pval_mec]=corr(mec_data(find(sig_mec),1),mec_data(find(sig_mec),6),'type','spearman');
[rho_pas,pval_pas]=corr(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),6),'type','spearman');
[rho_v1,pval_v1]=corr(v1_data(find(sig_v1),1),v1_data(find(sig_v1),6),'type','spearman');
pool_1=[mec_data(find(sig_mec),1);PaS_data(find(sig_pas),1);v1_data(find(sig_v1),1)];
pool_2=[mec_data(find(sig_mec),6);PaS_data(find(sig_pas),6);v1_data(find(sig_v1),6)];
[rho_tot,pval_tot]=corr(pool_1,pool_2,'type','spearman');

figure
scatter(mec_data(find(sig_mec),1),mec_data(find(sig_mec),8),'filled');
hold on
scatter(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),8),'filled');
scatter(v1_data(find(sig_v1),1),v1_data(find(sig_v1),8),'filled');
alpha 0.8
xlabel('Waveness');
ylabel('Maximum speed (cm/s)');
set(gca,'fontsize',16)
[rho_mec,pval_mec]=corr(mec_data(find(sig_mec),1),mec_data(find(sig_mec),8),'type','spearman');
[rho_pas,pval_pas]=corr(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),8),'type','spearman');
[rho_v1,pval_v1]=corr(v1_data(find(sig_v1),1),v1_data(find(sig_v1),8),'type','spearman');
pool_1=[mec_data(find(sig_mec),1);PaS_data(find(sig_pas),1);v1_data(find(sig_v1),1)];
pool_2=[mec_data(find(sig_mec),8);PaS_data(find(sig_pas),8);v1_data(find(sig_v1),8)];
[rho_tot,pval_tot]=corr(pool_1,pool_2,'type','spearman');


figure
scatter(mec_data(find(sig_mec),1),mec_data(find(sig_mec),10)-mec_data(find(sig_mec),9),60,cc_mec,'filled');
hold on
scatter(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),10)-PaS_data(find(sig_pas),9),60,cc_pas,'filled');
scatter(v1_data(find(sig_v1),1),v1_data(find(sig_v1),10)-v1_data(find(sig_v1),9),60,cc_v1,'filled');
alpha 0.8
xlabel('Significant Waveness');
ylabel('Max-Min acceleration (cm/s2)');
set(gca,'fontsize',16);
legend('MEC','PaS','V1');
legend boxoff


[rho_mec,pval_mec]=corr(mec_data(find(sig_mec),1),mec_data(find(sig_mec),10)-mec_data(find(sig_mec),9),'type','spearman');
[rho_pas,pval_pas]=corr(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),10)-PaS_data(find(sig_pas),9),'type','spearman');
[rho_v1,pval_v1]=corr(v1_data(find(sig_v1),1),v1_data(find(sig_v1),10)-v1_data(find(sig_v1),9),'type','spearman');

pool_1=[mec_data(find(sig_mec),1);PaS_data(find(sig_pas),1);v1_data(find(sig_v1),1)];
pool_2=[mec_data(find(sig_mec),10)-mec_data(find(sig_mec),9);PaS_data(find(sig_pas),10)-PaS_data(find(sig_pas),9);v1_data(find(sig_v1),10)-v1_data(find(sig_v1),9)];
[rho_tot,pval_tot]=corr(pool_1,pool_2,'type','spearman');


figure
scatter(mec_data(find(sig_mec),1),mec_data(find(sig_mec),11)/100,60,cc_mec,'filled');
hold on
scatter(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),11)/100,60,cc_pas,'filled');
scatter(v1_data(find(sig_v1),1),v1_data(find(sig_v1),11)/100,60,cc_v1,'filled');
alpha 0.8
xlabel('Significant Waveness');
ylabel('Run distance (m)');
set(gca,'fontsize',16);
legend('MEC','PaS','V1');
legend boxoff


[rho_mec,pval_mec]=corr(mec_data(find(sig_mec),1),mec_data(find(sig_mec),11),'type','spearman')
[rho_pas,pval_pas]=corr(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),11),'type','spearman')
[rho_v1,pval_v1]=corr(v1_data(find(sig_v1),1),v1_data(find(sig_v1),11),'type','spearman')
pool_1=[mec_data(find(sig_mec),1);PaS_data(find(sig_pas),1);v1_data(find(sig_v1),1)];
pool_2=[mec_data(find(sig_mec),11);PaS_data(find(sig_pas),11);v1_data(find(sig_v1),11)];
[rho_tot,pval_tot]=corr(pool_1,pool_2,'type','spearman');

figure
scatter(mec_data(find(sig_mec),1),mec_data(find(sig_mec),12),60,cc_mec,'filled');
hold on
scatter(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),12),60,cc_pas,'filled');
scatter(v1_data(find(sig_v1),1),v1_data(find(sig_v1),12),60,cc_v1,'filled');
alpha 0.8
xlabel('Significant waveness');
ylabel({'Fraction of session with','running behaviour'});
set(gca,'fontsize',16);
legend('MEC','PaS','V1');
legend boxoff

[rho_mec,pval_mec]=corr(mec_data(find(sig_mec),1),mec_data(find(sig_mec),12),'type','spearman');
[rho_pas,pval_pas]=corr(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),12),'type','spearman');
[rho_v1,pval_v1]=corr(v1_data(find(sig_v1),1),v1_data(find(sig_v1),12),'type','spearman');
pool_1=[mec_data(find(sig_mec),1);PaS_data(find(sig_pas),1);v1_data(find(sig_v1),1)];
pool_2=[mec_data(find(sig_mec),12);PaS_data(find(sig_pas),12);v1_data(find(sig_v1),12)];
[rho_tot,pval_tot]=corr(pool_1,pool_2,'type','spearman');


figure
scatter(mec_data(find(sig_mec),1),mec_data(find(sig_mec),13),'filled');
hold on
scatter(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),13),'filled');
scatter(v1_data(find(sig_v1),1),v1_data(find(sig_v1),13),'filled');
alpha 0.8
xlabel('Waveness');
ylabel('Mean speed larger than zero');
set(gca,'fontsize',16);
[rho_mec,pval_mec]=corr(mec_data(find(sig_mec),1),mec_data(find(sig_mec),13),'type','spearman');
[rho_pas,pval_pas]=corr(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),13),'type','spearman');
[rho_v1,pval_v1]=corr(v1_data(find(sig_v1),1),v1_data(find(sig_v1),13),'type','spearman');
pool_1=[mec_data(find(sig_mec),1);PaS_data(find(sig_pas),1);v1_data(find(sig_v1),1)];
pool_2=[mec_data(find(sig_mec),13);PaS_data(find(sig_pas),13);v1_data(find(sig_v1),13)];
[rho_tot,pval_tot]=corr(pool_1,pool_2,'type','spearman');


figure
scatter(mec_data(find(sig_mec),1),mec_data(find(sig_mec),14),'filled');
hold on
scatter(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),14),'filled');
scatter(v1_data(find(sig_v1),1),v1_data(find(sig_v1),14),'filled');
alpha 0.8
xlabel('Waveness');
ylabel('Median speed larger than zero');
set(gca,'fontsize',16);
[rho_mec,pval_mec]=corr(mec_data(find(sig_mec),1),mec_data(find(sig_mec),14),'type','spearman');
[rho_pas,pval_pas]=corr(PaS_data(find(sig_pas),1),PaS_data(find(sig_pas),14),'type','spearman');
[rho_v1,pval_v1]=corr(v1_data(find(sig_v1),1),v1_data(find(sig_v1),14),'type','spearman');
pool_1=[mec_data(find(sig_mec),1);PaS_data(find(sig_pas),1);v1_data(find(sig_v1),1)];
pool_2=[mec_data(find(sig_mec),14);PaS_data(find(sig_pas),14);v1_data(find(sig_v1),14)];
[rho_tot,pval_tot]=corr(pool_1,pool_2,'type','spearman');



%% Waveness and no waveness

clear total_data
total_data=[v1_data;PaS_data;mec_data];
sig=[sig_v1,sig_pas,sig_mec];
brain_area=[ones(1,19),2*ones(1,18),3*ones(1,length(mec_data))];

waveness=find(sig==1);
no_waveness=find(sig<1);

mean_speed_waveness=total_data(waveness,3);
mean_speed_nowaveness=total_data(no_waveness,3);

median_speed_waveness=total_data(waveness,4);
median_speed_nowaveness=total_data(no_waveness,4);

mean_ac_waveness=total_data(waveness,5);
mean_ac_nowaveness=total_data(no_waveness,5);

median_ac_waveness=total_data(waveness,6);
median_ac_nowaveness=total_data(no_waveness,6);

run_d_waveness=total_data(waveness,11);
run_d_nowaveness=total_data(no_waveness,11);

max_ac_waveness=total_data(waveness,10);
max_ac_nowaveness=total_data(no_waveness,10);

min_ac_waveness=total_data(waveness,9);
min_ac_nowaveness=total_data(no_waveness,9);

fraction_m_waveness=total_data(waveness,12);
fraction_m_nowaveness=total_data(no_waveness,12);

mean_speedlzerp_waveness=total_data(waveness,13);
mean_speedlzerp_nowaveness=total_data(no_waveness,13);

median_speedlzerp_waveness=total_data(waveness,14);
median_speedlzerp_nowaveness=total_data(no_waveness,14);

%----- Figures
%1.
% figure
% hold on
% %boxplot(mean_speed)
% [p_mean_speed_wnw]=ranksum(mean_speed_waveness,mean_speed_nowaveness);
% bar([1,2.5],[mean(mean_speed_waveness),mean(mean_speed_nowaveness)]);
% errorbar([1,2.5],[mean(mean_speed_waveness),mean(mean_speed_nowaveness)],...
%     [std(mean_speed_waveness)/sqrt(length(waveness)),std(mean_speed_nowaveness)/sqrt(length(no_waveness))],'k.','linewidth',1.5);
% ylabel('Mean speed (cm/s)');
% h=gca();
% h.XTick=[1 2.8];
% % xticklabels({'Sig. Waveness', 'Non-Sig. Waveness'})
% row1 = {'Significant' 'Non-significant'};
% row2 = {'Waveness' 'Waveness' };
% labelArray = [row1; row2]; 
% tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
% h.XTickLabels = tickLabels; 
% set(gca,'fontsize',16);
% [p,h,stat]=ranksum(mean_speed_waveness,mean_speed_nowaveness);

c_waveness=brain_area(waveness);
c_nowaveness=brain_area(no_waveness);
cc=plasma(5);
figure
bar([1,2],[mean(mean_speed_waveness),mean(mean_speed_nowaveness)],0.2,'FaceColor',[0.75 .75 .75]);
hold on
scatter(1,inf,[],cc(1,:),'filled');
scatter(1,inf,[],cc(2,:),'filled');
scatter(1,inf,[],cc(3,:),'filled');
legend('Mean','V1','PaS','MEC');
h=scatter(0.95*ones(length(waveness),1)+rand(length(waveness),1)/10,mean_speed_waveness,40,cc(c_waveness,:),'filled');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
h=scatter(1.95*ones(length(no_waveness),1)+rand(length(no_waveness),1)/10,mean_speed_nowaveness,40,cc(c_nowaveness,:),'filled');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend boxoff
axis([0 3 0 33]);
box off
ylabel('Mean speed (cm/s)');
h=gca();
h.XTick=[1 2.1];
% xticklabels({'Sig. Waveness', 'Non-Sig. Waveness'})
row1 = {'Significant' 'Non-significant'};
row2 = {'Waveness' 'Waveness' };
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
h.XTickLabels = tickLabels; 
set(gca,'fontsize',16);
legend boxoff
[p,h,stat]=ranksum(mean_speed_waveness,mean_speed_nowaveness);


%2.
% figure
% hold on
% bar([1,2.5],[mean(median_speed_waveness),mean(median_speed_nowaveness)]);
% errorbar([1,2.5],[mean(median_speed_waveness),mean(median_speed_nowaveness)],...
%     [std(median_speed_waveness)/sqrt(length(waveness)),std(median_speed_nowaveness)/sqrt(length(no_waveness))],'k.','linewidth',1.5);
% ylabel('Median speed (cm/s)');
% h=gca();
% h.XTick=[1 2.8];
% row1 = {'Significant' 'Non-significant'};
% row2 = {'Waveness' 'Waveness' };
% labelArray = [row1; row2]; 
% tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
% h.XTickLabels = tickLabels; 
% set(gca,'fontsize',16);
% set(gca,'fontsize',16);
% [p_median_speed_wnw]=ranksum(median_speed_waveness,median_speed_nowaveness);


c_waveness=brain_area(waveness);
c_nowaveness=brain_area(no_waveness);
cc=plasma(5);
figure
bar([1,2],[mean(median_speed_waveness),mean(median_speed_nowaveness)],0.2,'FaceColor',[0.75 .75 .75]);
hold on
scatter(1,inf,[],cc(1,:),'filled');
scatter(1,inf,[],cc(2,:),'filled');
scatter(1,inf,[],cc(3,:),'filled');
legend('Mean','V1','PaS','MEC');
h=scatter(0.95*ones(length(waveness),1)+rand(length(waveness),1)/10,median_speed_waveness,40,cc(c_waveness,:),'filled');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
h=scatter(1.95*ones(length(no_waveness),1)+rand(length(no_waveness),1)/10,median_speed_nowaveness,40,cc(c_nowaveness,:),'filled');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend boxoff
axis([0 3 0 33]);
box off
ylabel('Median speed (cm/s)');
h=gca();
h.XTick=[1 2.1];
% xticklabels({'Sig. Waveness', 'Non-Sig. Waveness'})
row1 = {'Significant' 'Non-significant'};
row2 = {'Waveness' 'Waveness' };
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
h.XTickLabels = tickLabels; 
set(gca,'fontsize',16);
legend boxoff
[p,h,stat]=ranksum(median_speed_waveness,median_speed_nowaveness);



%3.
% figure
% hold on
% bar([1,2.5],[mean(fraction_m_waveness),mean(fraction_m_nowaveness)]);
% errorbar([1,2.5],[mean(fraction_m_waveness),mean(fraction_m_nowaveness)],...
%     [std(fraction_m_waveness)/sqrt(length(waveness)),std(fraction_m_nowaveness)/sqrt(length(no_waveness))],'k.','linewidth',1.5);
% ylabel({'Fraction of the session'; 'with running behaviour'});
% h=gca();
% h.XTick=[1 2.8];
% row1 = {'Significant' 'Non-significant'};
% row2 = {'Waveness' 'Waveness' };
% labelArray = [row1; row2]; 
% tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
% h.XTickLabels = tickLabels; 
% set(gca,'fontsize',16);
% set(gca,'fontsize',16);
[p_median_fractionm_wnw]=ranksum(fraction_m_waveness,fraction_m_nowaveness);

c_waveness=brain_area(waveness);
c_nowaveness=brain_area(no_waveness);
cc=plasma(5);
figure
bar([1,2],[mean(fraction_m_waveness),mean(fraction_m_nowaveness)],0.2,'FaceColor',[0.75 .75 .75]);
hold on
scatter(1,inf,[],cc(1,:),'filled');
scatter(1,inf,[],cc(2,:),'filled');
scatter(1,inf,[],cc(3,:),'filled');
legend('Mean','V1','PaS','MEC');
h=scatter(0.95*ones(length(waveness),1)+rand(length(waveness),1)/10,fraction_m_waveness,40,cc(c_waveness,:),'filled');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
h=scatter(1.95*ones(length(no_waveness),1)+rand(length(no_waveness),1)/10,fraction_m_nowaveness,40,cc(c_nowaveness,:),'filled');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend boxoff
axis([0 3 0 1]);
box off
ylabel({'Fraction of session with';'running behaviour'});
h=gca();
h.XTick=[1 2.1];
% xticklabels({'Sig. Waveness', 'Non-Sig. Waveness'})
row1 = {'Significant' 'Non-significant'};
row2 = {'Waveness' 'Waveness' };
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
h.XTickLabels = tickLabels; 
set(gca,'fontsize',16);
legend boxoff
[p,h,stat]=ranksum(fraction_m_waveness,fraction_m_nowaveness);

%4.
% figure
% hold on
% bar([1,2.5],[mean(mean_speedlzerp_waveness),mean(mean_speedlzerp_nowaveness)]);
% errorbar([1,2.5],[mean(mean_speedlzerp_waveness),mean(mean_speedlzerp_nowaveness)],...
%     [std(mean_speedlzerp_waveness)/sqrt(length(waveness)),std(mean_speedlzerp_nowaveness)/sqrt(length(no_waveness))],'k.','linewidth',1.5);
% ylabel({'Mean speed >0 (cm/s)'});
% h=gca();
% h.XTick=[1 2.8];
% row1 = {'Significant' 'Non-significant'};
% row2 = {'Waveness' 'Waveness' };
% labelArray = [row1; row2]; 
% tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
% h.XTickLabels = tickLabels; 
% set(gca,'fontsize',16);
% set(gca,'fontsize',16);
% [p_median_fractionm_wnw]=ranksum(mean_speedlzerp_waveness,mean_speedlzerp_nowaveness);




%5.
% figure
% hold on
% bar([1,2.5],[mean(median_speedlzerp_waveness),mean(median_speedlzerp_nowaveness)]);
% errorbar([1,2.5],[mean(median_speedlzerp_waveness),mean(median_speedlzerp_nowaveness)],...
%     [std(median_speedlzerp_waveness)/sqrt(length(waveness)),std(median_speedlzerp_nowaveness)/sqrt(length(no_waveness))],'k.','linewidth',1.5);
% ylabel({'Median speed >0 (cm/s)'});
% h=gca();
% h.XTick=[1 2.8];
% row1 = {'Significant' 'Non-significant'};
% row2 = {'Waveness' 'Waveness' };
% labelArray = [row1; row2]; 
% tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
% h.XTickLabels = tickLabels; 
% set(gca,'fontsize',16);
% set(gca,'fontsize',16);
% [p_median_fractionm_wnw]=ranksum(median_speedlzerp_waveness,median_speedlzerp_nowaveness);



%6.
% figure
% hold on
% bar([1,2.5],[mean(run_d_waveness/100),mean(run_d_nowaveness/100)]);
% errorbar([1,2.5],[mean(run_d_waveness/100),mean(run_d_nowaveness/100)],...
%     [std(run_d_waveness/100)/sqrt(length(waveness)),std(run_d_nowaveness/100)/sqrt(length(no_waveness))],'k.','linewidth',1.5);
% ylabel({'Run distance(m)'});
% h=gca();
% h.XTick=[1 2.8];
% row1 = {'Significant' 'Non-significant'};
% row2 = {'Waveness' 'Waveness' };
% labelArray = [row1; row2]; 
% tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
% h.XTickLabels = tickLabels; 
% set(gca,'fontsize',16);
% set(gca,'fontsize',16);
% [p_median_fractionm_wnw]=ranksum(run_d_waveness/100,run_d_nowaveness/100);


c_waveness=brain_area(waveness);
c_nowaveness=brain_area(no_waveness);
cc=plasma(5);
figure
bar([1,2],[mean(run_d_waveness/100),mean(run_d_nowaveness/100)],0.2,'FaceColor',[0.75 .75 .75]);
hold on
scatter(1,inf,[],cc(1,:),'filled');
scatter(1,inf,[],cc(2,:),'filled');
scatter(1,inf,[],cc(3,:),'filled');
legend('Mean','V1','PaS','MEC');
h=scatter(0.95*ones(length(waveness),1)+rand(length(waveness),1)/10,run_d_waveness/100,40,cc(c_waveness,:),'filled');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
h=scatter(1.95*ones(length(no_waveness),1)+rand(length(no_waveness),1)/10,run_d_nowaveness/100,40,cc(c_nowaveness,:),'filled');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend boxoff
axis([0 3 0 1000]);
box off
ylabel({'Run distance (m)'});
h=gca();
h.XTick=[1 2.1];
% xticklabels({'Sig. Waveness', 'Non-Sig. Waveness'})
row1 = {'Significant' 'Non-significant'};
row2 = {'Waveness' 'Waveness' };
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
h.XTickLabels = tickLabels; 
set(gca,'fontsize',16);
legend boxoff
[p,h,stat]=ranksum(run_d_waveness,run_d_nowaveness);



%7.
% figure
% hold on
% bar([1,2.5],[mean(max_ac_waveness),mean(max_ac_nowaveness)]);
% errorbar([1,2.5],[mean(max_ac_waveness),mean(max_ac_nowaveness)],...
%     [std(max_ac_waveness)/sqrt(length(waveness)),std(max_ac_nowaveness)/sqrt(length(no_waveness))],'k.','linewidth',1.5);
% ylabel({'Max Acceleration(cm/s2)'});
% h=gca();
% h.XTick=[1 2.8];
% row1 = {'Significant' 'Non-significant'};
% row2 = {'Waveness' 'Waveness' };
% labelArray = [row1; row2]; 
% tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
% h.XTickLabels = tickLabels; 
% set(gca,'fontsize',16);
% set(gca,'fontsize',16);
% [p_median_fractionm_wnw]=ranksum(max_ac_waveness,max_ac_nowaveness);

%8.
% figure
% hold on
% bar([1,2.5],[mean(max_ac_waveness-min_ac_waveness),mean(max_ac_nowaveness-min_ac_nowaveness)]);
% errorbar([1,2.5],[mean(max_ac_waveness-min_ac_waveness),mean(max_ac_nowaveness-min_ac_nowaveness)],...
%     [std(max_ac_waveness-min_ac_waveness)/sqrt(length(waveness)),std(max_ac_nowaveness-min_ac_nowaveness)/sqrt(length(no_waveness))],'k.','linewidth',1.5);
% ylabel({'Delta Acceleration(cm/s2)'});
% h=gca();
% h.XTick=[1 2.8];
% row1 = {'Significant' 'Non-significant'};
% row2 = {'Waveness' 'Waveness' };
% labelArray = [row1; row2]; 
% tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
% h.XTickLabels = tickLabels; 
% set(gca,'fontsize',16);
% set(gca,'fontsize',16);
% [p_median_fractionm_wnw]=ranksum(max_ac_waveness-min_ac_waveness,max_ac_nowaveness-min_ac_nowaveness);


c_waveness=brain_area(waveness);
c_nowaveness=brain_area(no_waveness);
cc=plasma(5);
figure
bar([1,2],[mean(max_ac_waveness-min_ac_waveness),mean(max_ac_nowaveness-min_ac_nowaveness)],0.2,'FaceColor',[0.75 .75 .75]);
hold on
scatter(1,inf,[],cc(1,:),'filled');
scatter(1,inf,[],cc(2,:),'filled');
scatter(1,inf,[],cc(3,:),'filled');
legend('Mean','V1','PaS','MEC');
h=scatter(0.95*ones(length(waveness),1)+rand(length(waveness),1)/10,max_ac_waveness-min_ac_waveness,40,cc(c_waveness,:),'filled');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
h=scatter(1.95*ones(length(no_waveness),1)+rand(length(no_waveness),1)/10,max_ac_nowaveness-min_ac_nowaveness,40,cc(c_nowaveness,:),'filled');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend boxoff
axis([0 3 0 160]);
box off
ylabel({'Max-Min acceleration (cm/s2)'});
h=gca();
h.XTick=[1 2.1];
% xticklabels({'Sig. Waveness', 'Non-Sig. Waveness'})
row1 = {'Significant' 'Non-significant'};
row2 = {'Waveness' 'Waveness' };
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
h.XTickLabels = tickLabels; 
set(gca,'fontsize',16);
legend boxoff
[p,h,stat]=ranksum(max_ac_waveness-min_ac_waveness,max_ac_nowaveness-min_ac_nowaveness);

%% Waves VS no waves all brain areas
total_data=[v1_data;PaS_data;mec_data];
sig=[sig_v1,sig_pas,sig_mec];
 
 waves=find(total_data(:,2)>8/11);
 no_waves=find(total_data(:,2)<=8/11);
 
 mean_speed_waves=total_data(waves,3);
 mean_speed_nowaves=total_data(no_waves,3);
 
median_speed_waves=total_data(waves,4);
median_speed_nowaves=total_data(no_waves,4);

mean_ac_waves=total_data(waves,5);
mean_ac_nowaves=total_data(no_waves,5);

median_ac_waves=total_data(waves,6);
median_ac_nowaves=total_data(no_waves,6);

run_d_waves=total_data(waves,11);
run_d_nowaves=total_data(no_waves,11);

fraction_m_waves=total_data(waves,12);
fraction_m_nowaves=total_data(no_waves,12);


%----- Figures
%1.
mean_speed=nan(length(no_waves),2);
mean_speed(1:length(mean_speed_waves),1)=mean_speed_waves;
mean_speed(1:length(mean_speed_nowaves),2)=mean_speed_nowaves;

figure
boxplot(mean_speed)
[p_mean_speed_wnw]=ranksum(mean_speed_waves,mean_speed_nowaves);

%2.
median_speed=nan(length(no_waves),2);
median_speed(1:length(median_speed_waves),1)=median_speed_waves;
median_speed(1:length(median_speed_nowaves),2)=median_speed_nowaves;

figure
boxplot(median_speed)
[p_median_speed_wnw]=ranksum(median_speed_waves,median_speed_nowaves);

%3.
median_ac=nan(length(no_waves),2);
median_ac(1:length(median_ac_waves),1)=median_ac_waves;
median_ac(1:length(median_ac_nowaves),2)=median_ac_nowaves;

figure
boxplot(median_ac)
[p_median_ac_wnw]=ranksum(median_ac_waves,median_ac_nowaves);

%4.
median_ac=nan(length(no_waves),2);
median_ac(1:length(median_ac_waves),1)=median_ac_waves;
median_ac(1:length(median_ac_nowaves),2)=median_ac_nowaves;

figure
boxplot(median_ac)
[p_median_ac_wnw]=ranksum(median_ac_waves,median_ac_nowaves);


%% Waves VS no waves MEC

clear total_data
total_data=[mec_data];

waves=find(total_data(:,2)>8/11);
no_waves=find(total_data(:,2)<=8/11);

mean_speed_waves=total_data(waves,3);
mean_speed_nowaves=total_data(no_waves,3);

median_speed_waves=total_data(waves,4);
median_speed_nowaves=total_data(no_waves,4);

mean_ac_waves=total_data(waves,5);
mean_ac_nowaves=total_data(no_waves,5);

median_ac_waves=total_data(waves,6);
median_ac_nowaves=total_data(no_waves,6);

min_ac_waves=total_data(waves,9);
min_ac_nowaves=total_data(no_waves,9);

max_ac_waves=total_data(waves,10);
max_ac_nowaves=total_data(no_waves,10);

run_d_waves=total_data(waves,11);
run_d_nowaves=total_data(no_waves,11);

fraction_m_waves=total_data(waves,12);
fraction_m_nowaves=total_data(no_waves,12);


%----- Figures
%1.
figure
hold on
bar([1,2.5],[mean(mean_speed_waves),mean(mean_speed_nowaves)]);
errorbar([1,2.5],[mean(mean_speed_waves),mean(mean_speed_nowaves)],...
    [std(mean_speed_waves)/sqrt(length(waves)),std(mean_speed_nowaves)/sqrt(length(no_waves))],'k.','linewidth',1.5);
ylabel('Mean speed (cm/s)');
h=gca();
h.XTick=[1 2.8];
% xticklabels({'Sig. Waveness', 'Non-Sig. Waveness'})
row1 = {'Wave' 'No-wave'};
row2 = {'Sessions' 'Sessions' };
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
h.XTickLabels = tickLabels; 
set(gca,'fontsize',16);
[p,h,stat]=ranksum(mean_speed_waves,mean_speed_nowaves);



%2.
figure
hold on
bar([1,2.5],[mean(median_speed_waves),mean(median_speed_nowaves)]);
errorbar([1,2.5],[mean(median_speed_waves),mean(median_speed_nowaves)],...
    [std(median_speed_waves)/sqrt(length(waves)),std(median_speed_nowaves)/sqrt(length(no_waves))],'k.','linewidth',1.5);
ylabel('Median speed (cm/s)');
h=gca();
h.XTick=[1 2.8];
% xticklabels({'Sig. Waveness', 'Non-Sig. Waveness'})
row1 = {'Wave' 'No-wave'};
row2 = {'Sessions' 'Sessions' };
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
h.XTickLabels = tickLabels; 
set(gca,'fontsize',16);
[p,h,stat]=ranksum(median_speed_waves,median_speed_nowaves);

%3.
figure
hold on
bar([1,2.5],[mean(run_d_waves/100),mean(run_d_nowaves/100)]);
errorbar([1,2.5],[mean(run_d_waves/100),mean(run_d_nowaves/100)],...
    [std(run_d_waves/100)/sqrt(length(waves)),std(run_d_nowaves/100)/sqrt(length(no_waves))],'k.','linewidth',1.5);
ylabel('Run distance (m)');
h=gca();
h.XTick=[1 2.8];
% xticklabels({'Sig. Waveness', 'Non-Sig. Waveness'})
row1 = {'Wave' 'No-wave'};
row2 = {'Sessions' 'Sessions' };
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
h.XTickLabels = tickLabels; 
set(gca,'fontsize',16);
[p,h,stat]=ranksum(run_d_waves,run_d_nowaves);

%4.
figure
hold on
bar([1,2.5],[mean(fraction_m_waves),mean(fraction_m_nowaves)]);
errorbar([1,2.5],[mean(fraction_m_waves),mean(fraction_m_nowaves)],...
    [std(fraction_m_waves)/sqrt(length(waves)),std(fraction_m_nowaves)/sqrt(length(no_waves))],'k.','linewidth',1.5);
ylabel({'Fraction of the session'; 'with running behaviour'});
h=gca();
h.XTick=[1 2.8];
row1 = {'Wave' 'No-wave'};
row2 = {'Sessions' 'Sessions' };
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
h.XTickLabels = tickLabels; 
set(gca,'fontsize',16);
[p,h,stat]=ranksum(run_d_waves,run_d_nowaves);

%4.
figure
hold on
bar([1,2.5],[mean(max_ac_waves-min_ac_waves),mean(max_ac_nowaves-min_ac_nowaves)]);
errorbar([1,2.5],[mean(max_ac_waves-min_ac_waves),mean(max_ac_nowaves-min_ac_nowaves)],...
    [std(max_ac_waves-min_ac_waves)/sqrt(length(waves)),std(max_ac_nowaves-min_ac_nowaves)/sqrt(length(no_waves))],'k.','linewidth',1.5);
ylabel({'Delta acceleration (cm/s2)'});
h=gca();
h.XTick=[1 2.8];
row1 = {'Wave' 'No-wave'};
row2 = {'Sessions' 'Sessions' };
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
h.XTickLabels = tickLabels; 
set(gca,'fontsize',16);
[p,h,stat]=ranksum(run_d_waves,run_d_nowaves);



%% Waveness and wavescore for MEC

mec_data_ws(:,1)=big_table(mec_sessions,14); %waveness probability
mec_data_ws(:,2)=big_table(mec_sessions,6); %wavescore
mec_data_ws(:,3)=big_table(mec_sessions,13); %wavescore

figure
hold on
scatter(mec_data_ws(mec_data_ws(:,2)==0,3),mec_data_ws(mec_data_ws(:,2)==0,1),50,'filled');
scatter(mec_data_ws(mec_data_ws(:,2)==1,3),mec_data_ws(mec_data_ws(:,2)==1,1),50,'filled');
alpha 0.5
axis([-0.05 1.05 0.2 0.7])
set(gca,'fontsize',16);
ylabel('Waveness');
xlabel('Wave score');
legend('No-wave sessions','Wave sessions');
legend boxoff