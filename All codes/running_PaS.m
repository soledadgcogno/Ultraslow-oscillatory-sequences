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


%------------------------- PaS

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

%% Calculates properties of waves and behavioral variables

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
    
    
%     dt=floor(big_table(V1_sessions(w),8));
%     if isinteger(dt)
%     else
%         dt=floor(dt);
%     end
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

%     % Completed laps per wave
%     for ind_w=1:size(table_u,1)
% %         if ind_w==1
%             completed_laps_per_wave(ind_w + cont_w)=lap_index_d(table_u(ind_w,2)) - lap_index_d(table_u(ind_w,1));
%             distance_run_per_wave(ind_w + cont_w)=tot_position_d(table_u(ind_w,2)) - tot_position_d(table_u(ind_w,1));
%             completed_laps_per_wave_per_session{w}(ind_w)=lap_index_d(table_u(ind_w,2)) - lap_index_d(table_u(ind_w,1));
%             duration_waves{w}(ind_w)=(table_u(ind_w,2)-table_u(ind_w,1))/8;    
%             duration_waves_pooled(ind_w + cont_w)=(table_u(ind_w,2)-table_u(ind_w,1))/8;    
%             if (table_u(ind_w,1)>4*8)
%             speed_pre_post(ind_w + cont_w,1)=mean(speed_d2(table_u(ind_w,1)-4*8:table_u(ind_w,1)));
%             speed_pre_post(ind_w + cont_w,2)=mean(speed_d2(table_u(ind_w,1):table_u(ind_w,1)+4*8));
%             ac_pre_post(ind_w + cont_w,1)=mean(ac_d(table_u(ind_w,1)-4*8:table_u(ind_w,1)));
%             ac_pre_post(ind_w + cont_w,2)=mean(ac_d(table_u(ind_w,1):table_u(ind_w,1)+4*8));
%             end
%     end
%     cont_w=length(completed_laps_per_wave);
    

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
    meadian_speed(w)=median(speed_d2);
    run_distance(w)=tot_position_d(end);
    min_speed(w)=min(speed_d2);
    max_speed(w)=max(speed_d2);
    
    % Information about total number of laps in the session
    number_of_laps(w)=max(lap_index_d);
    
   
        
    clear spikes sorting sorting spikes_d_s table_u table_motion table_still motion
    clear alfa angle_d cells_d lap_index_d lap_position_d speed speed_d speed_d2 spikes spikes_d_osf spikes_d_s subset_waves times_s times_s_d times_tracking ...
        timestam tot_position tot_position_d ac_d motion subset_nowaves subset_waves subset_norunning subset_running wave_vec wave_epoch_speed_per_ses ... 
        wave_epoch_ac_per_ses T T_h T_osf table_u_copy angle angular_speed  lap_index lap_position motor motor_d fraction_still_waves epochs_motion epochs_still
    clear idx_iwi iwi
end
