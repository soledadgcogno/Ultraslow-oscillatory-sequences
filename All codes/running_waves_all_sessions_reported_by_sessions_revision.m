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
    ws_prob=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_',mice(m,:),'.mat']]);
 
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
                    big_table(count,11)=ws_prob.WS_stat.wave_score_prob(day,s); %WS - Ent
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

%% Calculates properties of waves and behavioral variables
sf=7.73; 

clus=10;
count=0;

speed_smoothing_kernel=15;
threshold_still=0;
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
full_speed_woTH=[]; %Concatenates speed across all sessions with waves but without applying a threshold
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
count_iwi_2=0;

for w=1:10 %We exclude litter 5 in this analysis
    row_w=waves(w);
    disp(w)
    
    %Parameters
    count=count+1;
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
    animal(w,:)=mouse;
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);    
    munit=big_table(waves(w),5); 
    dt=floor(big_table(waves(w),8));
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
    
    %Calculate phase 
    clear phase_r FRp;
    for i=1:N
        FRp(i,:)=full(fire_rate(spikes(i,:),1*dt,'g')); %smooth in the dt chosen for each session
    end    
    [~,scoret,~] = pca(FRp');
    phase_r=(atan2(scoret(:,2),scoret(:,1)))';
    phase_pooled=[phase_pooled,phase_r];
    
    clear session_i
    session_i=w*ones(1,length(phase_r));
    session_pooled=[session_pooled,session_i];
    

    % Idenfity waves. Generates spike matrix with only wave frames
    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);        
        
    table_u_copy=table_u;
    [~,b]=sort(table_u(:,1));
    table_u=table_u(b,:);   
  
    %Subset of frames with waves (subset_waves) and subset of frames without waves (subset_nowaves)
    subset_waves=[]; %Contains all frames consistent with waves
    for i=1:size(table_u,1)
        subset_waves=[subset_waves,table_u(i,1):table_u(i,2)]; %Frames with waves
    end
    subset_waves_per_session{w}=subset_waves;
    
    subset_nowaves=setxor(subset_waves,1:(T-1)); %Frames with no waves
    
    wave_vec=zeros(1,T);    %Dimension = 1xT - It has a one in the components that indicate a frame consistent with wave-like activity
    wave_vec(subset_waves)=1;
    full_wave_vec=[full_wave_vec,wave_vec];     %Concatenates wave_vec across all sessions with waves
        
    
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
    
    if w>3 
        %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
        tot_position=alfa(:,3)/10; %in cm
        lap_position=alfa(:,4)/10; %in cm
        lap_index=alfa(:,5);
        motor=alfa(:,6);
        
        %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
        R=max(lap_position)/(2*pi); %in cm
        speed=diff(tot_position)./diff(timestam); %cm/s % No smoothing
%         ac=diff(speed)./diff(timestam(1:end-1));
        angle= lap_position/R;
        angular_speed=angdiff(angle)./diff(timestam); 
        
        %%%%%%%%%%%%%%%%%%%%%% Imaging time steps and camera time steps
        sampling_rate_tracking=1/(alfa(3,1)-alfa(2,1)); %Hz
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
        speed_d=smooth(speed_d,speed_smoothing_kernel); %2 seconds of smoothing
        speed_d2=speed_d;
        speed_d2(speed_d2<threshold_still)=0;
        speed_d2(end+1)=speed_d2(end);
%         ac_d=diff(speed_d2)./diff(times_tracking(1:end));
%         ac_d(end+1)=ac_d(end);
        delta_t=diff(times_tracking(1:end));
        ac_d=diff(speed_d2)./delta_t(1);
        ac_d(end)=ac_d(end-1);
        ac_d(end+1)=ac_d(end);

        %         ac_d(end+2)=ac_d(end);
        
    elseif w<=3
        
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
        delta_t=diff(times_tracking(1:end));
        ac_d=diff(speed_d2)./delta_t(1);
        ac_d(end)=ac_d(end-1);
        ac_d(end+1)=ac_d(end);

    end
    
    full_tot_position=[full_tot_position;tot_position_d'];
    full_position=[full_position;lap_position_d'];
    full_speed=[full_speed;speed_d2];   %Concatenates speed across all sessions with waves
    full_speed_woTH=[full_speed_woTH;speed_d];   %Concatenates speed across all sessions with waves
    full_ac=[full_ac;ac_d];   %Concatenates ac across all sessions with waves
    
    min_speed(w)=min(speed_d2);
    max_speed(w)=max(speed_d2);
    min_ac(w)=min(ac_d);
    max_ac(w)=max(ac_d);
    
    % Identify IWI
    clear iwi idx_iwi
    for i=1:size(table_u,1)
        
        if i==1
            iwi(i)=(table_u(i,1))/sf;
            idx_iwi(i)=table_u(i,1); %indicates sequence onset
        else
            iwi(i)=(table_u(i,1)-table_u(i-1,2))/sf;
            idx_iwi(i)=table_u(i,1);
        end
    end

    %Speed and AC IWI for 10 seconds
    width_win_speed=floor(10*sf);% 10 in seconds
    
    for j=1:length(iwi)
        if (idx_iwi(j)-width_win_speed > 0 && idx_iwi(j)+width_win_speed <T)
        count_iwi=count_iwi+1;
        iwi_duration(count_iwi)=iwi(j); % in seconds
        speed_iwi(count_iwi,:)=speed_d2(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed);
        ac_iwi(count_iwi,:)=ac_d(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed);

        track_session_iwi(count_iwi)=w;
        end
%         speed_iwi_m(count_iwi,:)=speed_d2(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed)-mean(speed_d2(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed));
%         ac_iwi_m(count_iwi,:)=ac_d(table_u(idx_iwi(j),1)-8*4:table_u(idx_iwi(j),1)+8*4)-mean(ac_d(table_u(idx_iwi(j),1)-8*4:table_u(idx_iwi(j),1)+8*4));
%         
    end

    %Speed and AC IWI for 2 seconds
    width_win_speed=floor(2*sf);% 2 in seconds

    for j=1:length(iwi)
        if (idx_iwi(j)-width_win_speed > 0 && idx_iwi(j)+width_win_speed <T)
            count_iwi_2=count_iwi_2+1;
%             iwi_duration(count_iwi)=iwi(j); % in seconds
            speed_iwi_2(count_iwi,:)=speed_d2(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed);
            ac_iwi_2(count_iwi,:)=ac_d(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed);

            track_session_iwi_2(count_iwi)=w;
        end
        %         speed_iwi_m(count_iwi,:)=speed_d2(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed)-mean(speed_d2(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed));
        %         ac_iwi_m(count_iwi,:)=ac_d(table_u(idx_iwi(j),1)-8*4:table_u(idx_iwi(j),1)+8*4)-mean(ac_d(table_u(idx_iwi(j),1)-8*4:table_u(idx_iwi(j),1)+8*4));
        %
    end
                
    % Completed laps per wave
    for ind_w=1:size(table_u,1)
        %         if ind_w==1
        completed_laps_per_wave(ind_w + cont_w)=lap_index_d(table_u(ind_w,2)) - lap_index_d(table_u(ind_w,1));
        distance_run_per_wave(ind_w + cont_w)=tot_position_d(table_u(ind_w,2)) - tot_position_d(table_u(ind_w,1));
        completed_laps_per_wave_per_session{w}(ind_w)=lap_index_d(table_u(ind_w,2)) - lap_index_d(table_u(ind_w,1));
        duration_waves{w}(ind_w)=(table_u(ind_w,2)-table_u(ind_w,1))/sf;
        duration_waves_pooled(ind_w + cont_w)=(table_u(ind_w,2)-table_u(ind_w,1))/sf;
        animal_during_wave(ind_w + cont_w,:)=mouse;
        session_during_wave(ind_w + cont_w,:)=s;

        if (table_u(ind_w,1)>4*sf)
            speed_pre_post(ind_w + cont_w,1)=mean(speed_d2(table_u(ind_w,1)-floor(4*sf):table_u(ind_w,1)));
            speed_pre_post(ind_w + cont_w,2)=mean(speed_d2(table_u(ind_w,1):table_u(ind_w,1)+floor(4*sf)));
            ac_pre_post(ind_w + cont_w,1)=mean(ac_d(table_u(ind_w,1)-floor(4*sf):table_u(ind_w,1)));
            ac_pre_post(ind_w + cont_w,2)=mean(ac_d(table_u(ind_w,1):table_u(ind_w,1)+floor(4*sf)));
        end
    end
    cont_w=length(completed_laps_per_wave);


    % Idenfity running and stillness epochs    
    [table_motion,table_still,motion]=get_table_motion_updated_2(speed_d,threshold_still,gap_for_motion,gap_for_still);
    
    epochs_motion=(table_motion(:,2)-table_motion(:,1))/sf;
    epochs_still=(table_still(:,2)-table_still(:,1))/sf;
    epochs_motion_pooled=[epochs_motion_pooled;epochs_motion];
    epochs_still_pooled=[epochs_still_pooled;epochs_still];

    
    % Adds information about immobility epochs to the table "still_epochs".
    % Its collumns are: [Wave-session index ; Duration of still epoch in seconds; Fraction of still epochs for which there were waves]
    for s=1:size(table_still,1)
        ind_still=ind_still+1;
        subset_still=table_still(s,1):table_still(s,2); %Frames with the animal sitting still in still epoch "s"
        fraction_still_waves= length(intersect(subset_still,subset_waves))/length(subset_still); %fraction of the still epoch "s" consistent with waves
        
        still_epochs(ind_still,1)=w; % wave-session index
        still_epochs(ind_still,2)=length(subset_still)./sf; % duration in seconds of still epoch "s"
        still_epochs(ind_still,3)=fraction_still_waves; %fraction of the still epoch "s" consistent with waves
        
        still_epochs_frames(ind_still,1)=cont+table_still(s,1);
        still_epochs_frames(ind_still,2)=cont+table_still(s,2);

        clear subset_still fraction_still_waves
    end
    cont=length(session_pooled);   
    %Random version of subset_waves
    number_frames_with_waves=length(subset_waves);

    fraction_running(w)= length(find(motion==1))/length(motion);
    fraction_immobility(w)= length(find(motion==0))/length(motion);
    fraction_wave(w)= length(find(wave_vec==1))/length(wave_vec);
    fraction_no_wave(w)= length(find(wave_vec==0))/length(wave_vec);
    mean_speed(w)=mean(speed_d2);
    meadisn_speed(w)=median(speed_d2);
    run_distance(w)=tot_position_d(end);
    
    for s=1:size(table_still,1)
        ind_still_sh=ind_still_sh+1;
        subset_still=table_still(s,1):table_still(s,2); %Frames with the animal sitting still in still epoch "s"
        still_epochs_sh(ind_still_sh,1)=w; % wave-session index
        still_epochs_sh(ind_still_sh,2)=length(subset_still)./sf; %  in seconds of still epoch "s"
        
        for sh=1:500
            subset_waves_sh=randperm(T,number_frames_with_waves);
            fraction_still_waves= length(intersect(subset_still,subset_waves_sh))/length(subset_still); %fraction of the still epoch "s" consistent with waves
            
            still_epochs_sh(ind_still_sh,sh+3-1)=fraction_still_waves; %fraction of the still epoch "s" consistent with waves
            
            clear fraction_still_waves subset_waves_sh
        end
        clear subset_still
    end
    
    % Information about total number of laps in the session
    number_of_laps(w)=max(lap_index_d);
    
    % Information about total number of waves in the session
    number_waves(w)=size(table_u,1);
    session_duration(w)=T/sf; %in seconds
    wave_rate(w)=size(table_u,1)/(T/sf);
    
    % Acceleration and speed values in frames consistent with waves. Array for all wave-sessions, with info per wave.
    ind_we=0;
    for we=1:size(table_u,1) %loop on number of waves per session
        ind_we=ind_we+1;
        wave_epoch_speed{ind_we,w}=speed_d2(table_u(we,1):table_u(we,2));
        wave_epoch_ac{ind_we,w}=ac_d(table_u(we,1):table_u(we,2));
    end
    
    %Ratio between (# frames with waves) and (#number of frames with motion). Motion is obtained from get_table_motion_updated 
    locking_to_running(w) = length(subset_waves)./length(find(motion>0));
    
    % Extract information about wave epochs only for this session. One vector per variable (speed and acceleration) 
    wave_epoch_speed_per_ses=[]; %Contains all speed values consistent with waves in this session
    wave_epoch_ac_per_ses=[]; %Contains all acceleration values consistent with waves in this session
    for we=1: number_waves(w) %loop on number of waves per session
        wave_epoch_speed_per_ses=[wave_epoch_speed_per_ses;speed_d2(table_u(we,1):table_u(we,2))];
        wave_epoch_ac_per_ses=[wave_epoch_ac_per_ses;ac_d(table_u(we,1):table_u(we,2))];
    end
    
    %Fraction of frames consistent with waves for wich speed=0
    fraction_speedeqzero(w)=length(find(wave_epoch_speed_per_ses==0))./length(wave_epoch_speed_per_ses);
    %Fraction of frames consistent with waves for wich speed>0
    fraction_speedltzero(w)=length(find(wave_epoch_speed_per_ses>0))./length(wave_epoch_speed_per_ses);
    
    %Difference in behavior Half 1 VS Half 2
    T_h=floor(T/2);
    running_half1(w)=length(find(motion(1:T_h)>0))./T_h;
    running_half2(w)=length(find(motion(T_h+1:T)>0))./T_h;
    
    %Conditional probabilities
    subset_norunning=find(speed_d2==0);
    subset_running=find(speed_d2>0);
    
    p_w_r(w)=length(intersect(subset_waves,subset_running))/length(subset_running);
    p_w_nr(w)=length(intersect(subset_waves,subset_norunning))/length(subset_norunning);
    p_nw_r(w)=length(intersect(subset_nowaves,subset_running))/length(subset_running);
    p_nw_nr(w)=length(intersect(subset_nowaves,subset_norunning))/length(subset_norunning);
    
    p_r_w(w)=length(intersect(subset_waves,subset_running))/length(subset_waves);
    p_nr_w(w)=length(intersect(subset_waves,subset_norunning))/length(subset_waves);
    p_r_nw(w)=length(intersect(subset_nowaves,subset_running))/length(subset_nowaves);
    p_nr_nw(w)=length(intersect(subset_nowaves,subset_norunning))/length(subset_nowaves);
    
    subset_waves_pooled=[subset_waves_pooled,subset_waves];
    subset_nowaves_pooled=[subset_nowaves_pooled,subset_nowaves];
    subset_running_pooled=[subset_running_pooled;subset_running];
    subset_norunning_pooled=[subset_norunning_pooled;subset_norunning];
    
    %     full_speed=[full_speed,speed_d']; %Concatenate speed across sessions
        
    clear spikes sorting sorting spikes_d_s table_u table_motion table_still motion
    clear alfa angle_d cells_d lap_index_d lap_position_d speed speed_d speed_d2 spikes spikes_d_osf spikes_d_s subset_waves times_s times_s_d times_tracking ...
        timestam tot_position tot_position_d ac_d motion subset_nowaves subset_waves subset_norunning subset_running wave_vec wave_epoch_speed_per_ses ... 
        wave_epoch_ac_per_ses T T_h T_osf table_u_copy angle angular_speed  lap_index lap_position motor motor_d fraction_still_waves epochs_motion epochs_still
    clear idx_iwi iwi
    clear ac_d ac_iwi ac_iwi_2 ac_pre_post alfa angle angular_speed b delta_t distance_run_per_wave  FRp   ...
         iwi_dura wave_vec time_s time_s_d times_tracking  still_spochs_frames still_epochs_sh ...
        scoret phase_r session_during_wave session_i speed_iwi speed_iwi_2 still_epochs_frames still_epoch_frames speed_iwi speed_iwi2
end

%% Preparing the data for making the figures

%Speed and Ac during waves
 speed=[]; %concatenates all speed values for all frames with waves in all sessions
ac=[]; %concatenates all acceleration values for all frames with waves in all sessions

for w=1:10 % Speed and Ac during waves
    for i=1: number_waves(w)
        speed = [speed;wave_epoch_speed{i,w}];
        ac=[ac;wave_epoch_ac{i,w}];
    end
end
speed_0=length(find(speed==0))./length(speed); %fraction of frames with speed = 0 in ALL sessions
speed_lt0=length(find(speed>0))./length(speed); %fraction of frames with speed > 0 in ALL sessions

%% Statistics and simple characterization
%Figure epochs motion and epochs speed / Ac and Speed

immobility_bouts=length(epochs_still_pooled);
running_bouts=length(epochs_motion_pooled);

hist_mot=histogram(epochs_motion_pooled,0:2:160);
figure
bar(hist_mot.BinEdges(1:end-1),hist_mot.Values,'FaceColor','k');
set(gca,'fontsize',16)
ylabel('Counts - Running bouts');
xlabel('Time (s)')
set(gca,'YScale', 'log','fontsize',16,'YColor','k','XColor','k');
box off

figure
hist_still=histogram(epochs_still_pooled,0:2:160);
bar(hist_still.BinEdges(1:end-1),hist_still.Values,'FaceColor','k');
set(gca,'fontsize',16)
ylabel('Counts - Immobility bouts');
xlabel('Time (s)')
set(gca,'YScale', 'log','fontsize',16,'YColor','k','XColor','k');
box off

h_baseline=histogram(full_speed,2:2.5:100,'linewidth',0.8,'EdgeColor',[105 105 105]./255,'FaceColor',[255 255 255]./255,'FaceAlpha',0.8);
Val_speed_baseline=h_baseline.Values;
edges_speed_baseline=h_baseline.BinEdges;

figure
colororder({'#009999','#330000'})
fig=bar(edges_speed_baseline(1:end-1),Val_speed_baseline,'linewidth',0.8,'EdgeColor',[105 105 105]./255,'FaceColor',[0 153 153]./255,'FaceAlpha',0.6);
ylabel("Counts")
xlabel('Speed (cm/s)');
set(gca,"fontsize",18);
hold on
xticks([0 50 100])
box off


h_baseline_ac=histogram(full_ac,-50:5:50,'linewidth',0.8,'EdgeColor',[105 105 105]./255,'FaceColor',[255 255 255]./255,'FaceAlpha',0.8);
Val_ac_baseline=h_baseline_ac.Values;
Bins_ac_baseline=h_baseline_ac.BinEdges;

figure
colororder({'#009999','#330000'})
fig=bar(Bins_ac_baseline(1:end-1),Val_ac_baseline,'linewidth',0.8,'EdgeColor',[105 105 105]./255,'FaceColor',[0 153 153]./255,'FaceAlpha',0.6);
ylabel("Counts ")
xlabel('Acceleration (cm/s^2)');
set(gca,"fontsize",18);
box off

figure
histogram(full_speed);
ylabel('Counts');
xlabel('Speed (cm/s)');
set(gca,'fontsize',16)
axis([-5 50 0 inf]);
box off
median_full_speed_WTH=median(full_speed);

figure
histogram(full_speed_woTH);
ylabel('Counts');
xlabel('Speed (cm/s)');
set(gca,'fontsize',16)
axis([-5 50 0 inf]);
box off
median_full_speed_WOTH=median(full_speed_woTH);

%% Figure about immobility epochs with 100% waves
%REMOVE epochs of 0 seconds!!!!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POOL IMMOBILITY EPOCHS AND SPLIT INTO SEVERAL DURATION WINDOWS
edges=[3,5,10,15,20,25];

for w=1:10
    epochs=find(still_epochs(:,1)==w);
    durations = still_epochs(epochs,2);
    fraction_still_waves = still_epochs(epochs,3);
    
    disp(find(durations<0.0001))
    
    % edges=[prctile(durations,33),prctile(durations,66)];
    
    still_epochs_with_waves=find(fraction_still_waves==1); %still epochs that were 100% occupied with waves
    still_epochs_without_waves=find(fraction_still_waves~=1); %still epochs with no waves
    
    durations_with_waves=durations(still_epochs_with_waves); % duration of still epochs that were 100% occupied with waves
    durations_without_waves=durations(still_epochs_without_waves);  %duration of still epochs with no waves
    
    fraction_still_with_waves=fraction_still_waves(still_epochs_with_waves); %fraction of still epochs with waves
    fraction_still_without_waves=fraction_still_waves(still_epochs_without_waves); %fraction of still epochs without 100% of wave occupation
    
    D1=find(durations<edges(1));
    D2=find(durations<edges(2) & durations>edges(1));
    D3=find(durations<edges(3) & durations>edges(2));
    D4=find(durations<edges(4) & durations>edges(3));
    D5=find(durations<edges(5) & durations>edges(4));
    D6=find(durations<edges(6) & durations>edges(5));
    D7=find(durations>edges(6));
    
    D1_waves=find(durations_with_waves<edges(1));
    D2_waves=find(durations_with_waves<edges(2) & durations_with_waves>edges(1));
    D3_waves=find(durations_with_waves<edges(3) & durations_with_waves>edges(2));
    D4_waves=find(durations_with_waves<edges(4) & durations_with_waves>edges(3));
    D5_waves=find(durations_with_waves<edges(5) & durations_with_waves>edges(4));
    D6_waves=find(durations_with_waves<edges(6) & durations_with_waves>edges(5));
    D7_waves=find(durations_with_waves>edges(6));
    
    frac_D1(w)=length(D1_waves)/length(D1);
    frac_D2(w)=length(D2_waves)/length(D2);
    frac_D3(w)=length(D3_waves)/length(D3);
    frac_D4(w)=length(D4_waves)/length(D4);
    frac_D5(w)=length(D5_waves)/length(D5);
    frac_D6(w)=length(D6_waves)/length(D6);
    frac_D7(w)=length(D7_waves)/length(D7);
    
    % Now I do the same for the shuffle realizations
    
    for sh=1:500
        fraction_still_waves_sh = still_epochs_sh(epochs,sh+3-1);
        still_epochs_with_waves_sh=find(fraction_still_waves_sh==1); %still epochs that were 100% occupied with waves
        durations_with_waves_sh=durations(still_epochs_with_waves_sh); % duration of still epochs that were 100% occupied with waves
        
        D1_waves_sh=find(durations_with_waves_sh<edges(1));
        D2_waves_sh=find(durations_with_waves_sh<edges(2) & durations_with_waves_sh>edges(1));
        D3_waves_sh=find(durations_with_waves_sh<edges(3) & durations_with_waves_sh>edges(2));
        D4_waves_sh=find(durations_with_waves_sh<edges(4) & durations_with_waves_sh>edges(3));
        D5_waves_sh=find(durations_with_waves_sh<edges(5) & durations_with_waves_sh>edges(4));
        D6_waves_sh=find(durations_with_waves_sh<edges(6) & durations_with_waves_sh>edges(5));
        D7_waves_sh=find(durations_with_waves_sh>edges(6));
        
        frac_D1_sh(w,sh)=length(D1_waves_sh)/length(D1);
        frac_D2_sh(w,sh)=length(D2_waves_sh)/length(D2);
        frac_D3_sh(w,sh)=length(D3_waves_sh)/length(D3);
        frac_D4_sh(w,sh)=length(D4_waves_sh)/length(D4);
        frac_D5_sh(w,sh)=length(D5_waves_sh)/length(D5);
        frac_D6_sh(w,sh)=length(D6_waves_sh)/length(D6);
        frac_D7_sh(w,sh)=length(D7_waves_sh)/length(D7);
        
%         p_value1=find(frac_D1_sh>D1_waves_sh)/1000;
%         p_value2=find(frac_D2_sh>D2_waves_sh)/1000;
%         p_value3=find(frac_D3_sh>D3_waves_sh)/1000;
%         p_value4=find(frac_D4_sh>D4_waves_sh)/1000;
%         p_value5=find(frac_D5_sh>D5_waves_sh)/1000;
%         p_value6=find(frac_D6_sh>D6_waves_sh)/1000;
%         p_value7=find(frac_D7_sh>D7_waves_sh)/1000;
        
%         mean_frac_D1_sh=mean(frac_D1_sh);
%         mean_frac_D2_sh=mean(frac_D2_sh);
%         mean_frac_D3_sh=mean(frac_D3_sh);
%         mean_frac_D4_sh=mean(frac_D4_sh);
%         mean_frac_D5_sh=mean(frac_D5_sh);
%         mean_frac_D6_sh=mean(frac_D6_sh);
%         mean_frac_D7_sh=mean(frac_D7_sh);
%         
%         sem_frac_D1_sh=std(frac_D1_sh)/sqrt(1000);
%         sem_frac_D2_sh=std(frac_D2_sh)/sqrt(1000);
%         sem_frac_D3_sh=std(frac_D3_sh)/sqrt(1000);
%         sem_frac_D4_sh=std(frac_D4_sh)/sqrt(1000);
%         sem_frac_D5_sh=std(frac_D5_sh)/sqrt(1000);
%         sem_frac_D6_sh=std(frac_D6_sh)/sqrt(1000);
%         sem_frac_D7_sh=std(frac_D7_sh)/sqrt(1000);
       clear  fraction_still_waves_sh still_epochs_with_waves_sh durations_with_waves_sh D1_waves_sh D2_waves_sh D3_waves_sh D4_waves_sh ...
           D5_waves_sh D6_waves_sh D7_waves_sh still_epochs_with_waves_sh durations_with_waves_sh
    end    
    
    clear epochs durations  fraction_still_waves durations_with_waves fraction_still_with_waves still_epochs_without_waves durations_without_waves ...
        fraction_still_without_waves D1 D2 D3 D4 D5 D6 D7  ...
        still_epochs_with_waves D1_waves_sh D2_waves_sh D3_waves_sh D4_waves_sh D5_waves_sh D6_waves_sh D7_waves_sh  ...
        D1_waves D2_waves D3_waves D4_waves D5_waves D6_waves D7_waves 
end

%significance per session - using percentile
figure
for w=1:10
    
    %significance    
    sig_1(w,1)= frac_D1(w)>prctile(frac_D1_sh(w,:),95);
    sig_2(w,1)= frac_D2(w)>prctile(frac_D2_sh(w,:),95);
    sig_3(w,1)= frac_D3(w)>prctile(frac_D3_sh(w,:),95);
    sig_4(w,1)= frac_D4(w)>prctile(frac_D4_sh(w,:),95);
    sig_5(w,1)= frac_D5(w)>prctile(frac_D5_sh(w,:),95);
    sig_6(w,1)= frac_D6(w)>prctile(frac_D6_sh(w,:),95);
    sig_7(w,1)= frac_D7(w)>prctile(frac_D7_sh(w,:),95);
    
    sig_1(w,2)= frac_D1(w)>prctile(frac_D1_sh(w,:),99);
    sig_2(w,2)= frac_D2(w)>prctile(frac_D2_sh(w,:),99);
    sig_3(w,2)= frac_D3(w)>prctile(frac_D3_sh(w,:),99);
    sig_4(w,2)= frac_D4(w)>prctile(frac_D4_sh(w,:),99);
    sig_5(w,2)= frac_D5(w)>prctile(frac_D5_sh(w,:),99);
    sig_6(w,2)= frac_D6(w)>prctile(frac_D6_sh(w,:),99);
    sig_7(w,2)= frac_D7(w)>prctile(frac_D7_sh(w,:),99);
    
    sig_1(w,3)= frac_D1(w)>prctile(frac_D1_sh(w,:),99.9);
    sig_2(w,3)= frac_D2(w)>prctile(frac_D2_sh(w,:),99.9);
    sig_3(w,3)= frac_D3(w)>prctile(frac_D3_sh(w,:),99.9);
    sig_4(w,3)= frac_D4(w)>prctile(frac_D4_sh(w,:),99.9);
    sig_5(w,3)= frac_D5(w)>prctile(frac_D5_sh(w,:),99.9);
    sig_6(w,3)= frac_D6(w)>prctile(frac_D6_sh(w,:),99.9);
    sig_7(w,3)= frac_D7(w)>prctile(frac_D7_sh(w,:),99.9);
    
    mean_frac_D1_sh=mean(frac_D1_sh(w,:));
    mean_frac_D2_sh=mean(frac_D2_sh(w,:));
    mean_frac_D3_sh=mean(frac_D3_sh(w,:));
    mean_frac_D4_sh=mean(frac_D4_sh(w,:));
    mean_frac_D5_sh=mean(frac_D5_sh(w,:));
    mean_frac_D6_sh=mean(frac_D6_sh(w,:));
    mean_frac_D7_sh=mean(frac_D7_sh(w,:));
        
    sem_frac_D1_sh=std(frac_D1_sh(w,:))/sqrt(length(std(frac_D1_sh(w,:))));
    sem_frac_D2_sh=std(frac_D2_sh(w,:))/sqrt(length(std(frac_D2_sh(w,:))));
    sem_frac_D3_sh=std(frac_D3_sh(w,:))/sqrt(length(std(frac_D3_sh(w,:))));
    sem_frac_D4_sh=std(frac_D4_sh(w,:))/sqrt(length(std(frac_D4_sh(w,:))));
    sem_frac_D5_sh=std(frac_D5_sh(w,:))/sqrt(length(std(frac_D5_sh(w,:))));
    sem_frac_D6_sh=std(frac_D6_sh(w,:))/sqrt(length(std(frac_D6_sh(w,:))));
    sem_frac_D7_sh=std(frac_D7_sh(w,:))/sqrt(length(std(frac_D7_sh(w,:))));

    
    subplot(2,5,w)
    plot( [frac_D1(w), frac_D2(w), frac_D3(w), frac_D4(w), frac_D5(w), frac_D6(w), frac_D7(w)],'-*k','Linewidth',2)
    hold on
    errorbar([mean_frac_D1_sh,mean_frac_D2_sh,mean_frac_D3_sh,mean_frac_D4_sh,mean_frac_D5_sh,mean_frac_D6_sh,mean_frac_D7_sh],...
        [sem_frac_D1_sh,sem_frac_D2_sh,sem_frac_D3_sh,sem_frac_D4_sh,sem_frac_D5_sh,sem_frac_D6_sh,sem_frac_D7_sh],'-','Color',[169,169,169]/255,'Linewidth',2)
    axis([0 8 -0.2 1.2])
    xticks([1,2,3,4,5,6,7])
    xticklabels({'0-3','3-5','5-10','10-15','15-20','20-25','>25'})
    xtickangle(45)
    set(gca,'fontsize',12);
    box off
end

%pooling sessions - significance using ttest
[h_t(1),p_t(1),ci,stats_t(1)] =ttest2(frac_D1,frac_D1_sh(:));
[h_t(2),p_t(2),ci,stats_t(2)] =ttest2(frac_D2,frac_D2_sh(:));
[h_t(3),p_t(3),ci,stats_t(3)] =ttest2(frac_D3,frac_D3_sh(:));
[h_t(4),p_t(4),ci,stats_t(4)] =ttest2(frac_D4,frac_D4_sh(:));
[h_t(5),p_t(5),ci,stats_t(5)] =ttest2(frac_D5,frac_D5_sh(:));
[h_t(6),p_t(6),ci,stats_t(6)] =ttest2(frac_D6,frac_D6_sh(:));
[h_t(7),p_t(7),ci,stats_t(7)] =ttest2(frac_D7,frac_D7_sh(:));

[h(1),p(1),stats(1)] =ranksum(frac_D1',frac_D1_sh(:));
[h(2),p(2),stats(2)] =ranksum(frac_D2',frac_D2_sh(:));
[h(3),p(3),stats(3)] =ranksum(frac_D3',frac_D3_sh(:));
[h(4),p(4),stats(4)] =ranksum(frac_D4',frac_D4_sh(:));
[h(5),p(5),stats(5)] =ranksum(frac_D5',frac_D5_sh(:));
[h(6),p(6),stats(6)] =ranksum(frac_D6',frac_D6_sh(:));
[h(7),p(7),stats(7)] =ranksum(frac_D7',frac_D7_sh(:));

%Figure
mean_D1=nanmean(frac_D1);
mean_D2=nanmean(frac_D2);
mean_D3=nanmean(frac_D3);
mean_D4=nanmean(frac_D4);
mean_D5=nanmean(frac_D5);
mean_D6=nanmean(frac_D6);
mean_D7=nanmean(frac_D7);

sem_D1=nanstd(frac_D1)./sqrt(length(~isnan(frac_D1)));
sem_D2=nanstd(frac_D2)./sqrt(length(~isnan(frac_D2)));
sem_D3=nanstd(frac_D3)./sqrt(length(~isnan(frac_D3)));
sem_D4=nanstd(frac_D4)./sqrt(length(~isnan(frac_D4)));
sem_D5=nanstd(frac_D5)./sqrt(length(~isnan(frac_D5)));
sem_D6=nanstd(frac_D6)./sqrt(length(~isnan(frac_D6)));
sem_D7=nanstd(frac_D7)./sqrt(length(~isnan(frac_D7)));


sd_D1=nanstd(frac_D1);
sd_D2=nanstd(frac_D2);
sd_D3=nanstd(frac_D3);
sd_D4=nanstd(frac_D4);
sd_D5=nanstd(frac_D5);
sd_D6=nanstd(frac_D6);
sd_D7=nanstd(frac_D7);

mean_D1_sh=nanmean(frac_D1_sh(:));
mean_D2_sh=nanmean(frac_D2_sh(:));
mean_D3_sh=nanmean(frac_D3_sh(:));
mean_D4_sh=nanmean(frac_D4_sh(:));
mean_D5_sh=nanmean(frac_D5_sh(:));
mean_D6_sh=nanmean(frac_D6_sh(:));
mean_D7_sh=nanmean(frac_D7_sh(:));

sem_D1_sh=nanstd(frac_D1_sh(:))./sqrt(length(~isnan(frac_D1_sh(:))));
sem_D2_sh=nanstd(frac_D2_sh(:))./sqrt(length(~isnan(frac_D2_sh(:))));
sem_D3_sh=nanstd(frac_D3_sh(:))./sqrt(length(~isnan(frac_D3_sh(:))));
sem_D4_sh=nanstd(frac_D4_sh(:))./sqrt(length(~isnan(frac_D4_sh(:))));
sem_D5_sh=nanstd(frac_D5_sh(:))./sqrt(length(~isnan(frac_D5_sh(:))));
sem_D6_sh=nanstd(frac_D6_sh(:))./sqrt(length(~isnan(frac_D6_sh(:))));
sem_D7_sh=nanstd(frac_D7_sh(:))./sqrt(length(~isnan(frac_D7_sh(:))));

sd_D1_sh=nanstd(frac_D1_sh(:));
sd_D2_sh=nanstd(frac_D2_sh(:));
sd_D3_sh=nanstd(frac_D3_sh(:));
sd_D4_sh=nanstd(frac_D4_sh(:));
sd_D5_sh=nanstd(frac_D5_sh(:));
sd_D6_sh=nanstd(frac_D6_sh(:));
sd_D7_sh=nanstd(frac_D7_sh(:));

%Figure with SEM
figure
errorbar([mean_D1,mean_D2,mean_D3,mean_D4,mean_D5,mean_D6,mean_D7],[sem_D1,sem_D2,sem_D3,sem_D4,sem_D5,sem_D6,sem_D7],'k','linewidth',2.5)
hold on
errorbar([mean_D1_sh,mean_D2_sh,mean_D3_sh,mean_D4_sh,mean_D5_sh,mean_D6_sh,mean_D7_sh],...
    [sem_D1_sh,sem_D2_sh,sem_D3_sh,sem_D4_sh,sem_D5_sh,sem_D6_sh,sem_D7_sh],'k','linewidth',2.5,'Color',[169,169,169]/255)
% figure
% errorbar([mean_frac_D1_sh,mean_frac_D2_sh,mean_frac_D3_sh,mean_frac_D4_sh,mean_frac_D5_sh,mean_frac_D6_sh,mean_frac_D7_sh],...
%     [sem_frac_D1_sh,sem_frac_D2_sh,sem_frac_D3_sh,sem_frac_D4_sh,sem_frac_D5_sh,sem_frac_D6_sh,sem_frac_D7_sh],'--*','Color', [170 170 170]/255,'linewidth',2.5,'MarkerSize',10);
% hold on
% plot([frac_D1,frac_D2,frac_D3,frac_D4,frac_D5,frac_D6,frac_D7],'k-*','linewidth',2.5,'MarkerSize',10);
xticks([1,2,3,4,5,6,7])
xticklabels({'0-3','3-5','5-10','10-15','15-20','20-25','>25'})
axis([0.5 7.5 -0.05 1])
box off
ylabel({'Fraction of immobility epochs with ';'oscillations'});
xlabel('Immobility epoch duration (s)')
set(gca,'fontsize',16)

%Figure with SD

figure
errorbar([mean_D1,mean_D2,mean_D3,mean_D4,mean_D5,mean_D6,mean_D7],[sd_D1,sd_D2,sd_D3,sd_D4,sd_D5,sd_D6,sd_D7],'k','linewidth',2.5)
hold on
errorbar([mean_D1_sh,mean_D2_sh,mean_D3_sh,mean_D4_sh,mean_D5_sh,mean_D6_sh,mean_D7_sh],...
    [sd_D1_sh,sd_D2_sh,sd_D3_sh,sd_D4_sh,sd_D5_sh,sd_D6_sh,sd_D7_sh],'k','linewidth',2.5,'Color',[169,169,169]/255)
% figure
% errorbar([mean_frac_D1_sh,mean_frac_D2_sh,mean_frac_D3_sh,mean_frac_D4_sh,mean_frac_D5_sh,mean_frac_D6_sh,mean_frac_D7_sh],...
%     [sem_frac_D1_sh,sem_frac_D2_sh,sem_frac_D3_sh,sem_frac_D4_sh,sem_frac_D5_sh,sem_frac_D6_sh,sem_frac_D7_sh],'--*','Color', [170 170 170]/255,'linewidth',2.5,'MarkerSize',10);
% hold on
% plot([frac_D1,frac_D2,frac_D3,frac_D4,frac_D5,frac_D6,frac_D7],'k-*','linewidth',2.5,'MarkerSize',10);
xticks([1,2,3,4,5,6,7])
xticklabels({'0-3','3-5','5-10','10-15','15-20','20-25','>25'})
axis([0.5 7.5 -0.1 1.2])
box off
ylabel({'Fraction of immobility epochs with ';'oscillations'});
xlabel('Length (s)')
set(gca,'fontsize',16,'yColor','k','xColor','k')

[h_1,p_1,Z1]=ranksum(frac_D1,frac_D1_sh(:));
[h_2,p_2,Z2]=ranksum(frac_D2,frac_D2_sh(:));
[h_3,p_3,Z3]=ranksum(frac_D3,frac_D3_sh(:));
[h_4,p_4,Z4]=ranksum(frac_D4,frac_D4_sh(:));
[h_5,p_5,Z5]=ranksum(frac_D5,frac_D5_sh(:));
[h_5,p_6,Z6]=ranksum(frac_D6,frac_D6_sh(:));
[h_7,p_7,Z7]=ranksum(frac_D7,frac_D7_sh(:));


%% Distribution of fraction of still epoch with waves for
% epochs of more than 25 seconds-
edges=[3,5,10,15,20,25];

figure
for w=1:10
    epochs=find(still_epochs(:,1)==w);
    durations = still_epochs(epochs,2);
    fraction_still_waves = still_epochs(epochs,3);
    
    D7=find(durations>edges(6));
    D7_fraction_waves=fraction_still_waves(D7);
    
    subplot(2,5,w)
    histogram(D7_fraction_waves,[0:0.1:1],'Normalization','Probability');
    box off
    ylabel('Fraction of still epochs');
    xlabel('Fraction of epoch with waves');
    set(gca,'fontsize',16)
    axis([0 1 0 1])
    
    clear epochs durations fraction_still_waves
end

%% Simple statistics

mean_fraction_wave=mean(fraction_wave);
mean_fraction_no_wave=mean(fraction_no_wave);
mean_fraction_running=mean(fraction_running);
mean_fraction_immobility=mean(fraction_immobility);

sem_fraction_wave=std(fraction_wave)/sqrt(10);
sem_fraction_no_wave=std(fraction_no_wave)/sqrt(10);
sem_fraction_running=std(fraction_running)/sqrt(10);
sem_fraction_immobility=std(fraction_immobility)/sqrt(10);

%% Phase of the wave for immobility epochs

cont_epochs=0;
figure
for w=1:10
    epochs=find(still_epochs(:,1)==w);
    durations = still_epochs(epochs,2);
    fraction_still_waves = still_epochs(epochs,3);
    still_epochs_with_waves=find(fraction_still_waves==1); %still epochs that were 100% occupied with waves
    D7=find(durations>edges(6));  %still epochs that longer than a threshold
    D7_still_epochs_with_waves=intersect(still_epochs_with_waves,D7);
    
    for j=1:length(D7_still_epochs_with_waves) %Running over all immobility epochs longer than edges7 and that have a totality of waves 
        
        row_still_epochs_frames=epochs(D7_still_epochs_with_waves(j));
        phase=phase_pooled(still_epochs_frames(row_still_epochs_frames,1):still_epochs_frames(row_still_epochs_frames,2));        
        aux=diff(phase);
        a=length(find(aux>0));
        b=length(find(aux<0));
        
        if b>a
            phase=-phase;
        end

        subplot(22,2,cont_epochs+j)
        plot((1:length(phase))/sf,phase,'k','linewidth',1.5);
        ylabel('');
        yticks([]);
        axis([0 130 -3.14 3.14]);
        ax1 = gca;                   % gca = get current axis
        ax1.YAxis.Visible = 'off';   % remove y-axis        
        box off
        
        %disp(strcat('duration= ',num2str(length(phase)/sf)));
        duration_epochs(cont_epochs+j)=(length(phase)/sf);

%         if cont_epochs+j ~= 42 && cont_epochs+j ~= 43
%           xticks([]);   
%           ax1.XAxis.Visible = 'off';   % remove y-axis   
%         else
%           xticks([0,60,120]);
%           xlabel('Time (s)');
%           set(gca,'fontsize',10)
%         end
        
        clear phase
    end
    
    if isempty(D7_still_epochs_with_waves)==1
        disp(cont_epochs)
    else
        cont_epochs=cont_epochs+j;
        disp(cont_epochs)
    end

    
    clear epochs durations fraction_still_waves  epochs durations fraction_still_waves still_epochs_with_waves D7 ...
        D7_still_epochs_with_waves aux
end


% FM1=fraction_still_waves(D1); 
% FM2=fraction_still_waves(D2);
% FM3=fraction_still_waves(D3);
% FM4=fraction_still_waves(D4);
% FM5=fraction_still_waves(D5);
% FM6=fraction_still_waves(D6);
% FM7=fraction_still_waves(D7);
% 
% 
% FM1_f=FM1(find(FM1>0));
% FM2_f=FM2(find(FM2>0));
% FM3_f=FM3(find(FM3>0));
% FM4_f=FM3(find(FM4>0));
% FM5_f=FM3(find(FM5>0));
% FM6_f=FM3(find(FM6>0));
% FM7_f=FM3(find(FM7>0));

% FM1_Thr=length(find(FM1_f>0.9));
% FM2_Thr=length(find(FM2_f>0.9));
% FM3_Thr=length(find(FM3_f>0.9));
% FM4_Thr=length(find(FM4_f>0.9));
% FM5_Thr=length(find(FM5_f>0.9));
% FM6_Thr=length(find(FM6_f>0.9));
% FM7_Thr=length(find(FM7_f>0.9));

% FM1_Thr=length(find(FM1_f>0.9));
% FM2_Thr=length(find(FM2_f>0.9));
% FM3_Thr=length(find(FM3_f>0.9));
% FM4_Thr=length(find(FM4_f>0.9));
% FM5_Thr=length(find(FM5_f>0.9));
% FM6_Thr=length(find(FM6_f>0.9));
% FM7_Thr=length(find(FM7_f>0.9));

% FM1_Thr=length(FM1);
% FM2_Thr=length(FM2);
% FM3_Thr=length(FM3);
% FM4_Thr=length(FM4);
% FM5_Thr=length(FM5);
% FM6_Thr=length(FM6);
% FM7_Thr=length(FM7);

% figure
% plot([FM1_Thr,FM2_Thr,FM3_Thr,FM4_Thr,FM5_Thr,FM6_Thr,FM7_Thr],'k-*','linewidth',2.5,'MarkerSize',10);
% xticks([1,2,3,4,5,6,7])
% xticklabels({'0-3','3-5','5-10','10-15','15-20','20-25','>25'})
% box off
% ylabel('Immobility Epochs with > 90% of waves #');
% xlabel('Immobility epoch duration (s)')
% set(gca,'fontsize',16)

% figure
% plot([FM1_Thr,FM2_Thr,FM3_Thr,FM4_Thr,FM5_Thr,FM6_Thr,FM7_Thr],'k-*','linewidth',2.5,'MarkerSize',10);
% xticks([1,2,3,4,5,6,7])
% xticklabels({'0-3','3-5','5-10','10-15','15-20','20-25','>25'})
% box off
% ylabel('Immobility epochs with waves #');
% xlabel('Immobility epoch duration (s)')
% set(gca,'fontsize',16)
% 
% figure
% plot([FM1_Thr,FM2_Thr,FM3_Thr,FM4_Thr,FM5_Thr,FM6_Thr,FM7_Thr],'k-*','linewidth',2.5,'MarkerSize',10);
% xticks([1,2,3,4,5,6,7])
% xticklabels({'0-3','3-5','5-10','10-15','15-20','20-25','>25'})
% box off
% ylabel('Immobility epochs with waves #');
% xlabel('Immobility epoch duration (s)')
% set(gca,'fontsize',16)
% set(gca, 'YScale', 'log')
% figure
% subplot(3,1,1)
% % histogram(FM1_f,0:0.05:1,'Normalization','Probability');
% histogram(FM1_f,0:0.05:1,'FaceColor',[0    0.5000    0.4000]);
% ylabel('Counts')
% xticks([])
% title(['0 - ',num2str(edges(1)),' s'])
% % xlabel('Fraction of stillnes period with waves')
% axis([0 inf 0 inf])
% box off
% set(gca,'fontsize',23)
% subplot(3,1,2)
% % histogram(FM2_f,0:0.05:1,'Normalization','Probability');
% histogram(FM2_f,0:0.05:1,'FaceColor',[0    0.5000    0.4000]);
% ylabel('Counts')
% xticks([])
% % xlabel('Fraction of stillnes period with waves')
% title([num2str(edges(1)),' - ',num2str(edges(2)),' s'])
% axis([0 inf 0 inf])
% box off
% set(gca,'fontsize',23)
% subplot(3,1,3)
% % histogram(FM3_f,0:0.05:1,'Normalization','Probability');
% histogram(FM3_f,0:0.05:1,'FaceColor',[0    0.5000    0.4000]);
% xticks([0 0.5 1])
% xticklabels({0 50 100})
% ylabel('Counts')
% xlabel('Immobility epochs with waves %')
% title(['> ',num2str(edges(2)),' s'])
% axis([0 inf 0 inf])
% box off
% set(gca,'fontsize',23)


%%   Fraction of frames with waves and speed>0 or speed=0  - Pie charts
figure
labels = {'Speed=0','Speed>0'};
pie([speed_0,speed_lt0],labels)
colormap summer
set(gca,'fontsize',16)

figure
for w=1:10%length(waves)
 subplot(5,2,w)
labels = {'Speed=0','Speed>0'};
pie([fraction_speedeqzero(w),fraction_speedltzero(w)],labels)
colormap summer
set(gca,'fontsize',16)    
end
 

%%  Histogram of speed values>0 and histogram of acceleration - Session by session

for w=1:10
    disp(w)
    frames=find(session_pooled==w);
    wave_v=full_wave_vec(frames);
    phase_v=phase_pooled(frames);
    speed_v=full_speed(frames);
    ac_v=full_ac(frames);
    
    speed_positive=find(speed_v>2);
    speed_v_f=speed_v(speed_positive);
    ac_v_f=ac_v(speed_positive);
    phase_v_f=phase_v(speed_positive);
    wave_v_f=wave_v(speed_positive);    
    
    full_speed_wave_f=speed_v_f(find(wave_v_f>0));
    full_speed_no_wave_f=speed_v_f(find(wave_v_f==0));
  
    exp_eqn='a*exp(-b*x)';
    if (isempty(full_speed_no_wave_f)~=1)
        max_speed=max(speed_v_f);
        histo_speed_wave=histogram(full_speed_wave_f,[0:max_speed/10:max_speed],'Normalization','Probability');
        prob_speed_wave=histo_speed_wave.Values;
        histo_speed_no_wave=histogram(full_speed_no_wave_f,[0:max_speed/10:max_speed],'Normalization','Probability');
        prob_speed_no_wave=histo_speed_no_wave.Values;
%         kl(w)=kldiv(prob_speed_no_wave,prob_speed_wave);
        
        x=[1:10]';
        fit_wave=fit(x,prob_speed_wave',exp_eqn,'Start',[0.1,0.1]);
        fit_no_wave=fit(x,prob_speed_no_wave',exp_eqn,'Start',[0.1,0.1]);
        kl(w)=kldiv((fit_wave.a*exp(-fit_wave.b*x)')/sum(fit_wave.a*exp(-fit_wave.b*x)'),(fit_no_wave.a*exp(-fit_no_wave.b*x)')/sum(fit_no_wave.a*exp(-fit_no_wave.b*x)'));

        for sh=1:100
            wave_v_sh=shuffle(wave_v_f);
            full_speed_wave_sh=speed_v_f(wave_v_sh>0);
            full_speed_no_wave_sh=speed_v_f(wave_v_sh==0);
            
            histo_speed_wave_sh=histogram(full_speed_wave_sh,[0:max_speed/10:max_speed],'Normalization','Probability');
            prob_speed_wave_sh=histo_speed_wave_sh.Values;
            histo_speed_no_wave_sh=histogram(full_speed_no_wave_sh,[0:max_speed/10:max_speed],'Normalization','Probability');
            prob_speed_no_wave_sh=histo_speed_no_wave_sh.Values;
            
            x=[1:10]';
            fit_wave_sh=fit(x,prob_speed_wave_sh',exp_eqn,'Start',[0.1,0.1]);
            fit_no_wave_sh=fit(x,prob_speed_no_wave_sh',exp_eqn,'Start',[0.1,0.1]);
            
            kl_sh(w,sh)=kldiv((fit_wave_sh.a*exp(-fit_wave_sh.b*x)')/sum(fit_wave_sh.a*exp(-fit_wave_sh.b*x)'),...
                (fit_no_wave_sh.a*exp(-fit_no_wave_sh.b*x)')/sum(fit_no_wave_sh.a*exp(-fit_no_wave_sh.b*x)'));
            
            clear full_speed_wave_sh full_speed_no_wave_sh histo_speed_wave_sh prob_speed_wave_sh histo_speed_no_wave_sh prob_speed_no_wave_sh
        end

        
    else
       kl(w)=NaN;
       kl_sh(w,1:1000)=NaN;
       sig_kl(w,1:3)=-10;         
    end
    
 clear frames wave_v phase_v speed_v ac_v speed_positive speed_v_f ac_v_f phase_v_f wave_v_f full_speed_wave_f full_speed_no_wave_f max_speed histo_speed_wave ...
     prob_speed_wave histo_speed_no_wave prob_speed_no_wave wave_v_sh 
end

for w=1:10      
    sig_kl(w,1)=kl(w)>prctile(kl_sh(w,:),95);
    sig_kl(w,2)=kl(w)>prctile(kl_sh(w,:),99);
    sig_kl(w,3)=kl(w)>prctile(kl_sh(w,:),99.9);
    %         sig_kl(w,1:3)=-10;    
end


figure
bar(nanmean(kl),nanmean(kl_sh(:)))
[hkl,p_kl]=ttest2(kl,kl_sh(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Histogram of speed values>0 and histogram of acceleration    
figure
y=histogram(speed,2:2.5:100);
% bar(y.BinEdges(1:end-1),y.Values,'FaceColor',[192 192 192]./255)
bar(y.BinEdges(1:end-1),y.Values,'FaceColor',[0.97,0.97,0.49],'EdgeColor',[0.68,0.65,0.01],'LineWidth',1.5)
alpha 0.7
axis([0 60 0 inf])
box off    
ylabel('Counts');
xlabel('Speed (cm/s)');
yticks([0 5000])
set(gca,'fontsize',26)

h_baseline=histogram(full_speed,2:2.5:100,'linewidth',0.8,'EdgeColor',[105 105 105]./255,'FaceColor',[255 255 255]./255,'FaceAlpha',0.8);
Val_speed_baseline=h_baseline.Values;
edges_speed_baseline=h_baseline.BinEdges;
h_wave=histogram(speed,2:2.5:100,'linewidth',0.8,'EdgeColor',[105 105 105]./255,'FaceColor',[204 255 204]./255,'FaceAlpha',0.6);
Val_speed_wave=h_wave.Values;

figure
colororder({'#009999','#330000'})
yyaxis left
fig=bar(h_wave.BinEdges(1:end-1),Val_speed_wave,'linewidth',0.8,'EdgeColor',[105 105 105]./255,'FaceColor',[0 153 153]./255,'FaceAlpha',0.6);
ylabel("Counts - Speed|Sequences")
xlabel("Speed (cm/s)");
set(gca,"fontsize",16,'YColor','k','XColor','k');
hold on
yyaxis right
alpha 0.6
fig=plot(h_wave.BinEdges(1:end-1),Val_speed_baseline,'k','linewidth',2.5);
ylabel("Counts - Speed all values")
box off

figure
y=histogram(ac,-50:2.5:50,'FaceColor',[0 153 153]./255);
% bar(y.BinEdges(1:end-1),y.Values,'FaceColor',[0 153 153]./255)
bar(y.BinEdges(1:end-1),y.Values,'FaceColor',[0.97,0.97,0.49],'EdgeColor',[0.68,0.65,0.01],'LineWidth',1.5)
alpha 0.7
axis([-35 35 0 inf])
box off    
ylabel('Counts');
xlabel('Acceleration (cm/s^2)');
yticks([0 40000])
set(gca,'fontsize',26)

h_wave_ac=histogram(ac,-50:2.5:50,'linewidth',0.8,'EdgeColor',[105 105 105]./255,'FaceColor',[204 255 204]./255,'FaceAlpha',0.6);
Val_ac_wave=h_wave_ac.Values;
h_edges=h_wave_ac.BinEdges(1:end-1);
h_baseline_ac=histogram(full_ac,-50:2.5:50,'linewidth',0.8,'EdgeColor',[105 105 105]./255,'FaceColor',[255 255 255]./255,'FaceAlpha',0.8);
Val_ac_baseline=h_baseline_ac.Values;


figure
colororder({'#009999','#330000'})
yyaxis left
fig=bar(h_edges,Val_ac_wave,'linewidth',0.8,'EdgeColor',[105 105 105]./255,'FaceColor',[0 153 153]./255,'FaceAlpha',0.6);
ylabel("Counts - Acceleration|Sequences");
xlabel('Acceleration (cm/s^2)');
set(gca,"fontsize",16,'YColor','k','XColor','k');
hold on
yyaxis right
alpha 0.6
fig=plot(h_edges,Val_ac_baseline,'k','linewidth',2.5);
ylabel("Counts - Acceleration all values")
xticks([-40,-20,0,20,40]);
xlim([-40 40])
box off



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Number of laps per wave and as a function of duration (1 session)
ses=2;
figure
scatter(1:(number_waves(ses)),completed_laps_per_wave_per_session{ses},80,[0.8500, 0.3250, 0.0980],'filled');
xlabel('Wave #');
ylabel('Completed laps #');
set(gca,'fontsize',18)
axis([0 50 -0.5 3.5])
yticks([0 3])
xticks([1 50])

figure
scatter(duration_waves{ses},completed_laps_per_wave_per_session{ses},80,[0.8500, 0.3250, 0.0980],'filled');
xlabel('Wave duration (s)');
ylabel('Completed laps #');
set(gca,'fontsize',18)
axis([10 60 -0.5 3.5])
yticks([0 3])
xticks([10 35 60])



%% Number of waves as a function of number of completed laps


figure
scatter(duration_waves_pooled,completed_laps_per_wave,80,[0.8500, 0.3250, 0.0980],'filled');
alpha 0.7 
xlabel('Wave duration (s)');
ylabel('Completed laps #');
set(gca,'fontsize',18)
axis([0 220 -5 85])
[p_corr,h_corr]=corr(duration_waves_pooled',completed_laps_per_wave');

figure
scatter(duration_waves_pooled,distance_run_per_wave,80,[0.8500, 0.3250, 0.0980],'filled');
alpha 0.7 
xlabel('Wave duration (s)');
ylabel('Run distance (cm)');
set(gca,'fontsize',18)
axis([0 220 -5 5000])
[p_corr,h_corr]=corr(duration_waves_pooled',completed_laps_per_wave');

figure
scatter(number_waves,number_of_laps,100,[102 102 102]./255,'filled');
xlabel('Waves #');
ylabel('Laps #');
set(gca,'fontsize',18)
axis([0 70 0 1200])
yticks([0 600 1200])

for i=1:size(animal_during_wave,1)
    if animal_during_wave(i,:)=='L8M2'
        animal_wave(i)=    1;
    elseif animal_during_wave(i,:)=='L9M1'
        animal_wave(i)=    2;
    elseif animal_during_wave(i,:)=='L9M4'
        animal_wave(i)=    3;
    end
end


aux1=find(animal_wave==1);
aux2=find(animal_wave==2);
aux3=find(animal_wave==3);

cc=hsv(14);
color_mouse(1,:)=cc(9,:);
color_mouse(2,:)=cc(12,:);
color_mouse(3,:)=cc(14,:); 
 
figure
scatter(1:length(aux1),completed_laps_per_wave(aux1),80,color_mouse(1,:),'filled');
hold on
scatter(length(aux1)+1:length(aux1)+length(aux2),completed_laps_per_wave(length(aux1)+1:length(aux1)+length(aux2)),80,color_mouse(2,:),'filled');
scatter(length(aux1)+length(aux2)+1:length(aux1)+length(aux2)+length(aux3),completed_laps_per_wave(length(aux1)+length(aux2)+1:length(aux1)+length(aux2)+length(aux3)),80,color_mouse(3,:),'filled');
alpha 0.7 
xlabel('Wave #');
ylabel('Completed laps');
set(gca,'fontsize',18)
axis([0 350 0 90])
legend('Animal 1','Animal 2','Animal 3')
legend boxoff


median_animal1=median(completed_laps_per_wave(aux1));
meadian_animal2=median(completed_laps_per_wave(length(aux1)+1:length(aux1)+length(aux2)));
meadian_animal3=median(completed_laps_per_wave(length(aux1)+length(aux2)+1:length(aux1)+length(aux2)+length(aux3)));

figure
scatter(1:length(completed_laps_per_wave),completed_laps_per_wave,80,[0.8500, 0.3250, 0.0980],'filled');
alpha 0.7 
xlabel('Wave #');
ylabel('Completed laps');
set(gca,'fontsize',18)
axis([0 350 0 90])



% One figure per animal

line1=29.5;
line2=29+46+0.5;
figure
scatter(1:length(aux1),completed_laps_per_wave(aux1),80,'k','filled');
alpha 0.7 
xlabel('Cycle number');
ylabel('Completed laps');
set(gca,'fontsize',18,'ycolor','k','xcolor','k')
axis([0 105 0 6])
hold on
xline(line1,'r--','linewidth',2.5)
xline(line2,'r--','linewidth',2.5)
title('L8M2')

line1=28.5;
line2=28+76+0.5;
figure
scatter(1:length(aux2),completed_laps_per_wave(length(aux1)+1:length(aux1)+length(aux2)),80,'k','filled');
alpha 0.7 
xlabel('Cycle number');
ylabel('Completed laps');
set(gca,'fontsize',18)
axis([0 115 0 40])
hold on
xline(line1,'r--','linewidth',2.5)
xline(line2,'r--','linewidth',2.5)
title('L9M1');

line1=13.5;
line2=13+24+0.5;
line3=13+24+39+0.5;
figure
scatter(1:length(aux3),completed_laps_per_wave(length(aux1)+length(aux2)+1:length(aux1)+length(aux2)+length(aux3)),80,'k','filled');
alpha 0.7 
xlabel('Cycle number');
ylabel('Completed laps');
set(gca,'fontsize',18,'ycolor','k','xcolor','k');
axis([0 100 0 90])
hold on
xline(line1,'r--','linewidth',2.5)
xline(line2,'r--','linewidth',2.5)
xline(line3,'r--','linewidth',2.5)
title('L9M4');

%% Speed per IWI - 10 s

bins_number_epoch=floor(size(speed_iwi,2)/2);

for i=1:size(iwi_duration,2)
   speed_before(i)=mean(speed_iwi(i,1:bins_number_epoch));
   speed_after(i)=mean(speed_iwi(i,bins_number_epoch+1:end));
   delta_speed(i)= speed_after(i) - speed_before(i);
end


figure
scatter(speed_before,speed_after,40,'filled')
alpha 0.5
h=refline(1,0)
h.LineStyle='--';
h.Color=[105,105,105]/255;
xlabel({'Mean speed before'; 'sequence onset (cm/s)'});
ylabel({'Mean speed after'; 'sequence onset (cm/s)'});
set(gca,'fontsize',16);

min_delta_speed=min(delta_speed);
max_delta_speed=max(delta_speed);

[h,p,ci,stats]=ttest(delta_speed);
[hw,pw,statsw]=signrank(delta_speed');




aux_10=find(iwi_duration>10);
count=0;
for i=1:size(aux_10,2)
    count=count+1;
   speed_before_10(count)=mean(speed_iwi(aux_10(i),1:bins_number_epoch));
   speed_after_10(count)=mean(speed_iwi(aux_10(i),bins_number_epoch+1:end));
   delta_speed_10(count)= speed_after_10(count) - speed_before_10(count);
end


figure
scatter(speed_before_10,speed_after_10,40,'filled')
alpha 0.5
h=refline(1,0)
h.LineStyle='--';
h.Color=[105,105,105]/255;
xlabel({'Mean speed before'; 'wave onset (cm/s)'});
ylabel({'Mean speed after'; 'wave onset (cm/s)'});
set(gca,'fontsize',16);

min_delta_speed=min(delta_speed_10);
max_delta_speed=max(delta_speed_10);

[h,p,ci,stats]=ttest(delta_speed_10);
[h,p,stats]=signrank(delta_speed_10');







% aux_100=find(iwi_duration>100);
% count=0;
% for i=1:size(aux_100,2)
%     count=count+1;
%    speed_before_100(count)=mean(speed_iwi(aux_100(i),1:bins_number_epoch));
%    speed_after_100(count)=mean(speed_iwi(aux_100(i),bins_number_epoch+1:end));
%    delta_speed_100(count)= speed_after_100(count) - speed_before_100(count);
% end
% 
% figure
% scatter(speed_before_100,speed_after_100)
% refline(1,0)



%         iwi_duration(count_iwi)=(table_u(i,1)-table_u(i-1,2))/sf;
%         speed_iwi(count_iwi,:)=speed_d2(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed);
%         ac_iwi(count_iwi,:)=ac_d(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed);
%         track_session_iwi(count_iwi)=w;

%% Speed per IWI - 2 s

bins_number_epoch=floor(size(speed_iwi_2,2)/2);

for i=1:size(iwi_duration,2)
   speed_before(i)=mean(speed_iwi_2(i,1:bins_number_epoch));
   speed_after(i)=mean(speed_iwi_2(i,bins_number_epoch+1:end));
   delta_speed(i)= speed_after(i) - speed_before(i);
end


figure
scatter(speed_before,speed_after,40,'filled')
alpha 0.5
h=refline(1,0)
h.LineStyle='--';
h.Color=[105,105,105]/255;
xlabel({'Mean speed before'; 'sequence onset (cm/s)'});
ylabel({'Mean speed after'; 'sequence onset (cm/s)'});
set(gca,'fontsize',16);

min_delta_speed=min(delta_speed);
max_delta_speed=max(delta_speed);

[h,p,ci,stats]=ttest(delta_speed);




aux_10=find(iwi_duration>10);
count=0;
for i=1:size(aux_10,2)
    count=count+1;
   speed_before_10(count)=mean(speed_iwi(aux_10(i),1:bins_number_epoch));
   speed_after_10(count)=mean(speed_iwi(aux_10(i),bins_number_epoch+1:end));
   delta_speed_10(count)= speed_after_10(count) - speed_before_10(count);
end


figure
scatter(speed_before_10,speed_after_10,40,'filled')
alpha 0.5
h=refline(1,0)
h.LineStyle='--';
h.Color=[105,105,105]/255;
xlabel({'Mean speed before'; 'wave onset (cm/s)'});
ylabel({'Mean speed after'; 'wave onset (cm/s)'});
set(gca,'fontsize',16);

min_delta_speed=min(delta_speed_10);
max_delta_speed=max(delta_speed_10);

[h,p,ci,stats]=ttest(delta_speed_10);



% aux_100=find(iwi_duration>100);
% count=0;
% for i=1:size(aux_100,2)
%     count=count+1;
%    speed_before_100(count)=mean(speed_iwi(aux_100(i),1:bins_number_epoch));
%    speed_after_100(count)=mean(speed_iwi(aux_100(i),bins_number_epoch+1:end));
%    delta_speed_100(count)= speed_after_100(count) - speed_before_100(count);
% end
% 
% figure
% scatter(speed_before_100,speed_after_100)
% refline(1,0)



%         iwi_duration(count_iwi)=(table_u(i,1)-table_u(i-1,2))/sf;
%         speed_iwi(count_iwi,:)=speed_d2(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed);
%         ac_iwi(count_iwi,:)=ac_d(idx_iwi(j)-width_win_speed:idx_iwi(j)+width_win_speed);
%         track_session_iwi(count_iwi)=w;

%% Speed per wave and relationship between wave duration and speed


count_wave=0;
for i=1:10
    for j=1:size(wave_epoch_speed,1)
        if ~isempty(wave_epoch_speed{j,i})
        count_wave=count_wave+1;
        tracking_session(count_wave)=i;

        speed_mean_wave(count_wave)=mean(wave_epoch_speed{j,i});
        max_speed(count_wave)=max(wave_epoch_speed{j,i});
        speed_med_wave(count_wave)=median(wave_epoch_speed{j,i});
        speed_sem_wave(count_wave)=std(wave_epoch_speed{j,i})/sqrt(length(wave_epoch_speed{j,i}));
        speed_sd_wave(count_wave)=std(wave_epoch_speed{j,i});
        ac_mean_wave(count_wave)=mean(wave_epoch_ac{j,i});
        ac_median_wave(count_wave)=median(wave_epoch_ac{j,i});
        ac_sem_wave(count_wave)=std(wave_epoch_ac{j,i})/sqrt(length(wave_epoch_ac{j,i}));
        ac_sd_wave(count_wave)=std(wave_epoch_ac{j,i});

        end
    end
    
end

median(speed_mean_wave)
a=find(speed_mean_wave<prctile(speed_mean_wave,25));
b=find(speed_mean_wave>prctile(speed_mean_wave,75));

histogram_dur_low=histcounts(duration_waves_pooled(a),0:10:200);
histogram_dur_high=histcounts(duration_waves_pooled(b),0:10:200);

figure
plot(histogram_dur_low)
hold on
plot(histogram_dur_high)
legend('low','high')
ylabel('Counts');
xlabel('Wave duration (s)');
box off
set(gca,'fontsize',16);

figure
hold on
for i=1:10
    aux=find(tracking_session==i);
%     subplot(2,5,i)
hold on
%     scatter(duration_waves_pooled(aux),speed_sem_wave(aux));
     scatter(mean(duration_waves_pooled(aux)),mean(speed_med_wave(aux)),'filled');
end

figure
hold on
for i=1:10
    aux=find(tracking_session==i);
%     subplot(2,5,i)
hold on
%     scatter(duration_waves_pooled(aux),speed_sem_wave(aux));
     scatter(mean(duration_waves_pooled(aux)),mean(speed_mean_wave(aux)),'filled');
end
ylabel('Mean speed across waves per session (cm/s)');
xlabel('Mean wave duration per session (s)');
set(gca,'fontsize',16)
box off

% [r,p] = corr(duration_waves_pooled(aux)',speed_mean_wave(aux)','Type','Spearman');

%     for j=1:size(wave_epoch_speed,1)
%         if ~isempty(wave_epoch_speed{j,i})
%         count_wave=count_wave+1;
%         tracking_session(count_wave)=i;
% 
%         speed_mean_wave(count_wave)=mean(wave_epoch_speed{j,i});
%         speed_sem_wave(count_wave)=std(wave_epoch_speed{j,i})/sqrt(length(wave_epoch_speed{j,i}));
%         speed_sd_wave(count_wave)=std(wave_epoch_speed{j,i});
%         ac_mean_wave(count_wave)=mean(wave_epoch_ac{j,i});
%         ac_sem_wave(count_wave)=std(wave_epoch_ac{j,i})/sqrt(length(wave_epoch_ac{j,i}));
%         ac_sd_wave(count_wave)=std(wave_epoch_ac{j,i});
% 
%         end
%     end

% [r_1(i),p_1(i)] = corr(duration_waves_pooled(aux)',speed_mean_wave(aux)','Type','Spearman');
% [r_2(i),p_2(i)] = corr(duration_waves_pooled(aux)',speed_sd_wave(aux)','Type','Spearman');
% [r_3(i),p_3(i)] = corr(duration_waves_pooled(aux)',ac_mean_wave(aux)','Type','Spearman');
% [r_4(i),p_4(i)] = corr(duration_waves_pooled(aux)',ac_sd_wave(aux)','Type','Spearman');

    clear aux


figure
scatter(1:length(speed_mean_wave),speed_mean_wave)
ylabel('Mean speed per wave (cm/s)');
xlabel('Wave index');
set(gca,'fontsize',16)
box off


[r_total,p_total] = corr(duration_waves_pooled(aux)',speed_mean_wave(aux)','Type','Spearman');
[r_total,p_total] = corr(duration_waves_pooled',speed_med_wave','Type','Spearman');

figure
for i=1:10
    subplot(2,5,i)
    aux=find(tracking_session==i);
    scatter(duration_waves_pooled(aux),speed_mean_wave(aux),40,[0.8500, 0.3250, 0.0980],'filled');
    P = polyfit(duration_waves_pooled(aux),speed_mean_wave(aux),1);
    yfit = P(1)*duration_waves_pooled(aux)+P(2);
    hold on;
    plot(duration_waves_pooled(aux),yfit,'r-.');
    xlabel('Wave duration (s)');
    ylabel('Mean speed per wave (cm/s)');
end

figure
scatter(duration_waves_pooled,max_speed,80,[0.8500, 0.3250, 0.0980],'filled');
xlabel('Wave duration (s)');
ylabel('Max speed per wave (cm/s)');

figure
scatter(duration_waves_pooled,speed_med_wave,80,[0.8500, 0.3250, 0.0980],'filled');
xlabel('Wave duration (s)');
ylabel('Median speed per wave (cm/s)');

figure
scatter(duration_waves_pooled,speed_sem_wave,80,[0.8500, 0.3250, 0.0980],'filled');
xlabel('Wave duration (s)');
ylabel('SEM speed per wave (cm/s)');

figure
scatter(duration_waves_pooled,speed_sd_wave,80,[0.8500, 0.3250, 0.0980],'filled');
xlabel('Wave duration (s)');
ylabel('SD speed per wave (cm/s)');

figure
scatter(speed_mean_wave,speed_sd_wave,80,[0.8500, 0.3250, 0.0980],'filled');
xlabel('Mean speed per wave (cm/s)');
ylabel('SD speed per wave (cm/s)');

%% Probability of waves for different speed values

clear epochs
bins_speed=10;
figure
for w=1:10
    frames=find(session_pooled==w);
    wave_v=full_wave_vec(frames);
    phase_v=phase_pooled(frames);
    speed_v=full_speed(frames);
    ac_v=full_ac(frames);
   
    speed_edges=[2:(max(speed_v)-2)/bins_speed:max(speed_v)];
%     phase_pooled_wave=phase_v(find(wave_v>0));
%     full_speed_wave=speed_v(find(wave_v>0));
%     full_ac_wave=ac_v(find(wave_v>0));
    
    
   
    points=nan(bins_speed,length(phase_v));
    for i=1:bins_speed-1
        a=find(speed_v>speed_edges(i));
        b=find(speed_v<=speed_edges(i+1));
        points(i,1:length(intersect(a,b)))=intersect(a,b);        
        clear a b        
        if i==bins_speed-1
            points(bins_speed,1:length(find(speed_v>speed_edges(bins_speed))))=find(speed_v>speed_edges(bins_speed));
        end
        clear a b        
    end

    %speed
    speed_val=nan(length(phase_v),bins_speed);
    for i=1:bins_speed
        
        clear points_per_bin
        points_per_bin=points(i,~isnan(points(i,:)));
        speed_val(1:length(points_per_bin),i)=(wave_v(points_per_bin));
        
        fraction_waves_speed(w,i)=sum(wave_v(points_per_bin))/length(points_per_bin);
    end
    

    clear  frames wave_v phase_v speed_v ac_v phase_pooled_wave full_speed_wave full_ac_wave points speed_val full_ac_wave 
end
% fraction_waves_speed
anova1(fraction_waves_speed)

mean_fraction_waves_speed=mean(fraction_waves_speed);
sem_fraction_waves_speed=std(fraction_waves_speed)/(sqrt(10));

figure
errorbar(mean_fraction_waves_speed,sem_fraction_waves_speed)


%% Conditional probabilities session by session

%Conditional probabilities
% subset_norunning=find(speed_d2==0);
% subset_running=find(speed_d2>0);

matrix_dynamics_running=[p_w_r',p_nw_r'];
matrix_dynamics_no_running=[p_w_nr',p_nw_nr'];
matrix_sequence=[p_w_r',p_w_nr'];


p_w_r_m=mean(p_w_r);
p_w_nr_m=mean(p_w_nr);
p_nw_r_m=mean(p_nw_r);
p_nw_nr_m=mean(p_nw_nr);

p_r_w_m=mean(p_r_w);
p_nr_w_m=mean(p_nr_w);
p_r_nw_m=mean(p_r_nw);
p_nr_nw_m=mean(p_nr_nw);

p_w_r_sem=std(p_w_r)/sqrt(10);
p_w_nr_sem=std(p_w_nr)/sqrt(10);
p_nw_r_sem=std(p_nw_r)/sqrt(10);
p_nw_nr_sem=std(p_nw_nr)/sqrt(10);

p_r_w_sem=std(p_r_w)/sqrt(10);
p_nr_w_sem=std(p_nr_w)/sqrt(10);
p_r_nw_sem=std(p_r_nw)/sqrt(10);
p_nr_nw_sem=std(p_nr_nw)/sqrt(10);

mat_prob_full=[(p_r_w_m),(p_r_nw_m);(p_nr_w_m),(p_nr_nw_m)];
sem_mat_prob_full=[(p_r_w_sem),(p_r_nw_sem);(p_nr_w_sem),(p_nr_nw_sem)];

mat_prob2_full=[(p_w_r_m),(p_nw_r_m);(p_w_nr_m),(p_nw_nr_m)];
sem_mat_prob2_full=[(p_w_r_sem),(p_nw_r_sem);(p_w_nr_sem),(p_nw_nr_sem)];

[H_wr_nwr,P_wr_nwr]=ttest2(p_w_r,p_nw_r);
[H_wnr_nwnr,P_wnr_nwnr]=ttest2(p_w_nr,p_nw_nr);
[H_rw_nrw,P_rw_nrw]=ttest2(p_r_w,p_nr_w);
[H_rnw_nrnw,P_rnw_nrnw]=ttest2(p_r_nw,p_nr_nw);

figure
hold on
% col=jet(4);
col=summer(4);
b=bar([1 10],[p_r_w_m,p_r_nw_m],0.2,'FaceColor',col(4,:));
b=bar([3.5 12.5],[p_nr_w_m,p_nr_nw_m],0.2,'FaceColor',col(1,:));
errorbar([1,3.5,10,12.5],[p_r_w_m,p_nr_w_m,p_r_nw_m,p_nr_nw_m],[p_r_w_sem,p_nr_w_sem,p_r_nw_sem,p_nr_nw_sem],'k.');
axis([0.5 2.5 0 1])
xticks([2 11])
xticklabels({'P(Behavior|Wave)';'P(Behavior|No wave)'})
box off
set(gca,'fontsize',15)
axis([-2 15 0 1])
legend('Running','Still')
legend('boxoff')
ylabel('Probability')
    
figure
hold on
col=magma(4);
% col=summer(4);
b=bar([1 10],[p_w_r_m,p_w_nr_m],0.2,'FaceColor',col(4,:));
b=bar([3.5 12.5],[p_nw_r_m,p_nw_nr_m],0.2,'FaceColor',col(1,:));
errorbar([1,3.5,10,12.5],[p_w_r_m,p_nw_r_m,p_w_nr_m,p_nw_nr_m],[p_w_r_sem,p_nw_r_sem,p_w_nr_sem,p_nw_nr_sem],'k.');
axis([0.5 2.5 0 1])
xticks([2 11])
xticklabels({'P(Dynamics|Running)';'P(Dynamics|Still)'})
box off
set(gca,'fontsize',15)
axis([-2 15 0 1])
legend('Wave','No wave')
legend('boxoff')
ylabel('Probability')

figure
boxplot([matrix_dynamics_running,matrix_dynamics_no_running]);
ylabel('Probability');
set(gca,'fontsize',16,'YColor','k','XColor','k');
xticks([1.5,3.5])
xticklabels({'Dynamics|Running', 'Dynamics |Immobility'});
ylim([0 1])
yticks([])
box off

figure
boxplot([matrix_sequence]);
ylabel('Probability');
set(gca,'fontsize',16,'YColor','k','XColor','k');
xticks([1,2])
xticklabels({'Oscillation|Running', 'Oscillation|Immobility'});
ylim([0 1])
yticks([0 0.2 0.4 0.6 0.8 1])
box off

% Now I run the ttest on the differences
[H_wr_nwr_diff,P_wr_nwr_diff,stats_r_diff]=ttest((p_w_r-p_nw_r));
[H_wnr_nwnr_diff,P_wnr_nwnr_diff]=ttest((p_w_nr-p_nw_nr));
[H_rw_nrw_diff,P_rw_nrw_diff]=ttest((p_r_w-p_nr_w));
[H_rnw_nrnw_diff,P_rnw_nrnw_diff]=ttest((p_r_nw-p_nr_nw));

%probabilitidades cruzadas
[H_wr_nwr_diff,P_wr_nwr_diff,stats_r]=signrank((p_w_r-p_nw_r));
[H_wnr_nwnr_diff,P_wnr_nwnr_diff,stats_nr]=signrank((p_w_nr-p_nw_nr));
[H_rw_nrw_diff,P_rw_nrw_diff,stats_w]=signrank((p_r_w-p_nr_w));
[H_rnw_nrnw_diff,P_rnw_nrnw_diff,stats_nw]=signrank((p_r_nw-p_nr_nw));

[H_w,P_w,statsw]=signrank(p_w_r,p_w_nr);
[H_nw,P_nw,statsnw]=signrank(p_nw_r,p_nw_nr);



% %% Figure about Phase VS Speed
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Phase and speed,
% 
% clear epochs
% bins_phase=10;
% bins_speed=10;
% size_bin_phase=2*pi/bins_phase;
% phase_edges=-pi:2*pi/bins_phase:pi;
% speed_val2=[];
% speed_val2_norm=[];
% speed_val_sh_2=[];
% speed_val_sh_2_norm=[];
% cont=0;
% figure
% for w=1:10
%     disp(w);
%     frames=find(session_pooled==w);
%     wave_v=full_wave_vec(frames);
%     phase_v=phase_pooled(frames); %phase during session w
%     speed_v=full_speed(frames); %speed during session w
%     ac_v=full_ac(frames); %acceleration during during session w
%     
%     phase_pooled_wave=phase_v(find(wave_v>0));%phase with waves
%     full_speed_wave=speed_v(find(wave_v>0));%speed with waves
%     full_ac_wave=ac_v(find(wave_v>0));%ac with waves
%     
%     %MI
%     table(:,1)=phase_pooled_wave';
%     table(:,2)=full_speed_wave';    
%     edges_speed=0:max(full_speed_wave)/bins_speed:max(full_speed_wave);
%     MI_phase_speed(w)=compute_MI_SGC(table,phase_edges,edges_speed);    
%     table_sh(:,1)=table(:,1);
%     for sh=1:200
%          table_sh(:,2)=circshift(table(:,2),randperm(size(table,1),1));
%          MI_phase_speed_sh(w,sh)=compute_MI_SGC(table_sh,phase_edges,edges_speed);        
%     end    
%     bias(w)=mean(MI_phase_speed_sh(w,:));
%  
%     %time points into phase bins
%     points=nan(bins_phase,length(phase_pooled_wave));
%     for i=1:bins_phase-1
%         a=find(phase_pooled_wave>phase_edges(i));
%         b=find(phase_pooled_wave<=phase_edges(i+1));
%         points(i,1:length(intersect(a,b)))=intersect(a,b);        
%         clear a b        
%         if i==bins_phase-1
%             points(bins_phase,1:length(find(phase_pooled_wave>phase_edges(bins_phase))))=find(phase_pooled_wave>phase_edges(bins_phase));
%         end
%         clear a b        
%     end
% 
%     %Speed
%     speed_val=nan(length(phase_pooled_wave),bins_phase);
%     for i=1:bins_phase
%         
%         clear points_per_bin
%         points_per_bin=points(i,~isnan(points(i,:)));
%         speed_val(1:length(points_per_bin),i)=(full_speed_wave(points_per_bin));
%         
%         speed_stats(i,1,w)=nanmean(full_speed_wave(points_per_bin)./max(max(speed_val)));
%         speed_stats(i,2,w)=std(full_speed_wave(points_per_bin)./max(max(speed_val))/sqrt(length(points_per_bin)));
%     end
%     
%     %Shuffle speed
%     N_sh=500;
%     speed_val_sh=nan(length(phase_pooled_wave),bins_phase,N_sh);
%     for sh=1:N_sh
%         full_speed_wave_d_sh=circshift(full_speed_wave,randperm(size(table,1),1));
%         for i=1:bins_phase
%             
%             clear points_per_bin
%             points_per_bin=points(i,~isnan(points(i,:)));
%             speed_val_sh(1:length(points_per_bin),i,sh)=(full_speed_wave_d_sh(points_per_bin));
%             speed_stats_sh(i,1,sh+cont)=nanmean(full_speed_wave_d_sh(points_per_bin)./max(max(speed_val)));
%             speed_stats_sh(i,2,sh+cont)=std(full_speed_wave_d_sh(points_per_bin)./max(max(speed_val))/sqrt(length(points_per_bin)));
%         end
%     end
%     speed_stats_sh_w{1,w}=speed_stats_sh(:,:,cont+1:cont+N_sh);
% 
%     cont=cont+N_sh;
%     for i=1:bins_phase
%         for j=1:bins_phase            
%             [t_h_mat(i,j,w),p_h_mat(i,j,w)]=ttest2(speed_val(:,i),speed_val(:,j));            
%         end
%     end
%     
%     speed_val2=[speed_val2;speed_val];
%     speed_val2_norm=[speed_val2_norm;speed_val./max(max(speed_val))];
%   
% 
%     clear  frames wave_v phase_v speed_v ac_v phase_pooled_wave full_speed_wave full_ac_wave points speed_val full_ac_wave...
%         table table_sh phase_pooled_wave_d speed_pooled_wave_d
% end
% 
% 
% figure
% for w=1:10
%     
%     subplot(2,5,w);
%     bar(speed_stats(:,1,w));
%     hold on
%     errorbar(speed_stats(:,1,w),speed_stats(:,2,w),'*')  
% 
% end
% 
% % Sessions separately
% for w=1:10
%     for b=1:bins_phase
%         data_real=speed_stats(i,1,w);
%         data_shuffle=speed_stats_sh_w{1,w}(i,1,:);
%         [h_s(w,b),p_s(w,b)]=ttest2(data_real,data_shuffle);
%     end
%     
% end
% 
% for w=1:10    
%     for i=1:bins_phase
%         pooled_mean_s(w,i)=nanmean(speed_stats(i,1,w));
% %         pooled_sem_s(w,i)=nanstd(speed_stats(i,1,w))./sqrt(size(speed_stats,3));
%         
%         pooled_mean_sh_s(w,i)=nanmean(speed_stats_sh_w{1,w}(i,1,:));
%         pooled_sem_sh_s(w,i)=nanstd(speed_stats_sh_w{1,w}(i,1,:))./sqrt(size(speed_stats_sh_w{1,w}(i,1,:),3));
%     end
% end
% 
% figure 
% for i=1:10
%     subplot(2,5,i)
%     hold on
%     errorbar(pooled_mean_s(i,:),[],'k','linewidth',2.5)
%     errorbar(pooled_mean_sh_s(i,:),pooled_sem_sh_s(i,:),'r','linewidth',2.5)
%     axis([0.5,10.5,0,1])
%     ylabel('Normalized speed (cm/s)');
%     xlabel('Phase (rad)');
%     xticks([1,5.5,10]);
%     xticklabels({'-\pi','0','\pi'})
%     box off
%     set(gca,'fontsize',16);
% end
% 
% 
% %Pool sessions
% 
% for b=1:bins_phase
% data_real=speed_stats(b,1,:);
% data_shuffle=speed_stats_sh(b,1,:);
% [h(b),p(b),stast_speed(b)]=ranksum(data_real(:),data_shuffle(:));
% clear data_real data_shuffle
% end
% 
% for i=1:bins_phase
% pooled_mean(i)=nanmean(speed_stats(i,1,:));
% pooled_sem(i)=nanstd(speed_stats(i,1,:))./sqrt(size(speed_stats,3));
% 
% pooled_mean_sh(i)=nanmean(speed_stats_sh(i,1,:));
% pooled_sem_sh(i)=nanstd(speed_stats_sh(i,1,:))./sqrt(size(speed_stats_sh,3));
% end
% 
% 
% figure 
% hold on
% errorbar(pooled_mean,pooled_sem,'k','linewidth',2.5)
% errorbar(pooled_mean_sh,pooled_sem_sh,'r','linewidth',2.5)
% axis([0.5,10.5,0,1])
% ylabel('Normalized speed (cm/s)');
% xlabel('Phase (rad)');
% xticks([1,5.5,10]);
% xticklabels({'-\pi','0','\pi'})
% box off
% set(gca,'fontsize',16);
% axis([-inf inf 0 0.5]);
% 
% figure
% scatter(bias,MI_phase_speed,'o','filled')
% hline = refline([1,0]);
% hline.Color = 'r';
% ylabel('');
% xlabel('');
% 
% MI_corrected_phase_speed=MI_phase_speed-bias;
% 
% %% Phase and Acceleration
% 
% clear epochs table table_sh
% bins_phase=10;
% bins_ac=10;
% size_bin_phase=2*pi/bins_phase;
% phase_edges=-pi:2*pi/bins_phase:pi;
% ac_val2=[];
% ac_val2_norm=[];
% ac_val_sh_2=[];
% ac_val_sh_2_norm=[];
% cont=0;
% figure
% for w=1:10
%     disp(w);
%     frames=find(session_pooled==w);
%     wave_v=full_wave_vec(frames);
%     phase_v=phase_pooled(frames); %phase during session w
%     speed_v=full_speed(frames); %speed during session w
%     ac_v=full_ac(frames); %acceleration during during session w
%     
%     phase_pooled_wave=phase_v(find(wave_v>0));%phase with waves
%     full_speed_wave=speed_v(find(wave_v>0));%speed with waves
%     full_ac_wave=ac_v(find(wave_v>0));%ac with waves
%     
%     %MI   
%     table(:,1)=phase_pooled_wave';
%     table(:,2)=full_ac_wave';
%     
%     edges_ac=min(full_ac_wave):(max(full_ac_wave)-min(full_ac_wave))/bins_ac:max(full_ac_wave);
%     MI_phase_ac(w)=compute_MI_SGC(table,phase_edges,edges_ac);
%     
%     table_sh(:,1)=table(:,1);
%     for sh=1:200
%          table_sh(:,2)=circshift(table(:,2),randperm(size(table,1),1));
%          MI_phase_ac_sh(w,sh)=compute_MI_SGC(table_sh,phase_edges,edges_ac);        
%     end    
%     bias_ac(w)=mean(MI_phase_ac_sh(w,:));
%  
%     %time points into phase bins
%     points=nan(bins_phase,length(phase_pooled_wave));
%     for i=1:bins_phase-1
%         a=find(phase_pooled_wave>phase_edges(i));
%         b=find(phase_pooled_wave<=phase_edges(i+1));
%         points(i,1:length(intersect(a,b)))=intersect(a,b);        
%         clear a b        
%         if i==bins_phase-1
%             points(bins_phase,1:length(find(phase_pooled_wave>phase_edges(bins_phase))))=find(phase_pooled_wave>phase_edges(bins_phase));
%         end
%         clear a b        
%     end
% 
%     %Ac
%     ac_val=nan(length(phase_pooled_wave),bins_phase);
%     for i=1:bins_phase        
%         clear points_per_bin
%         points_per_bin=points(i,~isnan(points(i,:)));
%         ac_val(1:length(points_per_bin),i)=(full_ac_wave(points_per_bin));
%         
%         ac_stats(i,1,w)=nanmean(full_ac_wave(points_per_bin)./max(max(ac_val)));
%         ac_stats(i,2,w)=std(full_ac_wave(points_per_bin)./max(max(ac_val))/sqrt(length(points_per_bin)));
%     end
%     
%     %shuffle Ac
%     N_sh=500;
%     ac_val_sh=nan(length(phase_pooled_wave),bins_phase,N_sh);
%     for sh=1:N_sh
%         full_ac_wave_d_sh=circshift(full_ac_wave,randperm(size(table,1),1));
%         for i=1:bins_phase
%             
%             clear points_per_bin
%             points_per_bin=points(i,~isnan(points(i,:)));
%             ac_val_sh(1:length(points_per_bin),i,sh)=(full_ac_wave_d_sh(points_per_bin));
%             ac_stats_sh(i,1,sh+cont)=nanmean(full_ac_wave_d_sh(points_per_bin)./max(max(abs(ac_val))));
%             ac_stats_sh(i,2,sh+cont)=std(full_ac_wave_d_sh(points_per_bin)./max(max(abs(ac_val)))/sqrt(length(points_per_bin)));           
% 
%         end
% 
%     end
%     ac_stats_sh_w{1,w}=ac_stats_sh(:,:,cont+1:cont+N_sh);
%     cont=cont+N_sh;
%        
%     clear  frames wave_v phase_v speed_v ac_v phase_pooled_wave full_ac_wave full_ac_wave points speed_val full_ac_wave...
%         table table_sh phase_pooled_wave_d ac_pooled_wave_d
% end
% 
% 
% figure
% for w=1:10    
%     subplot(2,5,w);
%     bar(ac_stats(:,1,w));
%     hold on
%     errorbar(ac_stats(:,1,w),ac_stats(:,2,w),'*')  
% end
% 
% % Sessions separately
% for w=1:10
%     for b=1:bins_phase
%         data_real=ac_stats(i,1,w);
%         data_shuffle=ac_stats_sh_w{1,w}(i,1,:);
%         [h_s_ac(w,b),p_s_ac(w,b)]=ranksum(data_real(:),data_shuffle(:));
%     end
%     
% end
% 
% for w=1:10    
%     for i=1:bins_phase
%         pooled_mean_s_ac(w,i)=nanmean(ac_stats(i,1,w));
% %         pooled_sem_s(w,i)=nanstd(speed_stats(i,1,w))./sqrt(size(speed_stats,3));        
%         pooled_mean_sh_s_ac(w,i)=nanmean(ac_stats_sh_w{1,w}(i,1,:));
%         pooled_sem_sh_s_ac(w,i)=nanstd(ac_stats_sh_w{1,w}(i,1,:))./sqrt(size(ac_stats_sh_w{1,w}(i,1,:),3));
%     end
% end
% 
% figure 
% for i=1:10
%     subplot(2,5,i)
%     hold on
%     errorbar(pooled_mean_s_ac(i,:),[],'k','linewidth',2.5)
%     errorbar(pooled_mean_sh_s_ac(i,:),pooled_sem_sh_s_ac(i,:),'r','linewidth',2.5)
% %     axis([-1,1,0,1])
%     ylabel('Normalized Ac (cm2/s)');
%     xlabel('Phase (rad)');
%     xticks([1,5.5,10]);
%     xticklabels({'-\pi','0','\pi'})
%     box off
%     set(gca,'fontsize',16);
% end
% 
% 
% %Pool sessions % Here I am pooling shuffling. Perhaps I should pool means
% %of shufflings
% 
% for b=1:bins_phase
% data_real=ac_stats(b,1,:);
% data_shuffle=ac_stats_sh(b,1,:);
% [h_ac(b),p_ac(b),stats_ac(b)]=ranksum(data_real(:),data_shuffle(:));
% end
% 
% for i=1:bins_phase
% pooled_mean_ac(i)=nanmean(ac_stats(i,1,:));
% pooled_sem_ac(i)=nanstd(ac_stats(i,1,:))./sqrt(size(ac_stats,3));
% 
% pooled_mean_sh_ac(i)=nanmean(ac_stats_sh(i,1,:));
% pooled_sem_sh_ac(i)=nanstd(ac_stats_sh(i,1,:))./sqrt(size(ac_stats_sh,3));
% end
% 
% 
% figure 
% hold on
% errorbar(pooled_mean_ac,pooled_sem_ac,'k','linewidth',2.5)
% errorbar(pooled_mean_sh_ac,pooled_sem_sh_ac,'r','linewidth',2.5)
% axis([-inf inf -0.03 0.03])
% ylabel('Normalized acceleration (cm/s)');
% xlabel('Phase (rad)');
% xticks([1,5.5,10]);
% xticklabels({'-\pi','0','\pi'})
% box off
% set(gca,'fontsize',16);
% 
% figure
% scatter(bias_ac,MI_phase_ac,'o','filled')
% hline = refline([1,0]);
% hline.Color = 'r';
% ylabel('');
% xlabel('');
% 
% MI_corrected_phase_ac=MI_phase_ac-bias_ac;
% 
% 
% %% Phase and position
% 
% clear epochs table table_sh
% bins_phase=10;
% bins_position=10;
% size_bin_phase=2*pi/bins_phase;
% phase_edges=-pi:2*pi/bins_phase:pi;
% 
% cont=0;
% figure
% for w=1:10
%     disp(w);
%     frames=find(session_pooled==w);
%     wave_v=full_wave_vec(frames);
%     phase_v=phase_pooled(frames); %phase during session w
%     pos_v=full_position(frames); %pos during session w
%     
%     phase_pooled_wave=phase_v(find(wave_v>0));%phase with waves
%     full_pos_wave=pos_v(find(wave_v>0));%speed with waves
%     
%     %MI   
%     table(:,1)=phase_pooled_wave';
%     table(:,2)=full_pos_wave';
%     
%     edges_pos=0:54/bins_position:54;
%     MI_phase_pos(w)=compute_MI_SGC(table,phase_edges,edges_pos);
%     
%     table_sh(:,1)=table(:,1);
%     for sh=1:200
%          table_sh(:,2)=circshift(table(:,2),randperm(size(table,1),1));
%          MI_phase_pos_sh(w,sh)=compute_MI_SGC(table_sh,phase_edges,edges_pos);        
%     end    
%     bias_pos(w)=mean(MI_phase_pos_sh(w,:));
%  
%     %time points into phase bins
%     %points has, in each colum, the frames where the phase was within
%     %certain values
%     points=nan(bins_phase,length(phase_pooled_wave));
%     for i=1:bins_phase-1
%         a=find(phase_pooled_wave>phase_edges(i));
%         b=find(phase_pooled_wave<=phase_edges(i+1));
%         points(i,1:length(intersect(a,b)))=intersect(a,b);        
%         clear a b        
%         if i==bins_phase-1
%             points(bins_phase,1:length(find(phase_pooled_wave>phase_edges(bins_phase))))=find(phase_pooled_wave>phase_edges(bins_phase));
%         end
%         clear a b        
%     end
% 
%     %Pos
%     pos_val=nan(length(phase_pooled_wave),bins_phase);
%     for i=1:bins_phase        
%         clear points_per_bin
%         points_per_bin=points(i,~isnan(points(i,:)));
%         pos_val(1:length(points_per_bin),i)=(full_pos_wave(points_per_bin));
%         
%         pos_stats(i,1,w)=nanmean(full_pos_wave(points_per_bin)./max(max(pos_val)));
%         pos_stats(i,2,w)=std(full_pos_wave(points_per_bin)./max(max(pos_val))/sqrt(length(points_per_bin)));
%     end
%     
%     %shuffle Pos
%     N_sh=500;
%     pos_val_sh=nan(length(phase_pooled_wave),bins_phase,N_sh);
%     for sh=1:N_sh
%         full_pos_wave_d_sh=circshift(full_pos_wave,randperm(size(table,1),1));
%         for i=1:bins_phase
%             
%             clear points_per_bin
%             points_per_bin=points(i,~isnan(points(i,:)));
%             pos_val_sh(1:length(points_per_bin),i,sh)=(full_pos_wave_d_sh(points_per_bin));
%             pos_stats_sh(i,1,sh+cont)=nanmean(full_pos_wave_d_sh(points_per_bin)./max(max(abs(pos_val))));
%             pos_stats_sh(i,2,sh+cont)=std(full_pos_wave_d_sh(points_per_bin)./max(max(abs(pos_val)))/sqrt(length(points_per_bin)));           
% 
%         end
% 
%     end
%     pos_stats_sh_w{1,w}=pos_stats_sh(:,:,cont+1:cont+N_sh);
%     cont=cont+N_sh;
%        
%     clear  frames wave_v phase_v speed_v ac_v phase_pooled_wave full_ac_wave full_ac_wave points speed_val full_ac_wave...
%         table table_sh phase_pooled_wave_d ac_pooled_wave_d
% end
% 
% 
% figure
% for w=1:10    
%     subplot(2,5,w);
%     bar(pos_stats(:,1,w));
%     hold on
%     errorbar(pos_stats(:,1,w),pos_stats(:,2,w),'*')  
% end
% 
% % Sessions separately
% for w=1:10
%     for b=1:bins_phase
%         data_real=pos_stats(i,1,w);
%         data_shuffle=pos_stats_sh_w{1,w}(i,1,:);
%         [h_s_ac(w,b),p_s_ac(w,b)]=ranksum(data_real(:),data_shuffle(:));
%     end
%     
% end
% 
% for w=1:10    
%     for i=1:bins_phase
%         pooled_mean_s_pos(w,i)=nanmean(pos_stats(i,1,w));
% %         pooled_sem_s(w,i)=nanstd(speed_stats(i,1,w))./sqrt(size(speed_stats,3));        
%         pooled_mean_sh_s_pos(w,i)=nanmean(pos_stats_sh_w{1,w}(i,1,:));
%         pooled_sem_sh_s_pos(w,i)=nanstd(pos_stats_sh_w{1,w}(i,1,:))./sqrt(size(pos_stats_sh_w{1,w}(i,1,:),3));
%     end
% end
% 
% figure 
% for i=1:10
%     subplot(2,5,i)
%     hold on
%     errorbar(pooled_mean_s_pos(i,:),[],'k','linewidth',2.5)
%     errorbar(pooled_mean_sh_s_pos(i,:),pooled_sem_sh_s_pos(i,:),'r','linewidth',2.5)
% %     axis([-1,1,0,1])
%     ylabel('Normalized Pos (cm)');
%     xlabel('Phase (rad)');
%     xticks([1,5.5,10]);
%     xticklabels({'-\pi','0','\pi'})
%     box off
%     set(gca,'fontsize',16);
% end
% 
% 
% %Pool sessions
% 
% for b=1:bins_phase
% data_real=pos_stats(b,1,:);
% data_shuffle=pos_stats_sh(b,1,:);
% [h_pos(b),p_pos(b),stats_pos(b)]=ranksum(data_real(:),data_shuffle(:));
% end
% 
% for i=1:bins_phase
% pooled_mean_pos(i)=nanmean(pos_stats(i,1,:));
% pooled_sem_pos(i)=nanstd(pos_stats(i,1,:))./sqrt(size(pos_stats,3));
% 
% pooled_mean_sh_pos(i)=nanmean(pos_stats_sh(i,1,:));
% pooled_sem_sh_pos(i)=nanstd(pos_stats_sh(i,1,:))./sqrt(size(pos_stats_sh,3));
% end
% 
% 
% figure 
% hold on
% errorbar(pooled_mean_pos,pooled_sem_pos,'k','linewidth',2.5)
% errorbar(pooled_mean_sh_pos,pooled_sem_sh_pos,'r','linewidth',2.5)
% axis([-inf inf 0.45 0.6])
% ylabel('Normalized position (cm)');
% xlabel('Phase (rad)');
% xticks([1,5.5,10]);
% xticklabels({'-\pi','0','\pi'})
% box off
% set(gca,'fontsize',16);
% 
% figure
% scatter(bias_pos,MI_phase_pos,'o','filled')
% hline = refline([1,0]);
% hline.Color = 'r';
% ylabel('');
% xlabel('');
% 
% MI_corrected_phase_pos=MI_phase_pos-bias_pos;
% 
% %% MI
% 
% figure 
% hold on
% plot(MI_corrected_phase_pos,'-*');
% plot(MI_corrected_phase_speed,'-*');
% plot(MI_corrected_phase_ac,'-*');
% legend('Pos','Speed','Ac')
% ylabel('MI');
% xlabel('Session');
% 
% 
% figure
% hold on
% histogram(MI_corrected_phase_pos,[0:0.01:0.4]);
% histogram(MI_corrected_phase_speed,[0:0.01:0.4]);
% histogram(MI_corrected_phase_ac,[0:0.01:0.4]);
% 
% figure
% % plot([1,2,3];[MI_corrected_phase_pos';MI_corrected_phase_speed';MI_corrected_phase_ac'])