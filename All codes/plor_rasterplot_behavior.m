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
                    big_table(count,12)=size(spk.spikes_d_s,1); %N
                    
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


%% Data and running epochs for MEC

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


speed_smoothing_kernel=3;
threshold_still=2; %if speed<threshold_still then speed=0
R=8.5145; %Radius of the wheel in cm
distance_ran_MEC=zeros(1,length(mec_sessions));

for w=1:length(mec_sessions)
    row_w=mec_sessions(w);

    if big_table(row_w,1)>5
        disp(w)
        
        count=count+1;
        mouse=['L',num2str(big_table(row_w,1)),'M',num2str(big_table(row_w,2))];
        day=big_table(row_w,3);
        s=big_table(row_w,4);
        munit=big_table(row_w,5);
        
        clus=10;
        disc_phase=10;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
        
        load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
        file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (isfield(dates,'folder_name')==1)
            day_index=find (dates.actual_day==day);
            file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day_index),'_MUnit',num2str(munit),'.mat']];
        else
            file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        end
        file_name_osf=[dpath ['spikes_30ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if(exist(file_name_spk)==2)
            
            %     load(file_name_dff,'-mat'); %DFF
            load(file_name_spk,'-mat'); %Spike times
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            load(file_name_osf,'-mat');
            spikes_d_osf=full(spikes_d_s);
            
            [~,T]=size(spikes_d);
            [N,T_osf]=size(spikes_d_osf);
            
            [~,sorting,~]=get_sorting(spikes_d);
            
            if N>150
                
                %disp(day)
                if length(num2str(dates.days(day)))==7
                    st=num2str(dates.days(day));
                    st_new=['0',st];
                    file_name_beh=[dbeh_path,mouse,'\Flavio_2P_',mouse,'_',st_new,'_MUnit_',num2str(munit),'_TRACKING.csv'];
                else
                    file_name_beh=[dbeh_path,mouse,'\Flavio_2P_',mouse,'_',num2str(dates.days(day)),'_MUnit_',num2str(munit),'_TRACKING.csv'];
                end
                
                alfa=table2array(readtable(file_name_beh)); %Table with info about position and time
                timestam=alfa(:,1);
                
                if big_table(row_w,1)>8
                    %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
                    %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
                    tot_position=alfa(:,3)/10; %in cm
                    lap_position=alfa(:,4)/10; %in cm
                    lap_index=alfa(:,5);
                    motor=alfa(:,6);
                    
                    %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
                    R=max(lap_position)/(2*pi); %in mm
                    speed=diff(tot_position)./diff(timestam); %cm/s % No smoothing
                    %                 ac=diff(speed)./diff(timestam(1:end-1));
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
                    
                else
                    tot_position=alfa(:,2)/10; %in cm
                    
                    %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
                    
                    speed=diff(tot_position)./diff(timestam); %cm/s % Smoothing kernel of 2.3 seconds
                    
                    %%%%%%%%%%%%%%%%%%%%%% Imaging time stamps and camera time steps
                    
                    times_s=0.0323:0.0323:(T_osf*0.0323);
                    times_s_d=downsample(times_s,4);
                    times_tracking=times_s_d(1:T); %Imaging time points
                    
                    %%%%%%%%%%%%%%%%%%%%%% Interpolated and downsampled quantities
                    tot_position_d=interp1(timestam(1:end),tot_position,times_tracking(1:end)); %Interpolated to imagining time points
                    lap_index_d = ceil(abs(tot_position_d./(2*pi*R)));
                    lap_index_d(find(lap_index_d==0))=1;
                    lap_position_d=tot_position_d-(2*pi*R.*(lap_index_d-1));
                    %         angle_d=lap_position_d./R;
                    
                    speed_d=interp1(timestam(1:end-1),speed,times_tracking(1:end-1)); %Interpolated to imagining time points
                    speed_d=smooth(speed_d,speed_smoothing_kernel);
                    speed_d2=speed_d;
                    speed_d2(speed_d2<threshold_still)=0;
                    speed_d2(end+1)=speed_d2(end);
                    ac_d=diff(speed_d2)./diff(times_tracking(1:end));
                    ac_d(end+1)=ac_d(end);
                end
                
                %                 motion=ones(1,T);
                %                 motion(speed_d2==0)=0;
                
                %                 %Rasterplot
                %                 figure;
                %                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                %                 hold on
                %                 for i=1:N
                %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting(i),:),4,'k','filled')
                %                     alpha 0.3
                %                 end
                %                 box on
                %                 h=area((1:size(spikes_d,2))./8,N.*motion,'LineStyle','-');
                %                 h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
                %                 h(1).FaceColor = [0.2 0.8 0.8];
                %                 h(1).FaceAlpha = 0.4;
                %                 %title([mouse,' Day',num2str(day)]);
                %                 axis([-inf inf 1 inf]);
                %                 yticks([100 400])
                %                 xlabel('Time (s)');
                %                 ylabel('Neurons #');
                %                 set(gca,'fontsize',20);
                %                 box off
                %
%                 figure;
%                 subplot(2,1,1)
%                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.5]);
%                 hold on
%                 for i=1:N
%                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting(i),:),4,'k','filled')
%                     %                     alpha 0.15
%                     alpha 0.2
%                 end
%                 box on
%                 %                 h=area((1:size(spikes_d,2))./8,N.*motion,'LineStyle','-');
%                 %                 h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
%                 %                 h(1).FaceColor = [0.2 0.8 0.8];
%                 %                 h(1).FaceAlpha = 0.4;
%                 %title([mouse,' Day',num2str(day)]);
%                 axis([-inf inf 1 inf]);
%                 xticks([]);
%                 %                 xlabel('Time (s)');
%                 ylabel('Neurons #');
%                 yticks([100 400])
%                 set(gca,'fontsize',16);
%                 title([mouse,'    -    Day=',num2str(day),'    -    Session=',num2str(s),'    -    MUnit=',num2str(munit),'    -    Distance=',num2str(floor(tot_position(end))),'cm' ]);
%                 %                 subplot(3,1,2)
%                 %                 plot((1:size(lap_position_d,1)-1)./8,speed_d,'k','linewidth',2)
%                 %                 ylabel({'Speed';'(cm/s)'});
%                 %                 set(gca,'fontsize',30)
%                 % %                 axis([620 780 -2 10])
%                 %                 yticks([0 10]);
%                 %                 xticks([]);
%                 %                 box off
%                 subplot(2,1,2)
%                 plot((1:size(lap_position_d,1))./8,lap_position_d,'k','linewidth',2)
%                 ylabel({'Position';'on wheel (cm)'});
%                 axis([-inf inf 1 inf]);
%                 %                 xticks([620 700 780])
%                 yticks([0 50]);
%                 box off
%                 xlabel('Time (s)');
%                 set(gca,'fontsize',16);
%                 
%                 saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July 2021\Raster plots and behavior\',mouse,'_Day',num2str(day),'_Session',num2str(s),'.svg'));
%                 saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July 2021\Raster plots and behavior\',mouse,'_Day',num2str(day),'_Session',num2str(s),'.fig'));
%                 
%                 
                distance_ran_MEC(w)= floor(tot_position(end));
                mean_speed(w)=mean(speed_d2);
                std_speed(w)=std(speed_d2);
                mean_ac(w)=mean(ac_d);
                std_ac(w)=std(ac_d);
                ws_binary(w)=big_table(row_w,6);
                ws_entropy(w)=big_table(row_w,11);
                
                close all
                
                clear ac_d alfa cells_d lap_index_d lap_position_d N sortind_descend sorting sorting_0 sorting_ascend speed speed_s speed_d2 spikes_d ...
                    spikes_d_osf spikes_sh time_s simt_s_d times_tracking timestam tot_position tot_position_d times_s times_s_d spikes_d_s speed_d ...
                    motion
            end
        end
    end
end


%% Data and running epochs for PaS

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
for w=5:length(repeated_days)
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

PaS_sessions_waves=length(find(big_table(PaS_sessions,6)>0));
PaS_sessions_nowaves=length(find(big_table(PaS_sessions,6)==0));


speed_smoothing_kernel=5;
threshold_still=2; %if speed<threshold_still then speed=0
R=8.5145; %Radius of the wheel in cm
distance_ran_PaS=zeros(1,length(PaS_sessions));

for w=1:length(PaS_sessions)
    row_w=PaS_sessions(w);

    if big_table(row_w,1)>5
        disp(w)
        
        count=count+1;
        mouse=['L',num2str(big_table(row_w,1)),'M',num2str(big_table(row_w,2))];
        day=big_table(row_w,3);
        s=big_table(row_w,4);
        munit=big_table(row_w,5);
        
        clus=10;
        disc_phase=10;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
        
        load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
        file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (isfield(dates,'folder_name')==1)
            day_index=find (dates.actual_day==day);
            file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day_index),'_MUnit',num2str(munit),'.mat']];
        else
            file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        end
        file_name_osf=[dpath ['spikes_30ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if(exist(file_name_spk)==2)
            
            %     load(file_name_dff,'-mat'); %DFF
            load(file_name_spk,'-mat'); %Spike times
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            load(file_name_osf,'-mat');
            spikes_d_osf=full(spikes_d_s);
            
            [~,T]=size(spikes_d);
            [N,T_osf]=size(spikes_d_osf);
            
            [~,sorting,~]=get_sorting(spikes_d);
            
            if N>150
%                 day_i=find (dates.actual_day==day);
                %disp(day)
                if length(num2str(dates.days(day)))==7
                    st=num2str(dates.days(day));
                    st_new=['0',st];
                    file_name_beh=[dbeh_path,mouse,'\Flavio_2P_',mouse,'_',st_new,'_MUnit_',num2str(munit),'_TRACKING.csv'];
                else
                    file_name_beh=[dbeh_path,mouse,'\Flavio_2P_',mouse,'_',num2str(dates.days(day)),'_MUnit_',num2str(munit),'_TRACKING.csv'];
                end


                alfa=table2array(readtable(file_name_beh)); %Table with info about position and time
                timestam=alfa(:,1);
                
                if big_table(row_w,1)>8
                    %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
                    %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
                    tot_position=alfa(:,3)/10; %in cm
                    lap_position=alfa(:,4)/10; %in cm
                    lap_index=alfa(:,5);
                    motor=alfa(:,6);
                    
                    %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
                    R=max(lap_position)/(2*pi); %in mm
                    speed=diff(tot_position)./diff(timestam); %cm/s % No smoothing
                    %                 ac=diff(speed)./diff(timestam(1:end-1));
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
                    
                else
                    tot_position=alfa(:,2)/10; %in cm
                    
                    %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
                    
                    speed=diff(tot_position)./diff(timestam); %cm/s % Smoothing kernel of 2.3 seconds
                    
                    %%%%%%%%%%%%%%%%%%%%%% Imaging time stamps and camera time steps
                    
                    times_s=0.0323:0.0323:(T_osf*0.0323);
                    times_s_d=downsample(times_s,4);
                    times_tracking=times_s_d(1:T); %Imaging time points
                    
                    %%%%%%%%%%%%%%%%%%%%%% Interpolated and downsampled quantities
                    tot_position_d=interp1(timestam(1:end),tot_position,times_tracking(1:end)); %Interpolated to imagining time points
                    lap_index_d = ceil(abs(tot_position_d./(2*pi*R)));
                    lap_index_d(find(lap_index_d==0))=1;
                    lap_position_d=tot_position_d-(2*pi*R.*(lap_index_d-1));
                    %         angle_d=lap_position_d./R;
                    
                    speed_d=interp1(timestam(1:end-1),speed,times_tracking(1:end-1)); %Interpolated to imagining time points
                    speed_d=smooth(speed_d,speed_smoothing_kernel);
                    speed_d2=speed_d;
                    speed_d2(speed_d2<threshold_still)=0;
                    speed_d2(end+1)=speed_d2(end);
                    ac_d=diff(speed_d2)./diff(times_tracking(1:end));
                    ac_d(end+1)=ac_d(end);
                end
                
                %                 motion=ones(1,T);
                %                 motion(speed_d2==0)=0;
                
                %                 %Rasterplot
                %                 figure;
                %                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                %                 hold on
                %                 for i=1:N
                %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting(i),:),4,'k','filled')
                %                     alpha 0.3
                %                 end
                %                 box on
                %                 h=area((1:size(spikes_d,2))./8,N.*motion,'LineStyle','-');
                %                 h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
                %                 h(1).FaceColor = [0.2 0.8 0.8];
                %                 h(1).FaceAlpha = 0.4;
                %                 %title([mouse,' Day',num2str(day)]);
                %                 axis([-inf inf 1 inf]);
                %                 yticks([100 400])
                %                 xlabel('Time (s)');
                %                 ylabel('Neurons #');
                %                 set(gca,'fontsize',20);
                %                 box off
                %
% %                 figure;
% %                 subplot(2,1,1)
% %                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.5]);
% %                 hold on
% %                 for i=1:N
% %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting(i),:),4,'k','filled')
% %                     %                     alpha 0.15
% %                     alpha 0.2
% %                 end
% %                 box on
% %                 %                 h=area((1:size(spikes_d,2))./8,N.*motion,'LineStyle','-');
% %                 %                 h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
% %                 %                 h(1).FaceColor = [0.2 0.8 0.8];
% %                 %                 h(1).FaceAlpha = 0.4;
% %                 %title([mouse,' Day',num2str(day)]);
% %                 axis([-inf inf 1 inf]);
% %                 xticks([]);
% %                 %                 xlabel('Time (s)');
% %                 ylabel('Neurons #');
% %                 yticks([100 400])
% %                 set(gca,'fontsize',16);
% %                 title([mouse,'    -    Day=',num2str(day),'    -    Session=',num2str(s),'    -    MUnit=',num2str(munit),'    -    Distance=',num2str(floor(tot_position(end))),'cm' ]);
% %                 %                 subplot(3,1,2)
% %                 %                 plot((1:size(lap_position_d,1)-1)./8,speed_d,'k','linewidth',2)
% %                 %                 ylabel({'Speed';'(cm/s)'});
% %                 %                 set(gca,'fontsize',30)
% %                 % %                 axis([620 780 -2 10])
% %                 %                 yticks([0 10]);
% %                 %                 xticks([]);
% %                 %                 box off
% %                 subplot(2,1,2)
% %                 plot((1:size(lap_position_d,1))./8,lap_position_d,'k','linewidth',2)
% %                 ylabel({'Position';'on wheel (cm)'});
% %                 axis([-inf inf 1 inf]);
% %                 %                 xticks([620 700 780])
% %                 yticks([0 50]);
% %                 box off
% %                 xlabel('Time (s)');
% %                 set(gca,'fontsize',16);
% %                 
% %                 saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July 2021\Raster plots and behavior\',mouse,'_Day',num2str(day),'_Session',num2str(s),'.svg'));
% %                 saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July 2021\Raster plots and behavior\',mouse,'_Day',num2str(day),'_Session',num2str(s),'.fig'));
% %                 
% % %                 
                distance_ran_PaS(w)= floor(tot_position(end));
                mean_speed_PaS(w)=mean(speed_d2);
                std_speed_PaS(w)=std(speed_d2);
                mean_ac_PaS(w)=mean(ac_d);
                std_ac_PaS(w)=std(ac_d);
                ws_binary_PaS(w)=big_table(row_w,6);
                ws_entropy_PaS(w)=big_table(row_w,11);
                
                close all
                
                clear ac_d alfa cells_d lap_index_d lap_position_d N sortind_descend sorting sorting_0 sorting_ascend speed speed_s speed_d2 spikes_d ...
                    spikes_d_osf spikes_sh time_s simt_s_d times_tracking timestam tot_position tot_position_d times_s times_s_d spikes_d_s speed_d ...
                    motion
            end
        end
    end
end


%% Data and running epochs for V1

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


speed_smoothing_kernel=5;
threshold_still=2; %if speed<threshold_still then speed=0
R=8.5145; %Radius of the wheel in cm
distance_ran_V1=zeros(1,length(V1_sessions));

for w=1:length(V1_sessions)
    row_w=V1_sessions(w);

    if big_table(row_w,1)>0
        disp(w)
        
        count=count+1;
        mouse_i=big_table(row_w,2);
        if mouse_i==1
            mouse='92227';
        elseif mouse_i==2
            mouse='92229';
        else
            mouse='60961';
        end
        
        day=big_table(row_w,3);
        s=big_table(row_w,4);
        munit=big_table(row_w,5);
        
        clus=10;
        disc_phase=10;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
        
        load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
        file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (isfield(dates,'folder_name')==1)
            day_index=find (dates.actual_day==day);
            file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day_index),'_MUnit',num2str(munit),'.mat']];
        else
            file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day_index),'_MUnit',num2str(munit),'.mat']];
        end
        file_name_osf=[dpath ['spikes_30ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day_index),'_MUnit',num2str(munit),'.mat']];
        
        if(exist(file_name_spk)==2)
            
            %     load(file_name_dff,'-mat'); %DFF
            load(file_name_spk,'-mat'); %Spike times
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            load(file_name_osf,'-mat');
            spikes_d_osf=full(spikes_d_s);
            
            [~,T]=size(spikes_d);
            [N,T_osf]=size(spikes_d_osf);
            
            [~,sorting,~]=get_sorting(spikes_d);
            
            if N>150
                
                %disp(day)
                if length(num2str(dates.days(day_index)))==7
                    st=num2str(dates.days(day_index));
                    st_new=['0',st];
                    file_name_beh=[dbeh_path,mouse,'\Flavio_2P_V1_Spring 2020_',mouse,'_',st_new,'_MUnit_',num2str(munit),'_TRACKING.csv'];
                else
                    file_name_beh=[dbeh_path,mouse,'\Flavio_2P_V1_Spring 2020_',mouse,'_',num2str(dates.days(day_index)),'_MUnit_',num2str(munit),'_TRACKING.csv'];
                end
                
                alfa=table2array(readtable(file_name_beh)); %Table with info about position and time
                timestam=alfa(:,1);
                
                if big_table(row_w,1)<5
                    %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
                    %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
                    tot_position=alfa(:,3)/10; %in cm
                    lap_position=alfa(:,4)/10; %in cm
                    lap_index=alfa(:,5);
                    motor=alfa(:,6);
                    
                    %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
                    R=max(lap_position)/(2*pi); %in mm
                    speed=diff(tot_position)./diff(timestam); %cm/s % No smoothing
                    %                 ac=diff(speed)./diff(timestam(1:end-1));
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
                    
                else
                    tot_position=alfa(:,2)/10; %in cm
                    
                    %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
                    
                    speed=diff(tot_position)./diff(timestam); %cm/s % Smoothing kernel of 2.3 seconds
                    
                    %%%%%%%%%%%%%%%%%%%%%% Imaging time stamps and camera time steps
                    
                    times_s=0.0323:0.0323:(T_osf*0.0323);
                    times_s_d=downsample(times_s,4);
                    times_tracking=times_s_d(1:T); %Imaging time points
                    
                    %%%%%%%%%%%%%%%%%%%%%% Interpolated and downsampled quantities
                    tot_position_d=interp1(timestam(1:end),tot_position,times_tracking(1:end)); %Interpolated to imagining time points
                    lap_index_d = ceil(abs(tot_position_d./(2*pi*R)));
                    lap_index_d(find(lap_index_d==0))=1;
                    lap_position_d=tot_position_d-(2*pi*R.*(lap_index_d-1));
                    %         angle_d=lap_position_d./R;
                    
                    speed_d=interp1(timestam(1:end-1),speed,times_tracking(1:end-1)); %Interpolated to imagining time points
                    speed_d=smooth(speed_d,speed_smoothing_kernel);
                    speed_d2=speed_d;
                    speed_d2(speed_d2<threshold_still)=0;
                    speed_d2(end+1)=speed_d2(end);
                    ac_d=diff(speed_d2)./diff(times_tracking(1:end));
                    ac_d(end+1)=ac_d(end);
                end
                
                %                 motion=ones(1,T);
                %                 motion(speed_d2==0)=0;
                
                %                 %Rasterplot
                %                 figure;
                %                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                %                 hold on
                %                 for i=1:N
                %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting(i),:),4,'k','filled')
                %                     alpha 0.3
                %                 end
                %                 box on
                %                 h=area((1:size(spikes_d,2))./8,N.*motion,'LineStyle','-');
                %                 h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
                %                 h(1).FaceColor = [0.2 0.8 0.8];
                %                 h(1).FaceAlpha = 0.4;
                %                 %title([mouse,' Day',num2str(day)]);
                %                 axis([-inf inf 1 inf]);
                %                 yticks([100 400])
                %                 xlabel('Time (s)');
                %                 ylabel('Neurons #');
                %                 set(gca,'fontsize',20);
                %                 box off
                %
% % %                 figure;
% % %                 subplot(2,1,1)
% % %                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.5]);
% % %                 hold on
% % %                 for i=1:N
% % %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting(i),:),4,'k','filled')
% % %                     %                     alpha 0.15
% % %                     alpha 0.2
% % %                 end
% % %                 box on
% % %                 %                 h=area((1:size(spikes_d,2))./8,N.*motion,'LineStyle','-');
% % %                 %                 h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
% % %                 %                 h(1).FaceColor = [0.2 0.8 0.8];
% % %                 %                 h(1).FaceAlpha = 0.4;
% % %                 %title([mouse,' Day',num2str(day)]);
% % %                 axis([-inf inf 1 inf]);
% % %                 xticks([]);
% % %                 %                 xlabel('Time (s)');
% % %                 ylabel('Neurons #');
% % %                 yticks([100 200])
% % %                 set(gca,'fontsize',16);
% % %                 title([mouse,'    -    Day=',num2str(day),'    -    Session=',num2str(s),'    -    MUnit=',num2str(munit),'    -    Distance=',num2str(floor(tot_position(end))),'cm' ]);
% % %                 %                 subplot(3,1,2)
% % %                 %                 plot((1:size(lap_position_d,1)-1)./8,speed_d,'k','linewidth',2)
% % %                 %                 ylabel({'Speed';'(cm/s)'});
% % %                 %                 set(gca,'fontsize',30)
% % %                 % %                 axis([620 780 -2 10])
% % %                 %                 yticks([0 10]);
% % %                 %                 xticks([]);
% % %                 %                 box off
% % %                 subplot(2,1,2)
% % %                 plot((1:size(lap_position_d,1))./8,lap_position_d,'k','linewidth',2)
% % %                 ylabel({'Position';'on wheel (cm)'});
% % %                 axis([-inf inf 1 inf]);
% % %                 %                 xticks([620 700 780])
% % %                 yticks([0 50]);
% % %                 box off
% % %                 xlabel('Time (s)');
% % %                 set(gca,'fontsize',16);
% % %                 
% % %                 saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July 2021\Raster plots and behavior\',mouse,'_Day',num2str(day),'_Session',num2str(s),'.svg'));
% % %                 saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July 2021\Raster plots and behavior\',mouse,'_Day',num2str(day),'_Session',num2str(s),'.fig'));
% % %                 
%                 
                distance_ran_V1(w)= floor(tot_position(end));
                mean_speed_V1(w)=mean(speed_d2);
                std_speed_V1(w)=std(speed_d2);
                mean_ac_V1(w)=mean(ac_d);
                std_ac_V1(w)=std(ac_d);
                ws_binary_V1(w)=big_table(row_w,6);
                ws_entropy_V1(w)=big_table(row_w,11);
                
                close all
                
                clear ac_d alfa cells_d lap_index_d lap_position_d N sortind_descend sorting sorting_0 sorting_ascend speed speed_s speed_d2 spikes_d ...
                    spikes_d_osf spikes_sh time_s simt_s_d times_tracking timestam tot_position tot_position_d times_s times_s_d spikes_d_s speed_d ...
                    motion
            end
        end
    end
end

%% Figures

figure
scatter(distance_ran_V1./100,ws_entropy_V1,50,'filled')
alpha 0.5
ylabel('Waveness');
xlabel('Total distance on wheel (m)');
hold on 
scatter(nonzeros(distance_ran_PaS)./100,ws_entropy_PaS,50,'filled');
scatter(nonzeros(distance_ran_MEC)./100,ws_entropy,50,'filled');
legend('V1','PaS','MEC');
set(gca,'fontsize',18);

distance_ran_PaS_p=distance_ran_PaS(find(distance_ran_PaS>0));
distance_ran_MEC_p=distance_ran_MEC(find(distance_ran_MEC>0));

[rho_d_v1,p_d_v1]=corr((distance_ran_V1./100)',ws_entropy_V1');
[rho_d_pas,p_d_pas]=corr((distance_ran_PaS_p./100)',ws_entropy_PaS');
[rho_d_mec,p_d_mec]=corr((distance_ran_MEC_p./100)',ws_entropy');

figure
scatter(mean_speed_V1./100,ws_entropy_V1,50,'filled')
alpha 0.5
ylabel('Waveness');
xlabel('Mean speed (cm/s)');
hold on 
scatter(nonzeros(mean_speed_PaS)./100,ws_entropy_PaS,50,'filled');
scatter(nonzeros(mean_speed)./100,ws_entropy,50,'filled');
legend('V1','PaS','MEC');
set(gca,'fontsize',18);

distance_speed_PaS_p=mean_speed_PaS(find(mean_speed_PaS>0));
distance_speed_MEC_p=mean_speed(find(mean_speed>0));

[rho_s_v1,p_s_v1]=corr((mean_speed_V1./100)',ws_entropy_V1');
[rho_s_pas,p_s_pas]=corr((distance_speed_PaS_p./100)',ws_entropy_PaS');
[rho_s_mec,p_s_mec]=corr((distance_speed_MEC_p./100)',ws_entropy');

figure
scatter(mean_ac_V1./100,ws_entropy_V1,50,'filled')
alpha 0.5
ylabel('Waveness');
xlabel('Mean acceleration (cm/s2)');
hold on 
scatter(nonzeros(mean_ac_PaS)./100,ws_entropy_PaS,50,'filled');
scatter(nonzeros(mean_ac)./100,ws_entropy,50,'filled');
legend('V1','PaS','MEC');
set(gca,'fontsize',18);

distance_ac_PaS_p=mean_ac_PaS(find(abs(mean_ac_PaS)>0));
distance_ac_MEC_p=mean_ac(find(find(mean_ac)>0));

[rho_a_v1,p_a_v1]=corr((mean_ac_V1./100)',ws_entropy_V1');
[rho_a_pas,p_a_pas]=corr((distance_ac_PaS_p./100)',ws_entropy_PaS');
[rho_a_mec,p_a_mec]=corr((distance_ac_MEC_p./100)',ws_entropy');

