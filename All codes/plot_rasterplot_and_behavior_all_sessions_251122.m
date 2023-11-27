clear all
rec_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
spath = 'C:\Users\xscogno\Dropbox\SfN2019\Development\Rasterplots\L9M1\';

mouse='L11M2';  %Name of the animal
V1=0;
load([rec_path,strcat('recording_dates_',mouse,'.mat')]);
% Data and running epochs

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
sf=7.73; 

speed_smoothing_kernel=5;
threshold_still=2; %if speed<threshold_still then speed=0
R=8.5145; %Radius of the wheel in cm
m=5;
for day=1:dates.daysnum%18%16:19%:19
    
    for s=1:dates.sesnum(day)
        disp(s)
        munit=dates.ses{day}(s);
        file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        file_name_osf=[dpath ['spikes_30ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

%         file_name=[dpath ['spikes_120ms_Do_THR1_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
%         file_name_osf=[dpath ['spikes_30ms_Do_THR1_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

        if (exist(file_name) == 2)
            disp(day)
            
            load(file_name,'-mat');
            spikes_d=full(spikes_d_s);
            load(file_name_osf,'-mat');
            spikes_d_osf=full(spikes_d_s);
            
            [~,T]=size(spikes_d);
            [N,T_osf]=size(spikes_d_osf);


            [~,sorting,~]=get_sorting(spikes_d);
            
            if N>150
                
                disp(day)

                if V1==1
                    if length(num2str(dates.days(day)))==7
                        st=num2str(dates.days(day));
                        st_new=['0',st];
                        file_name_beh=[dbeh_path,mouse,'\Flavio_2P_V1_Spring 2020_',mouse,'_',st_new,'_MUnit_',num2str(munit),'_TRACKING.csv'];
                    else
                        file_name_beh=[dbeh_path,mouse,'\Flavio_2P_V1_Spring 2020_',mouse,'_',num2str(dates.days(day)),'_MUnit_',num2str(munit),'_TRACKING.csv'];
                    end
                else
                    if length(num2str(dates.days(day)))==7
                        st=num2str(dates.days(day));
                        st_new=['0',st];
                        file_name_beh=[dbeh_path,mouse,'\Flavio_2P_',mouse,'_',st_new,'_MUnit_',num2str(munit),'_TRACKING.csv'];
                    else
                        file_name_beh=[dbeh_path,mouse,'\Flavio_2P_',mouse,'_',num2str(dates.days(day)),'_MUnit_',num2str(munit),'_TRACKING.csv'];
                    end
                end
                
                alfa=table2array(readtable(file_name_beh)); %Table with info about position and time
                timestam=alfa(:,1);
                
                if m>4
                    %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
                    %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
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
% %                     ac_d(end+1)=ac_d(end);
%                     ac_d=diff(speed_d2)./diff(times_tracking(1:end));
%                     ac_d(end+1)=ac_d(end);
                    
                else
                    tot_position=alfa(:,2)/10; %in cm
                    
                    %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
                    
                    speed=diff(tot_position)./diff(timestam); %cm/s % Smoothing kernel of 2.3 seconds
                    
                    %%%%%%%%%%%%%%%%%%%%%% Imaging time steps and camera time steps
                    
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
%                     ac_d(end+1)=ac_d(end);
                end
                
% %                 motion=ones(1,T);
% %                 motion(speed_d2==0)=0;
% %                 
% %                 %Rasterplot
% %                 figure;
% %                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
% %                 hold on
% %                 for i=1:N
% %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting(i),:),4,'k','filled')
% %                     alpha 0.3
% %                 end
% %                 box on
% %                 h=area((1:size(spikes_d,2))./8,N.*motion,'LineStyle','-');
% %                 h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
% %                 h(1).FaceColor = [0.2 0.8 0.8];
% %                 h(1).FaceAlpha = 0.4;
% %                 %title([mouse,' Day',num2str(day)]);
% %                 axis([-inf inf 1 inf]);
% %                 yticks([100 400])
% %                 xlabel('Time (s)');
% %                 ylabel('Neurons #');
% %                 set(gca,'fontsize',20);
% %                 box off
%                sorting=flip(sorting);
                %Extended version of rasterplot
                figure;
                subplot(2,1,1)
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.7]);
                hold on
                for i=1:N
                    scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting(i),:),4,'k','filled')
                    %                     alpha 0.15
                    alpha 0.2
                end
                box on
%                 h=area((1:size(spikes_d,2))./8,N.*motion,'LineStyle','-');
%                 h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
%                 h(1).FaceColor = [0.2 0.8 0.8];
%                 h(1).FaceAlpha = 0.4;
                %title([mouse,' Day',num2str(day)]);
                axis([-inf inf 1 inf]);
                xticks([]);
                %                 xlabel('Time (s)');
                ylabel('Neurons #');
%                 yticks([100 400])
                set(gca,'fontsize',16);  
             %   title([mouse,'    -    Day=',num2str(day),'    -    Session=',num2str(s),'    -    MUnit=',num2str(munit),'    -    Distance=',num2str(floor(tot_position(end))),'cm' ]);
%                 subplot(3,1,2)
%                 plot((1:size(lap_position_d,1)-1)./8,speed_d,'k','linewidth',2)
%                 ylabel({'Speed';'(cm/s)'});
%                 set(gca,'fontsize',30)
% %                 axis([620 780 -2 10])
%                 yticks([0 10]);
%                 xticks([]);
%                 box off
                subplot(2,1,2)
                plot((1:size(speed_d2,1))./sf,speed_d2,'k','linewidth',2)
                ylabel({'Speed (cm/s)'});
                axis([-inf inf 1 inf]);
%                 xticks([620 700 780])
%                 yticks([0 20]);
                box off
                xlabel('Time (s)');
                set(gca,'fontsize',16);
                saveas(gcf,strcat('W:\soledad\rasterplots\calcium\',mouse,'\day_',num2str(day),'_ses',num2str(s),'.png'));
                close all
                
                clear ac_d alfa cells_d lap_index_d lap_position_d N sortind_descend sorting sorting_0 sorting_ascend speed speed_s speed_d2 spikes_d ...
                    spikes_d_osf spikes_sh time_s simt_s_d times_tracking timestam tot_position tot_position_d times_s times_s_d spikes_d_s speed_d ...
                motion lap_index lap_position motor motor_d angle angular_speed T T_osf 
            end
        end
    end
end