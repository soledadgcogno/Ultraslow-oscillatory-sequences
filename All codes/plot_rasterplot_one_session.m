% Sort neurons according to angle between PC1 and PC2. Plot rasterplot.
% Cell 1: Data.
% Cell 2: Shuffled data.
% Cell 3: Data and running epochs
% Cell 4: Raster plot for a V1 session

%%
clear all
rec_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
spath = 'C:\Users\xscogno\Dropbox\SfN2019\Development\Rasterplots\L9M1\';

mouse='L9M4';  %Name of the animal
load([rec_path,strcat('recording_dates_',mouse,'.mat')]);

%% Raster plot with sorted data

for day=17%:dates.daysnum
    
    for s=1:dates.sesnum(day)
        munit=dates.ses{day}(s);
        
        file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (exist(file_name) == 2)
            disp(day)
            
            load(file_name,'-mat');
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            if(N>150)                
                [sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
                
                fig=figure;
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                hold on
                for i=1:size(spikes_d,1)
%                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_descend(i),:),5,'k','filled')
                    scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_ascend(i),:),5,'k','filled')

                    alpha 0.3
                end
                axis([-inf inf 1 inf]);
                %title([mouse,' Day',num2str(day)]);
                xlabel('Time (s)');
                ylabel('Neurons #');
                set(gca,'fontsize',18);
                yticks([100 400])
                
            end
            %             clear fig
            %             close all
        end
%         clear spikes_d_aux spikes_d_2  mean_act  smooth_act c  spikes_d_aux_n mean_act smooth_act P index ...
%             ind coeff spikes_d coeff a1 b1 spikes_d
    end
end


%% Raster plot with sorted data
clear all
rec_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
spath = 'C:\Users\xscogno\Dropbox\SfN2019\Development\Rasterplots\L9M1\';
sf=7.73;
mouse='L8M2';  %Name of the animal
load([rec_path,strcat('recording_dates_',mouse,'.mat')]);

for day=16%:dates.daysnum
    
    for s=1%:dates.sesnum(day)
        munit=dates.ses{day}(s);
        
        file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (exist(file_name) == 2)
            disp(day)
            
            load(file_name,'-mat');
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            if(N>150)                
                [sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
                
                fig=figure;
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                hold on
                for i=1:size(spikes_d,1)
%                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_descend(i),:),5,'k','filled')
                    scatter((1:size(spikes_d,2))./sf,i*spikes_d(sorting_descend(i),:),5,'k','filled')

                    alpha 0.3
                end
                axis([-inf inf 1 inf]);
                %title([mouse,' Day',num2str(day)]);
                xlabel('Time (s)');
                ylabel('Neurons #');
                set(gca,'fontsize',18);
                yticks([100 400])
                
            end
            %             clear fig
            %             close all
        end
%         clear spikes_d_aux spikes_d_2  mean_act  smooth_act c  spikes_d_aux_n mean_act smooth_act P index ...
%             ind coeff spikes_d coeff a1 b1 spikes_d
    end
end

%% Raster plot with sorted data

for day=17%:dates.daysnum
    
    for s=1:dates.sesnum(day)
        munit=dates.ses{day}(s);
        
        file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (exist(file_name) == 2)
            disp(day)
            
            load(file_name,'-mat');
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            if(N>150)                
                [sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
                
                fig=figure;
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                hold on
                for i=1:size(spikes_d,1)
%                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_descend(i),:),5,'k','filled')
                    scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_ascend(i),:),5,'k','filled')

                    alpha 0.3
                end
                axis([-inf inf 1 inf]);
                %title([mouse,' Day',num2str(day)]);
                xlabel('Time (s)');
                ylabel('Neurons #');
                set(gca,'fontsize',18);
                yticks([100 400])
                
            end
            %             clear fig
            %             close all
        end
%         clear spikes_d_aux spikes_d_2  mean_act  smooth_act c  spikes_d_aux_n mean_act smooth_act P index ...
%             ind coeff spikes_d coeff a1 b1 spikes_d
    end
end


%% Raster plot after shuffling spike matrix

for day=17%:dates.daysnum
    
    for s=1:dates.sesnum(day)
        munit=dates.ses{day}(s);
        
        file_name=[dpath ['spikes_120ms_Do_THR2_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (exist(file_name) == 2)
            disp(day)
            
            load(file_name,'-mat');
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            if(N>150)                
                spikes_sh=shuffle(spikes_d')';
                [sorting_ascend,sortind_descend,sorting_0]=get_sorting(spikes_sh);
                
                fig=figure;
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                hold on
                for i=1:size(spikes_d,1)
                    scatter((1:size(spikes_d,2))./8,i*spikes_sh(sortind_descend(i),:),5,'k','filled')
                    alpha 0.2
                end
                axis([-inf inf 1 inf]);
                %title([mouse,' Day',num2str(day)]);
                xlabel('Time (s)');
                ylabel('Neurons #');
                set(gca,'fontsize',18);
                yticks([100 400])

                
            end
            %             clear fig
            %             close all
        end
        clear spikes_d_aux spikes_d_2  mean_act  smooth_act c  spikes_d_aux_n mean_act smooth_act P index ...
            ind coeff spikes_d coeff a1 b1 spikes_d
    end
end

%% Raster plot with random sorting

for day=17%:dates.daysnum
    
    for s=1:dates.sesnum(day)
        munit=dates.ses{day}(s);
        
        file_name=[dpath ['spikes_120ms_Do_THR2_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (exist(file_name) == 2)
            disp(day)
            
            load(file_name,'-mat');
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            if(N>150)                
                [sortind_descend]=randperm(N);
                
                fig=figure;
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                hold on
                for i=1:size(spikes_d,1)
                    scatter((1:size(spikes_d,2))./8,i*spikes_d(sortind_descend(i),:),5,'k','filled')
                    alpha 0.3
                end
                axis([-inf inf 1 inf]);
                %title([mouse,' Day',num2str(day)]);
                xlabel('Time (s)');
                ylabel('Neurons #');
                set(gca,'fontsize',18);
                yticks([100 400])
 
            end
            %             clear fig
            %             close all
        end
        clear spikes_d_aux spikes_d_2  mean_act  smooth_act c  spikes_d_aux_n mean_act smooth_act P index ...
            ind coeff spikes_d coeff a1 b1 spikes_d
    end
end

%% Data and running epochs

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

speed_smoothing_kernel=5;
threshold_still=2; %if speed<threshold_still then speed=0
R=8.5145; %Radius of the wheel in cm
m=2;
for day=19
    
    for s=1:dates.sesnum(day)
        disp(s)
        munit=dates.ses{day}(s);
        file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        file_name_osf=[dpath ['spikes_30ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
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
                if length(num2str(dates.days(day)))==7
                    st=num2str(dates.days(day));
                    st_new=['0',st];
                    file_name_beh=[dbeh_path,mouse,'\Flavio_2P_',mouse,'_',st_new,'_MUnit_',num2str(munit),'_TRACKING.csv'];
                else
                    file_name_beh=[dbeh_path,mouse,'\Flavio_2P_',mouse,'_',num2str(dates.days(day)),'_MUnit_',num2str(munit),'_TRACKING.csv'];
                end
                
                alfa=table2array(readtable(file_name_beh)); %Table with info about position and time
                timestam=alfa(:,1);
                
                if m>4
                    %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
                    %Order: Timestamps	Clock	Position	Lap_position	Lap_index	Motor
                    tot_position=alfa(:,3)/10; %in cm
                    lap_position=alfa(:,4)/10; %in cm
                    lap_index=alfa(:,5);
                    motor=alfa(:,6);
                    
                    %%%%%%%%%%%%%%%%%%%%%% Quantities from the table
                    R=max(lap_position)/(2*pi); %in mm
                    speed=diff(tot_position)./diff(timestam); %cm/s % No smoothing
                    ac=diff(speed)./diff(timestam(1:end-1));
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
                    %         lap_position_d=interp1(timestam,lap_position,times_tracking);
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
                    ac_d(end+1)=ac_d(end);
                end
                
                motion=ones(1,T);
                motion(speed_d2==0)=0;
                
                %Rasterplot
                figure;
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                hold on
                for i=1:N
                    scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting(i),:),4,'k','filled')
                    alpha 0.3
                end
                box on
                h=area((1:size(spikes_d,2))./8,N.*motion,'LineStyle','-');
                h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
                h(1).FaceColor = [0.2 0.8 0.8];
                h(1).FaceAlpha = 0.4;
                %title([mouse,' Day',num2str(day)]);
                axis([-inf inf 1 inf]);
                yticks([100 400])
                xlabel('Time (s)');
                ylabel('Neurons #');
                set(gca,'fontsize',20);
                box off
                
                %Extended version of rasterplot
                figure;
                subplot(3,1,1)
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                hold on
                for i=1:N
                    scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting(i),:),4,'k','filled')
                    %                     alpha 0.15
                    alpha 0.6
                end
                box on
                h=area((1:size(spikes_d,2))./8,N.*motion,'LineStyle','-');
                h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
                h(1).FaceColor = [0.2 0.8 0.8];
                h(1).FaceAlpha = 0.4;
                %title([mouse,' Day',num2str(day)]);
                axis([620 780 1 inf]);
                xticks([]);
                %                 xlabel('Time (s)');
                ylabel('Neurons #');
                yticks([100 400])
                set(gca,'fontsize',30);
                subplot(3,1,2)
                plot((1:size(lap_position_d,1)-1)./8,speed_d,'k','linewidth',2)
                ylabel({'Speed';'(cm/s)'});
                set(gca,'fontsize',30)
                axis([620 780 -2 10])
                yticks([0 10]);
                xticks([]);
                box off
                subplot(3,1,3)
                plot((1:size(lap_position_d,1))./8,lap_position_d,'k','linewidth',2)
                ylabel({'Position';'on wheel (cm)'});
                axis([620 780 -inf inf])
                xticks([620 700 780])
                yticks([0 50]);
                box off
                xlabel('Time (s)');
                set(gca,'fontsize',30)
                
            end
        end
    end
end

%% Rasterplot for DFF signal 

file_name_dff=[dpath ['DFF_120ms_Do_SNRH','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_dff,'-mat');
dff=signal_dff;

for i=1:N
    dff_z(i,:)=zscore(dff(i,:));
end

for i=1:N
    aux=find(dff_z(i,:)<1);
    dff_z_t(i,:)=dff_z(i,:);
    dff_z_t(i,aux)=0;
    clear aux
end

fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
imagesc(dff_z((sorting_descend),:));
colormap jet
axis([-inf inf 1 inf]);
%title([mouse,' Day',num2str(day)]);
xlabel('Time (s)');
ylabel('Neurons #');
set(gca,'fontsize',18);
yticks([84 384])
yticklabels([400 100])
caxis([0 4])
set(gca,'fontsize',24);
xticks([8*500 8*1000 8*1500 8*2000 8*2500 8*3000])
xticklabels([500 1000 1500 2000 2500 3000])
colorbar 

fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
imagesc(dff_z_t((sorting_ascend),:));
colormap jet
caxis([0 5])
axis([-inf inf 1 inf]);
xticks([8*200 8*400 8*600 8*800 8*1000 8*1200 8*1400 8*1600])
yticklabels([200 400 600 800 1000 1200 1400 1600])
%title([mouse,' Day',num2str(day)]);
xlabel('Time (s)');
ylabel('Neurons #');
set(gca,'fontsize',18);
% yticks([20 420])
% yticklabels([500 100])
% yticks([120 420])
yticks([84 384])
yticklabels([400 100])
colorbar


%% Rasterplot for one V1 session

% clear all
rec_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
spath = 'C:\Users\xscogno\Dropbox\SfN2019\Development\Rasterplots\L9M1\';

mouse='92229';
load([rec_path,strcat('recording_dates_',mouse,'.mat')]);

% Raster plot with sorted data

for day=7%:dates.daysnum
    
    for s=1:dates.sesnum(day)
        munit=dates.ses{day}(s);
        
        file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (exist(file_name) == 2)
            disp(day)
            
            load(file_name,'-mat');
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            if(N>150)                
                [sorting_ascend,sortind_descend,sorting_0]=get_sorting(spikes_d);
                
                fig=figure;
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                hold on
                for i=1:size(spikes_d,1)
                    scatter((1:size(spikes_d,2))./8,i*spikes_d(sortind_descend(i),:),5,'k','filled')
                    alpha 0.3
                end
                axis([-inf inf 1 inf]);
                %title([mouse,' Day',num2str(day)]);
                axis([-inf inf 1 inf])
                ylabel('Neurons #');
                xlabel('Time (s)');
                set(gca,'fontsize',18)
                yticks([100 200])
                
            end
            %             clear fig
            %             close all
        end
        clear spikes_d_aux spikes_d_2  mean_act  smooth_act c  spikes_d_aux_n mean_act smooth_act P index ...
            ind coeff spikes_d coeff a1 b1 spikes_d
    end
end

%% Raster plot after shuffling spike matrix

for day=7%:dates.daysnum
    
    for s=1:dates.sesnum(day)
        munit=dates.ses{day}(s);
        
        file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (exist(file_name) == 2)
            disp(day)
            
            load(file_name,'-mat');
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            if(N>150)                
                spikes_sh=shuffle(spikes_d')';
                [sorting_ascend,sortind_descend,sorting_0]=get_sorting(spikes_sh);
                
                fig=figure;
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                hold on
                for i=1:size(spikes_d,1)
                    scatter((1:size(spikes_d,2))./8,i*spikes_sh(sortind_descend(i),:),5,'k','filled')
                    alpha 0.2
                end
                axis([-inf 13926/8 1 inf]);
                %title([mouse,' Day',num2str(day)]);
                xlabel('Time (s)');
                ylabel('Neurons #');
                set(gca,'fontsize',18);
                yticks([100 200])

                
            end
            %             clear fig
            %             close all
        end
        clear spikes_d_aux spikes_d_2  mean_act  smooth_act c  spikes_d_aux_n mean_act smooth_act P index ...
            ind coeff spikes_d coeff a1 b1 spikes_d
    end
end


%% Rasterplot for one MEC medial session

clear all
rec_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
spath = 'C:\Users\xscogno\Dropbox\SfN2019\Development\Rasterplots\L9M1\';

% mouse='L8M3';
mouse='L8M4';

load([rec_path,strcat('recording_dates_',mouse,'.mat')]);

% Raster plot with sorted data

for day=17%19:dates.daysnum
    
    for s=2%:dates.sesnum(day)
        munit=dates.ses{day}(s);
        
        file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (exist(file_name) == 2)
            disp(day)
            
            load(file_name,'-mat');
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            if(N>150)                
                [sorting_ascend,sortind_descend,sorting_0]=get_sorting(spikes_d);
                
                fig=figure;
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                hold on
                for i=1:size(spikes_d,1)
                    scatter((1:size(spikes_d,2))./8,i*spikes_d(sortind_descend(i),:),5,'k','filled')
                    alpha 0.3
                end
                axis([-inf inf 1 inf]);
                %title([mouse,' Day',num2str(day)]);
                axis([-inf inf 1 inf])
                ylabel('Neurons #');
                xlabel('Time (s)');
                set(gca,'fontsize',18)
                yticks([100 400])
                
            end
            %             clear fig
            %             close all
        end
        clear spikes_d_aux spikes_d_2  mean_act  smooth_act c  spikes_d_aux_n mean_act smooth_act P index ...
            ind coeff spikes_d coeff a1 b1 spikes_d
    end
end
%% Plot raster plot and phase of the wave


for day=17%:dates.daysnum
    
    for s=1:dates.sesnum(day)
        munit=dates.ses{day}(s);        
        file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        
        if (exist(file_name) == 2)
            disp(day)
            
            load(file_name,'-mat');
            spikes_d=full(spikes_d_s);
            [N,T]=size(spikes_d);
            
            if(N>150)                
                [sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
                [~,score] = pca(spikes_d');
                phase=atan2(smooth(score(:,2),80),smooth(score(:,1),80));

                
                window=[2690:3150];
                window_2=window*8;
                fig=figure;
                subplot(2,1,1)
                plot(window_2/8,phase(window_2),'K','linewidth',2);
                ylim([-3.14 3.14]);
                yticks([-3.14 0 3.14]);
                xlim([-inf inf])
                xticks([])
                ylabel('Phase (rad)');
                set(gca,'fontsize',16);
                box off
                subplot(2,1,2)
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.4 0.4]);
                hold on
                for i=1:size(spikes_d,1)
                    %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_descend(i),:),5,'k','filled')
                    scatter(window_2./8,i*spikes_d(sorting_ascend(i),window_2),5,'k','filled');
                    
                    alpha 0.3
                end
                axis([-inf inf 1 inf]);
                %title([mouse,' Day',num2str(day)]);
                xlabel('Time (s)');
                ylabel('Neurons #');
                set(gca,'fontsize',18);
                yticks([100 400])
                box off
                
                window=[2690:3150];
                window_2=window*8;
                phase_d=discretize(phase,N);
                fig=figure;
                colororder({'k','r'});
                yyaxis left
                hold on
                for i=1:size(spikes_d,1)
                    %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_descend(i),:),5,'k','filled')
                    scatter(window_2./8,i*spikes_d(sorting_ascend(i),window_2),5,'k','filled');
                    
                    alpha 0.4
                end
                ylim([2 N])
                xlim([-inf inf])
                ylabel('Neurons #');
                yticks([100 400])
                set(gca,'fontsize',16)
                xlabel('Time (s)');
                yyaxis right
                plot(window_2/8,phase(window_2),'r','linewidth',2);
                ylim([-3.14 3.14]);
                yticks([-3.14 0 3.14]);
                yticklabels({'-\pi','0','\pi'});
                ylabel('Phase (rad)');
                xticks([2800 3000]);
                set(gca,'fontsize',16)
                
                window=[2690:3150];
                window_2=window*8;
                phase_d=discretize(phase,N);
                fig=figure;
                colororder({'k','r'});
                yyaxis left
                hold on
                for i=1:size(spikes_d,1)
                    %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_descend(i),:),5,'k','filled')
                    scatter(window_2./8,i*spikes_d(sorting_ascend(i),window_2),5,'k','filled');
                    
                    alpha 0.4
                end                
                ylim([2 N])
                xlim([-inf inf])
                ylabel('Neurons #');
                yticks([100 400])
                set(gca,'fontsize',16)
                xlabel('Time (s)');
                yyaxis right
                plot(window_2/8,phase(window_2),'r','linewidth',2);
                ylim([-3.14 3.14]);
                yticks([-3.14 0 3.14]);
                yticklabels({'-\pi','0','\pi'});
                ylabel('Phase (rad)');
                xticks([2800 3000]);
                set(gca,'fontsize',16)
                
                
                
            end
            %             clear fig
            %             close all
        end
        clear spikes_d_aux spikes_d_2  mean_act  smooth_act c  spikes_d_aux_n mean_act smooth_act P index ...
            ind coeff spikes_d coeff a1 b1 spikes_d
    end
end

%% Plot raster plot and phase of the cycle, plus individual cycles

mouse='L9M4';
day=17;
s=1;
munit=0;
clus=10;
disc_phase=10;

file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name,'-mat');
spikes_d=full(spikes_d_s);
[N,T]=size(spikes_d);
num_clus_discr=10;
make_fig=1;
dt=117;
sf=7.73;

[table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes_d);

for i=1:N
      FRp(i,:)=smoothdata(spikes_d(i,:),'gaussian',1*dt); %smooth in the dt chosen for each session
end

[~,sorting_FRp,~] = get_sorting(FRp); %gets sorting based on smoothed spike matrix

[~,scoret,~] = pca(FRp');
phase_r=(atan2(scoret(:,2),scoret(:,1)));
phase_f=(atan2(smooth(scoret(:,2),floor(1*dt)),smooth(scoret(:,1),floor(1*dt))));
phase=phase_f;
sorting_FRp=flip(sorting_FRp);
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
    hold on
    for i=1:size(spikes_d,1)
        scatter((1:size(spikes_d,2))./sf,i*spikes_d(sorting_FRp(i),:),5,'k','filled')
        alpha 0.2
    end
    axis([-inf inf 1 inf]);
    ylabel('Neurons #');
    set(gca,'fontsize',16);
    xticks([]);
    yticks([100 400]);

    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
    plot((1:size(spikes_d,2))./sf,phase,'k','linewidth',1.5);
    hold on
    for i=1:size(table_u,1)
        plot((table_u(i,1):table_u(i,2))./sf,phase(table_u(i,1):table_u(i,2)),'linewidth',2,'color','c');
    end
    pbaspect([24 1 5])
    axis([-inf inf -inf inf])
    yticks([-3.14,0,3.14])
    yticklabels({'-\pi','0','\pi'})
    set(gca,'fontsize',16)
    xlabel('Time (s)');
    ylabel({'Phase of the';'cycle (rad)'})
    box off

