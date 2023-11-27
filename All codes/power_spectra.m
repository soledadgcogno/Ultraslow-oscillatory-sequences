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

% figpath='C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Raster Plots all MEC sessions\';
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

%% Calculates the spectograms for L9M4 Day 17

% For example session use w=8

% file_name_mvl='C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\MVL_all_sessions';
% load(file_name_mvl,'-mat');

for w=8
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
    
    
    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    file_name_spk_120=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_spk_30=[dpath ['spikes_30ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_spk_dff=[dpath ['DFF_120ms_Do_SNRH','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    
    aux=load(file_name_spk_120,'-mat');
    spikes_120=full(aux.spikes_d_s);
    clear aux
    aux=load(file_name_spk_30,'-mat');
    spikes_30=full(aux.spikes_d_s);
    clear aux
    aux=load(file_name_spk_dff,'-mat');
    spikes_dff=full(aux.signal_dff);
    
      
    fs_120 = 7.73;
    fs_30 = 30.95;
    %sampling frequency in tetrode data was 400Hz
    
    maxf = 0.5; %100; %Hz
    window_size = 4096;%1024;%80;%1024;
    overlap = window_size/2;%512;%512;%40;%512;
    fonts = 16;
    lines = 2;

    [~,sorting_w,~]=get_sorting(spikes_120);

    figure
    spy(spikes_120(sorting_w,:));
    pbaspect([8 1 1]);

    cells=[60,210,418];

    for i=cells
        x_30 = spikes_30(sorting_w(i),:);
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
        [Y,F,T,PSD] = spectrogram([x_30],window_size,floor(window_size*0.9),0:0.001:maxf,fs_30,'yaxis');
        TT = linspace(0,length(x_30)/fs_30,length(x_30));
        imagesc(TT,F,(abs(PSD).^0.25))
        set(gca,'YDir','normal'), colormap(jet)
        axis xy
        axis tight
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        axis([-inf inf 0 0.4])
        %     set(gca,'XTick',[0 1200 2400 3600],...
        %         'YTick',[0 2 4 6 8 10],'FontSize',fonts)
        set(gca,'fontsize',16)
        box on
        colorbar
        caxis([0 1])
    end

    
%     sum_cells=sum(spikes_30(sorting_w(1:10),:));
%     figure
%     [Y,F,T,PSD] = spectrogram([sum_cells],window_size/8,overlap/8,0:0.001:maxf,fs_30,'yaxis');
%     TT = linspace(0,length(x_30)/fs_30,length(x_30));
%     imagesc(TT,F,(abs(PSD).^0.25))
%     set(gca,'YDir','normal'), colormap(jet)
%     axis xy
%     axis tight
%     xlabel('Time [s]');
%     ylabel('Frequency [Hz]');
%     %     set(gca,'XTick',[0 1200 2400 3600],...
%     %         'YTick',[0 2 4 6 8 10],'FontSize',fonts)
%     set(gca,'fontsize',16)
%     box on
    
    %         figure
    %     plot(sum_cells)
    
    
   
%     figure
%     plot(x_30)
    
end


%% Calculates the spectograms for L8M2 Day 19
% 
% file_name_mvl='C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\MVL_all_sessions';
% load(file_name_mvl,'-mat');

for w=2
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
    file_name_spk_120=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_spk_30=[dpath ['spikes_30ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_spk_dff=[dpath ['DFF_120ms_Do_SNRH','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    
    aux=load(file_name_spk_120,'-mat');
    spikes_120=full(aux.spikes_d_s);
    clear aux
    aux=load(file_name_spk_30,'-mat');
    spikes_30=full(aux.spikes_d_s);
    clear aux
    aux=load(file_name_spk_dff,'-mat');
    spikes_dff=full(aux.signal_dff);
    
    
    %  load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    
    fs_120 = 7.73;
    fs_30 = 30.95;
    %sampling frequency in tetrode data was 400Hz
    
    maxf = 0.5; %100; %Hz
    window_size = 4096;%1024;%80;%1024;
    overlap = window_size/2;%512;%512;%40;%512;
    fonts = 16;
    lines = 2;

    [~,sorting_w,~]=get_sorting(spikes_120);

    figure
    spy(spikes_120(sorting_w,:));
    pbaspect([8 1 1]);
    

    
    [~,sorting_w,~]=get_sorting(spikes_120);
    x_30 = spikes_30(sorting_w(3),:);

    cells=[20,400];

    for i=400
        x_30 = spikes_30(sorting_w(i),:);
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
        [Y,F,T,PSD] = spectrogram([x_30],window_size/8,overlap/8,0:0.001:maxf,fs_30,'yaxis');
        TT = linspace(0,length(x_30)/fs_30,length(x_30));
        imagesc(TT,F,(abs(PSD).^0.25))
        set(gca,'YDir','normal'), colormap(jet)
        axis xy
        axis tight
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        axis([-inf inf 0 0.4])
        %     set(gca,'XTick',[0 1200 2400 3600],...
        %         'YTick',[0 2 4 6 8 10],'FontSize',fonts)
        set(gca,'fontsize',16)
        box on
        colorbar
        caxis([0 1])
    end



    
%     sum_cells=sum(spikes_30(sorting_w(1:10),:));
%     figure
%     [Y,F,T,PSD] = spectrogram([sum_cells],window_size/8,overlap/8,0:0.001:maxf,fs_30,'yaxis');
%     TT = linspace(0,length(x_30)/fs_30,length(x_30));
%     imagesc(TT,F,(abs(PSD).^0.25))
%     set(gca,'YDir','normal'), colormap(jet)
%     axis xy
%     axis tight
%     xlabel('Time [s]','FontSize',fonts);
%     ylabel('Frequency [Hz]','FontSize',fonts);
%     %     set(gca,'XTick',[0 1200 2400 3600],...
%     %         'YTick',[0 2 4 6 8 10],'FontSize',fonts)
%     set(gca,'fontsize',16)
%     box on
    
    %         figure
    %     plot(sum_cells)
    
    
   
%     figure
%     plot(x_30)
    
end

