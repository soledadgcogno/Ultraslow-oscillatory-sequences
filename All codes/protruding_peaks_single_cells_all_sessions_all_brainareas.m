%% Frequency of single cell oscillations - All sessions

clear all
close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

[big_table, waves] = get_big_table();

%Params for analysis
downsampling_factor=4;
fs_120 = 7.73;
fs_30 = 30.95;
fs = fs_120/downsampling_factor;
maxlag=floor(560*fs_120);
N_sh=200;
count_c=0;
epoch_length=20;%1/8; %in seconds
win_size=ceil(fs_120)*epoch_length;


%% Wave identification - Wave sessions in MEC

downsampling_factor=4;
fs_120 = 7.73;
fs_30 = 30.95;
fs = fs_120/downsampling_factor;
maxlag=floor(560*fs_120);%240;%250 - 520;
N_sh=200;
count_c=0;

epoch_length=1;%1/8; %in seconds
win_size=ceil(fs_120)*epoch_length;

%Identifying waves
for w=1:length(waves)
    disp(w);
    row_w=waves(w);
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);
    munit=big_table(waves(w),5);
    dt=floor(big_table(waves(w),8));
    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_spk,'-mat');
    spikes=full(spikes_d_s);

    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);
    number_of_cells(w)=N;
    number_of_sequences(w)=size(table_u,1);

    %Prepares new spike matrix by keeping time bins with sequences only
    spikes_w=[];
    for i=1:size(table_u,1)
        spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
    end
    full_length(w)=size(spikes,2);
    wave_length(w)=size(spikes_w,2);

    %Calculates PSDs  
    FR=spikes;
    count=0;

    for i=1:N
        count_c=count_c+1;
       
        [Powerspec2, Powerfreq2,peak_freq,peak_psd,Powerspec2_sh,alfa] = get_PSD_significance(spikes_w(i,:), maxlag,fs_120,win_size,N_sh);

        clear aux; aux=find(Powerfreq2==peak_freq);
        thr=prctile(Powerspec2_sh,95,1);
        cell_ind_win1(count_c)=0;
        if (Powerspec2(aux)>thr(aux)) cell_ind_win1(count_c)=1;    end 
      

        clear  Powerspec2 P2 P1  Y Powerfreq2 signal laf alfa val thr spk Powerspec2_sh Powerfreq2_sh spk_d
    end

    clear table_u spikes dt spikes_d_s row_w mat_b spikes_W Ens_b sorting_before FR alfa val lags ...
        spikes_w table_u count N Powerspec2 Powerfreq2 spk cells_d
end

%% Visual cortex sessions

%Identify Visual Cortex sessions
V1_ses=find(big_table(:,10)<0); 
adults=find(big_table(:,3)>15); 
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
table_cells=[];
%Identify Visual Cortex sessions

for w=1:length(V1_sessions)    
    row_w=V1_sessions(w);
    disp(w)    
%     count=count+1;
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
    ws=big_table(row_w,6);
    ws_ent=big_table(row_w,11); 
    clus=10;
    disc_phase=10;    
    % Load files    
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];    
    if (isfield(dates,'folder_name')==1)
        day_index=find (dates.actual_day==day);
        file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day_index),'_MUnit',num2str(munit),'.mat']];
    else
        file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    end
    load(file_name_spk,'-mat'); %Spike times
    spikes_d=full(spikes_d_s);
    [N,T]=size(spikes_d);

    for i=1:N
        count_c=count_c+1;
        table_cells(count_c,1)=mouse_i;
        table_cells(count_c,2)=s;
        [Powerspec2, Powerfreq2,peak_freq,peak_psd,Powerspec2_sh,alfa] = get_PSD_significance(spikes_d(i,:), maxlag,fs_120,win_size,N_sh);
        clear aux; aux=find(Powerfreq2==peak_freq);
        thr95=prctile(Powerspec2_sh,95,1);
        thr99=prctile(Powerspec2_sh,99,1);

        cell_ind_95(count_c)=0;
        cell_ind_99(count_c)=0;
        cell_ind_win1(count_c)=0;

        if (Powerspec2(aux)>thr95(aux)) cell_ind_95(count_c)=1;    end
        if (Powerspec2(aux)>thr99(aux)) cell_ind_99(count_c)=1;    end
        if (peak_freq>0.1) cell_ind_win1(count_c)=1;    end

        clear  Powerspec2 P2 P1  Y Powerfreq2 signal laf alfa val thr spk Powerspec2_sh Powerfreq2_sh spk_d
    end

     clear spikes_d_s spikes_d sorting_descend sorting_ascend sorting_0 N mouse cells_d
     clear table_u spikes dt spikes_d_s row_w mat_b spikes_W Ens_b sorting_before FR alfa val lags ...
        spikes_w table_u count N Powerspec2 Powerfreq2 spk cells_d
end

% epoch_length=20;
% win_size=ceil(fs_120)*epoch_length;
