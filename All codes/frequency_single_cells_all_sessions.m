%% Frequency of single cell oscillations - All sessions

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

clus=10;
fs_120=7.73;

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
            %disp(s)
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

% Wave sessions

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


%% Wave identification

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


    downsampling_factor=4;
    fs_120 = 7.73;
    fs_30 = 30.95;
    fs = fs_120/downsampling_factor;
    maxlag=floor(560*fs_120);%240;%250 - 520;

    FR=spikes;
    count=0;

    for i=1:N
        [alfa(i,:),lags]=xcorr(spikes_w(i,:),spikes_w(i,:),maxlag,'coeff');
        [val(i,:),~]=zscore(alfa(i,:)); %Check whether I need to zscore
        signal=alfa(i,:);

        clear Powerspec2 Powerfreq2

        %Calculate spectogram using pwelch method
        [Powerspec2,Powerfreq2] = doPwelch(signal,fs_120,2*4096);

        %Check if the cell is oscillatory using the pwelch
        %method
        [peak_freq{w}(i),peak_psd{w}(i),quality{w}(i)]=check_peak_quality_3c_f(Powerspec2,Powerfreq2,0.04);

        if quality{w}(i)==1
            count=count+1;
            cells_osc{w}(count)=i;
        else     
        end
        clear  Powerspec2 P2 P1  Y Powerfreq2 signal
    end

      % Frequency of the global oscillation
%     full_length_oscillation=sum(table_u(:,2)-table_u(:,1))/sf;
    full_length_oscillation2=size(spikes_w,2)/fs_120;

%     freq1(w)=size(table_u,1)/full_length_oscillation;
    freq(w)=size(table_u,1)/full_length_oscillation2;

    rel_freq{w}=peak_freq{w}./freq(w);


    clear table_u spikes dt spikes_d_s row_w mat_b spikes_W Ens_b sorting_before FR alfa val lags ...
        spikes_w table_u count N Powerspec2 Powerfreq2

end

mouse_idx=[1,1,1,2,2,2,3,3,3,3,4,4,4,4,4];

%% figures

%Example session
w=8;
number_cells=size(peak_freq{w},2);
figure
[h,edges]=histcounts(peak_freq{w},[0:0.001:0.1]);
bar(edges(1:end-1),h);

cells_less_0p01=length(find(peak_freq{w}<0.01))/number_of_cells(w);


for w=1:15
    min_freq(w)=min(peak_freq{w});
    max_freq(w)=max(peak_freq{w});
    median_freq(w)=median(peak_freq{w});

    prob_freq(w,:)=histcounts(peak_freq{w},[0:0.001:0.15],'Normalization','probability');
    oscillatory_fraction_full(w)=size(cells_osc{w},2)/size(peak_freq{w},2);
    [b,max_prob(w)]=max(prob_freq(w,:));

end

total_cells=sum(number_of_cells);
freqs=(0:0.001:0.15)+(0.001/2);
max_prob=freqs(max_prob);

rate_waves_tot=number_of_sequences./(full_length./fs_120);
rate_waves_wave=number_of_sequences./(wave_length./fs_120);

figure
scatter(max_prob,rate_waves_tot)
hold on
refline(1,0)

[p,h]=ttest(max_prob-rate_waves_tot);

figure
scatter(max_prob,rate_waves_wave)
hold on
refline(1,0)

w=8;
number_cells=size(peak_freq{w},2);
figure
[h,edges]=histcounts(rel_freq{w},[0:0.5:20]);
bar(edges(1:end-1),h);

%% final figures All sessions

rel_freqs=[];
for i=1:15
    rel_freqs=[rel_freqs;rel_freq{i}'];
end

%all sessions
figure
[h,edges]=histcounts(rel_freqs,[-0.25:0.5:20]);
bar(edges(1:end-1)+0.25,h,'FaceColor','k');
set(gca,'fontsize',16,'YColor','k','XColor','k');
axis([-0.5 5 0 3000]);
xticks([0 1 2 3 4])
box off
ylabel('Counts');
xlabel({'Single cell frequency relative to';'ultraslow population oscillation frequency'});
hold on
xline(prctile(rel_freqs,75),'r--','linewidth',2);
% text(prctile(rel_freqs,75),3000,'  75','Color','red')
% xline(prctile(rel_freqs,50),'r--','linewidth',2);
% text(prctile(rel_freqs,75),3000,'  50','Color','red')
xline(prctile(rel_freqs,25),'r--','linewidth',2);
% text(prctile(rel_freqs,75),3000,'  25','Color','red')

figure
[h,edges]=histcounts(rel_freqs,[-0.125:0.25:20]);
bar(edges(1:end-1)+0.125,h,'FaceColor','k');
set(gca,'fontsize',16,'YColor','k','XColor','k');
axis([-0.5 5 0 1700]);
xticks([0 1 2 3 4])
box off
ylabel('Counts');
xlabel({'Single cell frequency relative to';'ultraslow population oscillation frequency'});
hold on
xline(prctile(rel_freqs,75),'r--','linewidth',2);
% text(prctile(rel_freqs,75),3000,'  75','Color','red')
% xline(prctile(rel_freqs,50),'r--','linewidth',2);
% text(prctile(rel_freqs,75),3000,'  50','Color','red')
xline(prctile(rel_freqs,25),'r--','linewidth',2);
% text(prctile(rel_freqs,75),3000,'  25','Color','red')



lb=find(rel_freqs>0.5);
ub=find(rel_freqs<1.5);
cells_with_wimilar_freq=length(intersect(lb,ub));

%one session
% figure
% [h,edges]=histcounts(rel_freq{w},[-0.25:0.5:20]);
% bar(edges(1:end-1)+0.25,h,'FaceColor','k');
% set(gca,'fontsize',16,'YColor','k','XColor','k');
% axis([0 5 0 350]);
% xticks([0 1 2 3 4])
% box off
% ylabel('Counts');
% xlabel({'Single cell frequency relative to';'population ultraslow oscillation frequency'});

