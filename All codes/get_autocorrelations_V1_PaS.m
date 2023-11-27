% Autocorrelations for PaS example session

clear all
close all

rec_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

mouse='L8M4';
load([rec_path,strcat('recording_dates_',mouse,'.mat')]);
day=17;
s=2;
munit=dates.ses{day}(s);
file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name,'-mat');
spikes=full(spikes_d_s);
[N,T]=size(spikes);
sf=7.73;
sf_d=7.73/4;
T_PaS=T;

index_c=0;
index_t=0;
count=0;
clear period quality

% spikes=spikes(:,1:8000);

downsampling_factor=4;
fs_120 = 7.73;
fs_30 = 30.95;
fs = fs_120/downsampling_factor;
maxlag=floor(560*fs_120);%240;%250 - 520;

[N,T]=size(spikes);

FR=spikes;
count=0;
for i=1:N
    [alfa(i,:),lags]=xcorr(FR(i,:),FR(i,:),maxlag,'coeff');
    [val(i,:),~]=zscore(alfa(i,:)); %Check whether I need to zscore
    signal=alfa(i,:);

    clear Powerspec2 Powerfreq2

    %Calculate spectogram using pwelch method
    [Powerspec2,Powerfreq2] = doPwelch(signal,fs_120,2*4096);
    pwelch_fft(i,:)=Powerspec2;
    pwelch_fft_section(i,:)=Powerspec2(5:100);


    %Check if the cell is oscillatory using the pwelch
    %method
    [peak_freq(i),peak_psd(i),quality(i)]=check_peak_quality_3c_f(Powerspec2,Powerfreq2,0.04);
    %         period_fft(i)=1/peak_fft(i);
    period_ps(i)=1/peak_freq(i);
    if quality(i)==1
        count=count+1;
        cells_osc(count)=i;
%         if count<50
%             figure
%             plot(Powerfreq2(1:100),Powerspec2(1:100))
%             title(i)
%         end
    else
%         figure
%         plot(Powerfreq2(1:100),Powerspec2(1:100))
%         title(i)
    end
    clear  Powerspec2 P2 P1  Y
end

[val_sorting_ps2,sorting_ps2]=sort(peak_psd,'descend');
[val_sorting_freq2,sorting_freq2]=sort(peak_freq,'descend');
[sorting_ascend,sortind_descend,sorting_0] = get_sorting(spikes);

figure
spy(spikes(sortind_descend,:))
pbaspect([8 1 1])

figure
histogram(peak_freq,20)

length(find(peak_freq>0.008))

figure
[h,edges]=histcounts(peak_freq,[0:0.001:0.1]);

bar(edges(1:end-1),h);

% Figures autocorrelograms for all cells


%PCA sorting
figure
imagesc(val(sorting_ascend,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 400])
axis square
colorbar
caxis([0 1.5])

%Frequency sorting
figure
imagesc(val(sorting_freq2,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 400])
axis square
colorbar
title('PaS - Freq sorting')


%Power sorting
figure
imagesc(val(sorting_ps2,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 400])
axis square
colorbar
title('PaS - Power sorting')

%PCA sorting
figure
imagesc(val(sortind_descend,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 400])
axis square
colorbar
title('PaS - PCA sorting')

%% Example cells from PaS

for j=[2,5,10,18]%1:20
    i=sorting_ps2(j);
    figure
    subplot(1,2,1)
    plot((alfa(i,:)),'k','linewidth',3.5)
    xticks([1 maxlag maxlag*2])
    xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
    xlabel('Time lag (s)');
    ylabel('Autocorrelation');
    axis([-inf inf 0 0.5])
    set(gca,'fontsize',30);
    % h = gca; h.YAxis.Visible = 'off';
    axis square
    box off
    subplot(1,2,2)
    plot(Powerfreq2(1:65),(pwelch_fft(i,1:65)),'k','linewidth',3.5);
    hold on
    % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
    % title(MVL(1))
    hold on
    %  l=xline(peak_freq(i));
    %         l.LineStyle='--';
    %         l.Color=[170 170 170]/280;
    %         l.LineWidth=1.5;
    % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
    % title(MVL(1))
    ylabel('PSD (mV^2/Hz)');% ylabel('PSD ');
    xlabel('Frequency (Hz)');
    set(gca,'fontsize',30,'YColor','k','XColor','k');
    box off
    % h=xline(freq_wave,'--r','linewidth',2.5);
    % h=xline(freq_wave/2,'--b','linewidth',2.5);
    axis([0 inf 0 2])
    axis square
    title(['cell = ', num2str(sorting_ps2(j))])
end


%% Autocorrelations for V1 example session


clear X Y spikes spikes_d_s score mouse i FRp FR coeff cells_d T alfa cells_osc h lags peak_freq peak_psd period_ps ...
    Powerspec2 Powerfreq2 pwelch_fft_section pwelch_fft quality signal sorting_0 sorting_ascend sortind_descend sorting_freq2 ...
    sorting_ps2 val_sorting_freq2 val_sorting_ps2


%V1
mouse='92229';
load([rec_path,strcat('recording_dates_',mouse,'.mat')]);
day=7;
s=1;
munit=dates.ses{day}(s);
file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name,'-mat');
spikes=full(spikes_d_s);
[N,T]=size(spikes);
sf=7.73;
sf_d=7.73/4;

spikes_full=spikes;
clear spikes
spikes=spikes_full(:,1:T_PaS);


index_c=0;
index_t=0;
count=0;
clear period quality


downsampling_factor=4;
fs_120 = 7.73;
fs_30 = 30.95;
fs = fs_120/downsampling_factor;
maxlag=floor(560*fs_120);%240;%250 - 520;

[N,T]=size(spikes);

FR=spikes;
count=0;
for i=1:N
    [alfa(i,:),lags]=xcorr(FR(i,:),FR(i,:),maxlag,'coeff');
    [val(i,:),~]=zscore(alfa(i,:)); %Check whether I need to zscore
    signal=alfa(i,:);

    clear Powerspec2 Powerfreq2

    %Calculate spectogram using pwelch method
    [Powerspec2,Powerfreq2] = doPwelch(signal,fs_120,2*4096);
    pwelch_fft(i,:)=Powerspec2;
    pwelch_fft_section(i,:)=Powerspec2(5:100);


    %Check if the cell is oscillatory using the pwelch
    %method
    [peak_freq(i),peak_psd(i),quality(i)]=check_peak_quality_3c_f(Powerspec2,Powerfreq2,0.04);
    %         period_fft(i)=1/peak_fft(i);
    period_ps(i)=1/peak_freq(i);
    if quality(i)==1
        count=count+1;
        cells_osc(count)=i;
%         if count<50
%             figure
%             plot(Powerfreq2(1:100),Powerspec2(1:100))
%             title(i)
%         end
    else
%         figure
%         plot(Powerfreq2(1:100),Powerspec2(1:100))
%         title(i)
    end
    clear  Powerspec2 P2 P1  Y
end

[val_sorting_ps2,sorting_ps2]=sort(peak_psd,'descend');
[val_sorting_freq2,sorting_freq2]=sort(peak_freq,'descend');
[sorting_ascend,sortind_descend,sorting_0] = get_sorting(spikes);


figure
histogram(peak_freq,20)

length(find(peak_freq>0.008))

figure
[h,edges]=histcounts(peak_freq,[0:0.001:0.1]);

bar(edges(1:end-1),h);

% Figures autocorrelograms for all cells


% %PCA sorting
% figure
% imagesc(val(sortind_descend,:))
% caxis([0 1])
% colormap cividis;
% set(gca,'fontsize',18)
% xticks([1 maxlag maxlag*2])
% xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
% xlabel('Time (s)');
% ylabel('Neurons #');
% yticks([100 400])
% axis square
% colorbar
% caxis([0 1.5])

%Frequency sorting
figure
imagesc(val(sorting_freq2,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 200])
axis square
colorbar
title('Visual cortex - Freq sorting')


%Power sorting
figure
imagesc(val(sorting_ps2,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 200])
axis square
colorbar
title('Visual cortex - Power sorting')


%% Example cells from Vis

for j=[8,10]
    i=sorting_ps2(j);
    figure
    subplot(1,2,1)
    plot((alfa(i,:)),'k','linewidth',3.5)
    xticks([1 maxlag maxlag*2])
    xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
    xlabel('Time lag (s)');
    ylabel('Autocorrelation');
    axis([-inf inf 0 0.5])
    set(gca,'fontsize',30);
    % h = gca; h.YAxis.Visible = 'off';
    axis square
    box off
    subplot(1,2,2)
    plot(Powerfreq2(1:65),(pwelch_fft(i,1:65)),'k','linewidth',3.5);
    hold on
    % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
    % title(MVL(1))
    hold on
    %  l=xline(peak_freq(i));
    %         l.LineStyle='--';
    %         l.Color=[170 170 170]/280;
    %         l.LineWidth=1.5;
    % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
    % title(MVL(1))
    ylabel('PSD (mV^2/Hz)');% ylabel('PSD ');
    xlabel('Frequency (Hz)');
    set(gca,'fontsize',30,'YColor','k','XColor','k');
    box off
    % h=xline(freq_wave,'--r','linewidth',2.5);
    % h=xline(freq_wave/2,'--b','linewidth',2.5);
    axis([0 inf 0 2])
    axis square
    title(['cell = ', num2str(sorting_ps2(j))])
end
