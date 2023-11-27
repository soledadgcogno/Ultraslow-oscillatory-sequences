%% Autocorrelations and PSDs of NPX example cells - Both for theta and USO

clear all
npx_init;
mouseid='104638';
list = dir(strcat('W:\npxwaves\MEC\',mouseid));
% list = dir(strcat('W:\lechctime\for Soledad (temporary) 104640\'));
% addpath(genpath('C:\Users\xscogno\MATLAB\neuropixels'));

%Parameters
type_of_unit='good';
probe_num=1;


%Load data and get spike times
d=8;
disp(d)
dd=list(d).name;
pathf=strcat('W:\npxwaves\MEC\',mouseid,'\',dd);
[spk_times,gen,start,stop] = get_spike_times_Npx(pathf,probe_num,type_of_unit); %Check with Ane
spk_times_o=spk_times;
gen_o=gen;

%% USO
%Generate spike matrix
kernel_size=5;
thr=1.5;
% bin_size=0.01;
bin_size=0.12;
bin_size_uso=0.12;
[spk_g,bin_g]=gauss_spike_train_full_b(spk_times,start,stop,kernel_size,thr,bin_size);
[spk,bin] = spike_train_full_b(spk_times,start,stop,bin_size,thr);
N=size(spk,1);
T=size(spk,2);
fs_120=1/bin_size;

%Remove the first 300 s
time_stamps=(1:size(bin,2))*bin_size;
bin_o=bin;
bin(:,time_stamps<300)=[];

%Autocorrelation and PSD calculated on autocorrelation
maxlag=560;
clear autocorrelations
% [autocorrelations,lag_uso,PSD,Freq,peak_freq,peak_psd,quality]=calculate_PSD(bin,bin_size,maxlag,'autocorrelation');

maxlag=floor(maxlag*fs_120);
%PSD calculated on the autororrelations
for i=1:size(bin,1)    
    [alfa_uso(i,:),lags_uso]=xcorr(bin(i,:),bin(i,:),floor(maxlag),'coeff');
    [Powerspec3_uso(i,:),Powerfreq3_uso] = doPwelch(alfa_uso(i,:),fs_120,2*4096);
    [peak_freq(i),peak_psd(i),quality(i)]=check_peak_quality_3c_f(Powerspec3_uso(i,:)',Powerfreq3_uso,0.04);


    if i==1        
        ini3=find(Powerfreq3_uso>0.003,1);
        endi1=find(Powerfreq3_uso<0.1);
        endi3=endi1(end);
    end
    peak_uso_3(i)=max(Powerspec3_uso(i,ini3:endi3));
    clear  Powerspec2 P2 P1  Y
end

[psd_uso_local,psd_uso_local]=sort(peak_uso_3,'descend');
[psd_uso_gen,psd_uso_gen]=sort(peak_psd,'descend');


%% Final figures

%%%%%%%%%%%% Cell 
i=14;
%%%%%%%%%%%%%%%%%

%USO
figure
plot(lags_uso*bin_size_uso,alfa_uso(psd_uso_local(i),:),'k','linewidth',2)
% set(gca, 'YScale', 'log')
    axis([-200 200 0.00 0.2])
xticks([-200 0 200]);
xlabel('Time lag (s)');
ylabel('Autocorrelation');
% axis([-inf inf 0 0.2])
set(gca,'fontsize',20);
% h = gca; h.YAxis.Visible = 'off';
axis square
box off

close all

figure
plot(Powerfreq3_uso,Powerspec3_uso(psd_uso_local(i),:),'k','linewidth',2.5);
axis([0 0.06 0 0.1])
hold on
l=xline(peak_freq(psd_uso_local(i)));
l.LineStyle='--';
l.Color=[170 170 170]/255;
l.LineWidth=1.5;
% plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
% title(MVL(1))
ylabel('PSD (mV^2/Hz)');
%         ylabel('PSD ');
xlabel('Frequency (Hz)');
set(gca,'fontsize',20,'YColor','k','XColor','k');
box off
