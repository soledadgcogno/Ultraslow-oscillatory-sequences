% close all
clear all
npx_init;
clear all
npx_init;

%General parameters
color_shanks=0;
starting_point_plot=300; %in seconds

%%%% KIMI
% % mouseid='104638';
% % list = dir(strcat('W:\npxwaves\MEC\',mouseid));
% % % Parameters
% % type_of_unit='good';
% % probe_num=1;
% % dt=1;
% % thr=2.5;
% % chunk_ini=1200;
% % chunk_end=1700;
% % % session
% % d=8;
% % dd=list(d).name;
% % pathf=strcat('W:\npxwaves\MEC\',mouseid,'\',dd);


%%%%% KNUD
mouseid='102335';
list = dir(strcat('W:\npxwaves\MEC_HC\',mouseid));
%parameters
type_of_unit='all';
probe_num=1;
dt=1;
thr=2.5;
% chunk_ini=940; %Try with 1000
% chunk_end=1318;
chunk_ini=1100; %Try with 1000
chunk_end=1400;
%session
d=18;
dd=list(d).name;
pathf=strcat('W:\npxwaves\MEC_HC\',mouseid,'\',dd);


%Generates all matrices of activity with Gaussian kernel
[spk_times,~,start,stop] = get_spike_times_Npx(pathf,probe_num,type_of_unit);
% for i=1:size(spk_times,1)
%     FR(i)=length(spk_times{i,1})/(stop-start);
% end

% shanks_p=0;
% [spk_times,gen] = remove_shanks(shanks_p,gen,spk_times);
kernel_size=5;
thr=1;
downsample_new_bin=0.12;%0.12; %number of samples per second
sf=1/downsample_new_bin; %sampling frequency
[~,bin_s]=gauss_spike_train_full_b(spk_times,start,stop,kernel_size,thr,downsample_new_bin);


% Sequence identification
starting_point_plot_b=floor(starting_point_plot*(1/downsample_new_bin));
duration_waves=nan(1,100);
% spikes=spk_trace_aux;%spk_trace_aux_full(:,starting_point_plot:end);
num_clus_discr=10;
dt_s=5;%in seconds 5 seconds for Kimi
make_fig=1;
[table_u,spk_chunk,spikes,sorting,cells2keep]=identify_waves_latestversion_6_f_forNpx(num_clus_discr,downsample_new_bin,dt_s,make_fig,bin_s,chunk_ini,chunk_end);


%Concatenate sequences
spikes_w=[];
for i=1:size(table_u,1)
    spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
end

total_events=sum(spikes_w,2);
spikes_w_f=spikes_w;
spikes_w_f(total_events<50,:)=[];

figure
spy(spikes_w(sorting,:));
pbaspect([8 1 1])

number_waves=size(table_u,1);
duration_waves=(table_u(:,2)-table_u(:,1))/sf;
spikes_wo5min=spikes(:,starting_point_plot_b:end);

%Oscillation score
[fraction,dt_osci] = get_oscillation_score_npx_b(spk_chunk,spikes_w,downsample_new_bin);


%Fraction of oscillatory cells
epoch_length=20;
ptl=95; %percentile to determine significance
[fraction_osc,osc,tot] = get_oscillatory_cells_npx(spikes_w,downsample_new_bin,epoch_length,ptl);


clear table_u spikes dt spikes_d_s row_w spikes_w

%% Make figures Rasterplots
make_rasterplot(spikes(:,starting_point_plot_b:end),flip(sorting),downsample_new_bin,0); %alpha 0.1
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Neuropixels\Rasterplots_Npx\Knud_Bin_Gaussian.png'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Neuropixels\Rasterplots_Npx\Knud_Bin_Gaussian.svg'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Neuropixels\Rasterplots_Npx\Knud_Bin_Gaussian.fig'));
close all

%Raster plot on matrix obtained *without* Gaussian convolution -
bin_size=1;
starting_point_plot_c=floor(starting_point_plot*(1/bin_size));
thr=1.5;
[spk_from_full,bin_from_full] = spike_train_full_b(spk_times,start,stop,bin_size,thr);
mat=bin_from_full(cells2keep,:);
make_rasterplot(mat(:,starting_point_plot_c:end),flip(sorting),bin_size,0); %alpha 0.5
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Neuropixels\Rasterplots_Npx\Knud_Bin.png'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Neuropixels\Rasterplots_Npx\Knud_Bin.svg'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Neuropixels\Rasterplots_Npx\Knud_Bin.fig'));
close all

%% Make figures Autocorrelations (only for Kimi)

thr=1.5;
bin_size=0.12;
[spk,bin_uto] = spike_train_full_b(spk_times,start,stop,bin_size,thr);
N=size(spk,1);
T=size(spk,2);
fs_120=1/bin_size;

%Remove the first starting_point_plot seconds
time_stamps=(1:size(bin_uto,2))*bin_size;
bin_o=bin_uto;
bin_uto(:,time_stamps<starting_point_plot)=[];

%Autocorrelation and PSD calculated on autocorrelation
maxlag=250;
clear autocorrelations
[autocorrelations,lag_uso,PSD,Freq,peak_freq,peak_psd,quality]=calculate_PSD(bin_uto,bin_size,maxlag,'autocorrelation');

%PSD calculated on spike train
for i=1:size(bin,1)    
    [Powerspec3_uso,Powerfreq3_uso] = doPwelch(bin(i,:),fs_120,2*4096);
    if i==1        
        ini3=find(Powerfreq3_uso>0.003,1);
        endi1=find(Powerfreq3_uso<0.1);
        endi3=endi1(end);
    end
    peak_uso_3(i)=max(Powerspec3_uso(ini3:endi3));
    clear  Powerspec2 P2 P1  Y
end


%i=192 (panel a)
figure
plot(lag_uso*bin_size,autocorrelations(i,:),'k','linewidth',2)
axis([-200 200 0.02 0.2])
xticks([-200 0 200]);
xlabel('Time lag (s)');
ylabel('Autocorrelation');
set(gca,'fontsize',20);
axis square
box off
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Neuropixels\USOs in single cells\USO_Autocorr_cell_',num2str(psd_uso_2_auto(i)),'.fig'));        
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Neuropixels\USOs in single cells\USO_Autocorr_cell_',num2str(psd_uso_2_auto(i)),'.svg'));        
% close all

figure
plot(Freq,PSD(psd_uso_2_auto(i),:),'k','linewidth',2.5);
axis([0 0.1 0 0.1])
hold on
l=xline(peak_freq(psd_uso_2_auto(i)));
l.LineStyle='--';
l.Color=[170 170 170]/255;
l.LineWidth=1.5;
ylabel('PSD (mV^2/Hz)');
xlabel('Frequency (Hz)');
set(gca,'fontsize',20,'YColor','k','XColor','k');
box off
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Neuropixels\USOs in single cells\USO_PSD_cell_',num2str(psd_uso_2_auto(i)),'.fig'));        
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Neuropixels\USOs in single cells\USO_PSD_cell_',num2str(psd_uso_2_auto(i)),'.svg'));        
