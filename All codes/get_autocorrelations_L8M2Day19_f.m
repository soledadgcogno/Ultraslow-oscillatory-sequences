%% Autocorrelations
clear all
close all

rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

clusters=10;
mice_number=14;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L10M1';'L10M2';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];

index_c=0;
index_t=0;
count=0;
clear period quality

file_name_mvl='C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\MVL_all_sessions';
a=load(file_name_mvl,'-mat');
mvl=a.MVL_allsessions{8};

m=2;%9:mice_number

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

day=19;%:dates.daysnum
s=1;%:dates.sesnum(day)
disp(s)
munit=dates.ses{day}(s);
file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];


load([rec_data_path ['WS_Osc_14_','L08M2','.mat']]);
% dt
dt=WS_stat.dt(day,s);
if isinteger(dt)
else
    dt=floor(dt);
end

load(file_name,'-mat');
spikes=full(spikes_d_s);
downsampling_factor=4;
fs_120 = 7.73;
fs_30 = 30.95;
fs = fs_120/downsampling_factor;
maxlag=floor(560*fs_120);%240;%250 - 520;

[N,T]=size(spikes);
FRp = spikes_downsample(spikes,N,downsampling_factor);

FR=spikes;
count=0;
for i=1:N
    [alfa(i,:),~]=xcorr(FR(i,:),FR(i,:),maxlag,'coeff');
    [val(i,:),~]=zscore(alfa(i,:)); %Check whether I need to zscore
    signal=alfa(i,:);

    clear Powerspec2 Powerfreq2

    %Calculate spectogram using pwelch method
    [Powerspec2,Powerfreq2] = doPwelch(signal,fs_120,2*4096);
    pwelch_fft(i,:)=Powerspec2;

    pwelch_fft_section(i,:)=Powerspec2(5:100);


    %         P2_max=max(Powerspec2(3:end));
    %         ind2=find(Powerspec2==P2_max);
    %         peak_ps(i)=Powerfreq2(ind2);

    %Check if the cell is oscillatory using the pwelch
    %method
    [peak_freq(i),peak_psd(i),quality(i)]=check_peak_quality_3c_f(Powerspec2,Powerfreq2,0.04);
    %         period_fft(i)=1/peak_fft(i);
    period_ps(i)=1/peak_freq(i);
    if quality(i)==1
        count=count+1;
        cells_osc(count)=i;

%         if count>50 && count<150
%             figure
%             plot(Powerfreq2(1:100),Powerspec2(1:100))
%             title(i)
%         end

    else
        %                                 figure
        %                                 plot(Powerfreq2(1:100),Powerspec2(1:100))
        %                                 title(i)
    end
    clear  Powerspec2 P2 P1  Y
end


[val_sorting_ps2,sorting_ps2]=sort(peak_freq,'descend');



%% Figures autocorrelograms for all cells


%PSD for all cells
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



%PSD for oscillatory cells

aux=find(quality==0);
sorting_ps2_c=sorting_ps2;
sorting_ps2_c(aux)=[];

figure
imagesc(val(sorting_ps2_c,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 300])
axis square
hcb = colorbar;
hcb.Title.String = 'zscore(autocorrelation)';
title('PS sorting'); 

%PSD for non-oscillatory cells

aux2=find(quality==1);
sorting_ps2_c2=sorting_ps2;
sorting_ps2_c2(aux2)=[];

figure
imagesc(val(sorting_ps2_c2,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([10 30])
axis square
hcb = colorbar;

%% Mean autocorrelation

%Mean autocorrelation
mean_auto=mean(alfa(find(quality>0),:));
std_auto=std(alfa(find(quality>0),:));
std_auto=std(alfa(find(quality>0),:))/sqrt(count);

figure
x=1:size(alfa,2);
x2=[x,fliplr(x)];
inBetween = [mean_auto+std_auto, fliplr(mean_auto-std_auto)];
fill(x2, inBetween, [252 215 215]./255);
hold on
plot(mean_auto+std_auto,'w','linewidth',1.5);
hold on
plot(mean_auto-std_auto,'w','linewidth',1.5);
plot((mean(alfa(find(quality>0),:))),'k','linewidth',1);
xticks([1,length(mean_auto)/2,length(mean_auto)]);
xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
axis([-inf inf -inf 1])
xlabel('Time')
ylabel('Mean autocorrelation');
set(gca, 'YScale', 'log')
set(gca,'fontsize',16);


%Fraction of oscillatory cells
osc=count/N;
non_osc=(N-count)/N;

figure
bar([1,2],[osc,non_osc]*100);
axis([0.5 2.5 0 100]);
xticks([1,2]);
xticklabels({'Osc','Non-Osc'});
ylabel('Cell percentage %')
box off
set(gca,'fontsize',16)

%Period distribution for cells that oscillate
periods=1./peak_freq(cells_osc);
figure
histogram(periods,0:10:240);
ylabel('Counts');
xlabel('Cell period (s)');
set(gca,'fontsize',16);


% figure
% subplot(1,2,1)
% plot(peak_freq)
% title('Freq')
% subplot(1,2,2)
% plot(1./peak_freq)
% title('Period')
% 



%% Examples of individual cells 

[val_peak_psd,peak_psd_ind]=sort(peak_psd, 'descend');

N_non_osc=50;
N_osc=100;
for i=1:N_non_osc
    cell_non_osc(i)=peak_psd_ind(i); 
end

non_osc=N-count;
peak_psd_ind_c=peak_psd_ind;
peak_psd_ind_c(1:non_osc)=[];
val_peak_psd_c=val_peak_psd;
val_peak_psd_c(1:non_osc)=[];

for i=1:N_osc
    cell_osc(i)=peak_psd_ind_c(i); 
end

countnonosc=0;
for i=(find(quality==0))
    countnonosc=countnonosc+1;
    if countnonosc<10
        figure
        plot((alfa(i,:)),'k','linewidth',2)
        % set(gca, 'YScale', 'log')
        xticks([1 maxlag maxlag*2])
        xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
        xlabel('Time (s)');
        ylabel('Autocorrelation');
        axis([-inf inf 0 0.2])
        set(gca,'fontsize',20);
        % h = gca; h.YAxis.Visible = 'off';
        axis square
        box off
        figure
        plot(Powerfreq2(1:65),(pwelch_fft(i,1:65)),'k','linewidth',2.5);
        hold on
        % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
        % title(MVL(1))
        % ylabel('PSD (mV^2/Hz)');
        ylabel('PSD ');
        xlabel('Frequency (Hz)');
        set(gca,'fontsize',20);
        box off
        % h=xline(freq_wave,'--r','linewidth',2.5);
        % h=xline(freq_wave/2,'--b','linewidth',2.5);
        axis square
    end
end


for i=71:80%31:50%40:N_osc
figure
% subplot(1,2,1)
plot((alfa(cell_osc(i),:)),'k','linewidth',2)
xticks([1 maxlag maxlag*2])
xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Autocorrelation');
axis([-inf inf 0 0.3])
set(gca,'fontsize',20);
% h = gca; h.YAxis.Visible = 'off';
axis square
box off
figure
plot(Powerfreq2(1:65),(pwelch_fft(cell_osc(i),1:65)),'k','linewidth',2.5);
hold on
% plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
% title(MVL(1))
ylabel('PSD (mV^2/Hz)');
% ylabel('PSD ');
xlabel('Frequency (Hz)');
set(gca,'fontsize',20);
box off
% h=xline(freq_wave,'--r','linewidth',2.5);
% h=xline(freq_wave/2,'--b','linewidth',2.5);
axis square
% title(round(mvl(cell_osc(i)),2))
end

%cell=2
figure
plot((alfa(2,:)),'k','linewidth',2)
xticks([1 maxlag maxlag*2])
xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Autocorrelation');
axis([-inf inf 0 0.3])
set(gca,'fontsize',20);
axis square
box off
figure
plot(Powerfreq2(1:65),(pwelch_fft(2,1:65)),'k','linewidth',2.5);
hold on
ylabel('PSD (mV^2/Hz)');
xlabel('Frequency (Hz)');
set(gca,'fontsize',20);
box off
axis square

%% Anatomical distribution of frequencies

mouse='L8M2';
day=19;
munit=0;

file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_anat,'-mat');
pixel_size_new=1.18185;


b=discretize(peak_freq,[0.0017,3.6*0.0018,4.1*0.0018,5.1*0.0018,7*0.0018,10]);
cc=copper(5);
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_anat,'-mat');
figure
hold on
for i=1:N
    x=Anat.pos{1,i}(:,1);
    y=Anat.pos{1,i}(:,2);
    ov=Anat.overlap{1,i};
    
    if peak_freq(i)==Inf
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
    else
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(b(i),:),'filled');
    end
    clear x y ov
end
box on
xlabel('X [\mum]');
ylabel('Y [\mum]');
axis([0 580 20 600])
yticks([120 420]);
yticklabels([100 400]);
xticks([100 400]);
set(gca,'fontsize',16)
axis square
colormap copper(5)
% co=colorbar('XTick',0:1,'XTickLabel',{'555','23'});
co.Label.String = 'Freq (Hz)';



%% Anatomical distribution of periods

mouse='L8M2';
day=19;
munit=0;

file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_anat,'-mat');
pixel_size_new=1.18185;


b=discretize(period_ps,[10,30,50,70,90,120,160,500]);
cc=hsv(9);
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_anat,'-mat');
figure
hold on
for i=1:N
    x=Anat.pos{1,i}(:,1);
    y=Anat.pos{1,i}(:,2);
    ov=Anat.overlap{1,i};
    
    if period_ps(i)==0
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
    else
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(b(i),:),'filled');
    end
    clear x y ov
end
box on
xlabel('X [\mum]');
ylabel('Y [\mum]');
axis([0 580 20 600])
yticks([120 420]);
yticklabels([100 400]);
xticks([100 400]);
set(gca,'fontsize',16)
axis square
colormap hsv(9)
% co=colorbar('XTick',0:1,'XTickLabel',{'555','23'});
co.Label.String = 'Period (s)';



% figure
% plot(period)
% val_p=val(cells,:);
% [~,sorted]=sort(period_ps,'ascend');
% figure
% imagesc(val_p(sorted,:))
% caxis([0 1])
% colormap cividis;
% xticklabels({-240 0 240});
% xlabel('Time (s)');
% ylabel('Neurons #');
% yticks([100 300])
% axis square
% colorbar
% 
% figure
% plot(period)
% val_p=val(cells,:);
% [~,sorted]=sort(period_fft,'ascend');
% figure
% imagesc(val_p(sorted,:))
% caxis([0 1])
% colormap cividis;
% xticklabels({-240 0 240});
% xlabel('Time (s)');
% ylabel('Neurons #');
% yticks([100 300])
% axis square
% colorbar

osc=cells_osc;
not_osc=setdiff(1:N,cells_osc);


% figure
% labels = {'Oscillatory','Not oscillatory'};
% pie([length(osc)/N,length(not_osc)/N],labels)
% colormap jet(4)
% set(gca,'fontsize',16)
% 
% figure
% plot(val(osc(1),:),'k')
% xticks([1 maxlag maxlag*2])
% xticklabels({-280 0 280});
% box off
% figure
% plot(val(osc(2),:),'k')
% xticks([1 maxlag maxlag*2])
% xticklabels({-280 0 280});
% box off
% figure
% plot(val(not_osc(2),:),'k')
% xticks([1 maxlag maxlag*2])
% xticklabels({-280 0 280});
% box off

%% Locking

disc_phase=10;

% Compute radius and MI
[coeff,score,latent_pca] = pca(spikes');
phase=atan2(score(:,2),score(:,1));
radius=sqrt(coeff(:,2).*coeff(:,2)+coeff(:,1).*coeff(:,1));

% Compute tuning curves and locking using deconvolved spikes
phase_d=downsample(phase,4);
phase_di=discretize(phase_d,-pi:2*pi/disc_phase:pi);
% phase_d=discretize(phase,-pi:2*pi/disc_phase:pi);

new_bin = 4;
for i=1:length(phase_d)
    sp_do(:,i)=sum(spikes(:,(i-1)*new_bin+1:i*new_bin),2);
end

for i=1:N
    p=phase(find(spikes(i,:)));
    var_p(i)=circ_var(p);
    mean_p(i)=circ_mean(p);
    std_p(i)=circ_std(p);
    MVL(i) = circ_r(p);
    FR_tot(i)=length(find(spikes(i,:)));
    H=histogram(p,-pi:2*pi/40:pi,'Normalization','Probability');
    prob_phase_firing(i,:)=H.Values;
    clear p H
end

i=1;
p=phase(find(spikes(i,:)));
H=histogram(p,-pi:2*pi/40:pi,'Normalization','Probability');
phase_TC=H.BinEdges;

for i=1:N
    disp(i)
    for sh=1:100
        p_sh=phase(floor(T*rand(1,length(find(spikes(i,:)))))+1);
        var_p__sh(i,sh)=circ_var(p_sh);
        std_p__sh(i,sh)=circ_std(p_sh);
        mean_p_sh(i,sh)=circ_mean(p_sh);
        MVL_sh(i,sh) = circ_r(p_sh);
%         if sh==1
%             H=histogram(p_sh,-pi:2*pi/40:pi,'Normalization','Probability');
%             prob_phase_firing_sh(i,:)=H.Values;
%             clear H
%         end
%         H=histogram(p_sh,-pi:2*pi/40:pi,'Normalization','Probability');
%         prob_phase_firing_sh_tot(i,:,sh)=H.Values;
%         clear H
        
        clear p_sh
    end
end

locking=1-var_p;
locking_sh=1-var_p__sh;
locking_sh_mean=mean(MVL_sh,2);
locking_sh_99=prctile(MVL_sh,99,2);
locking_sh_1=prctile(MVL_sh,1,2);

cells_mi_calc=find(locking>locking_sh_99');
fraction_mi_calc=length(cells_mi_calc)./N;
fraction_mi_calc2=length(find(locking<=locking_sh_99'))./N;
tot_cells=1:N;
not_locked=find(locking<=locking_sh_99');
locked=setdiff(tot_cells,not_locked);

%% Look at overlaps


cells_locked_osc=(intersect(locked,osc));
cells_locked_notosc=(intersect(locked,not_osc));
cells_notlocked_osc=(intersect(not_locked,osc));
cells_notlocked_notosc=(intersect(not_locked,not_osc));

locked_osc=length(intersect(locked,osc))/N;
locked_notosc=length(intersect(locked,not_osc))/N;
notlocked_osc=length(intersect(not_locked,osc))/N;
notlocked_notosc=length(intersect(not_locked,not_osc))/N;


figure
labels = {'Locked - Osc','Locked - Not Osc','Not Locked - Osc','Not locked - Not Osc'};
pie([locked_osc,locked_notosc,notlocked_osc,notlocked_notosc],labels)
colormap jet(4)
set(gca,'fontsize',16)


figure
subplot(2,1,1)
plot(val(cells_locked_osc(1),:),'k')
xticks([1 maxlag maxlag*2])
xticklabels({-280 0 280});
box off
subplot(2,1,2)
[Powerspec2,Powerfreq2] = doPwelch(val(cells_locked_osc(1),:),2,512*2);
plot(Powerfreq2(1:50),Powerspec2(1:50),'k');

figure
subplot(2,1,1)
plot(val(cells_locked_notosc(1),:),'k')
xticks([1 maxlag maxlag*2])
xticklabels({-280 0 280});
box off
subplot(2,1,2)
[Powerspec2,Powerfreq2] = doPwelch(val(cells_locked_notosc(1),:),2,512*2);
plot(Powerfreq2(1:50),Powerspec2(1:50),'k');

figure
subplot(2,1,1)
plot(val(cells_notlocked_osc(1),:),'k')
xticks([1 maxlag maxlag*2])
xticklabels({-280 0 280});
box off
subplot(2,1,2)
[Powerspec2,Powerfreq2] = doPwelch(val(cells_notlocked_osc(1),:),2,512*2);
plot(Powerfreq2(1:50),Powerspec2(1:50),'k');


figure
subplot(2,1,1)
plot(val(cells_notlocked_notosc(1),:),'k')
xticks([1 maxlag maxlag*2])
xticklabels({-280 0 280});
box off
subplot(2,1,2)
[Powerspec2,Powerfreq2] = doPwelch(val(cells_notlocked_notosc(1),:),2,512*2);
plot(Powerfreq2(1:50),Powerspec2(1:50),'k');
% 
% 
% figure
% plot(val(cells_locked_osc(1),:),'k')
% xticks([1 maxlag maxlag*2])
% xticklabels({-280 0 280});
% box off
% 
% figure
% plot(val(cells_locked_osc(1),:),'k')
% xticks([1 maxlag maxlag*2])
% xticklabels({-280 0 280});
% box off






figure
plot(MVL(not_osc),'*')

% aux=find(MVL(not_osc)>0.7);
% 
% figure
% plot(val(27,:))
% 
% [Powerspec2,Powerfreq2] = doPwelch(val(79,:),2,512);
% P2_max=max(Powerspec2(4:end));
% ind2=find(Powerspec2==P2_max);
% peak_ps(i)=Powerfreq2(ind2);
% quality(i)=check_peak_quality_3b(Powerspec2);
% period_fft(i)=1/peak_fft(i);
% period_ps(i)=1/peak_ps(i);

%% Figures locking

figure
labels = {'Locked','Not locked'};
pie([fraction_mi_calc,fraction_mi_calc2],labels)
colormap jet(4)
set(gca,'fontsize',16)

figure
plot(1:N,MVL,'*');
hold on
plot(not_osc,MVL(not_osc),'*');

figure
plot(FR_tot,MVL,'*')
ylabel('MVL');
xlabel('Total # frames of 0.5 seconds with spikes');
box off
set(gca,'fontsize',16)

figure
scatter(1:N,1./peak_fft((sorting_fft)),45,'o','filled');
alpha 0.4
xlabel('Neuron #');
ylabel('Period (s)');
set(gca,'fontsize',16);
% axis([-2 400 0 200])
not_locked_fft=1./(peak_fft(not_locked));

figure
colormap jet
scatter(1:N,1./peak_fft,45,'o','filled');
alpha 0.4
xlabel('Neuron #');
ylabel('Period (s)');
set(gca,'fontsize',16);
% axis([-2 400 0 200])
not_locked_fft=1./(peak_fft(not_locked));
hold on
scatter(not_locked,not_locked_fft,'filled')
not_osc_fft=1./(peak_fft(not_osc));
hold on
scatter(not_osc,not_osc_fft,'filled')
legend({'All','Not locked','Not Osc'})


% for i=1:length(locked)
%     aux(i)=find(sorting_pca==locked(i));
% end
% figure
% plot(PI(aux),peak_fft(sorting_pca(aux)),'*')

% figure
% plot(MVL,peak_fft,'*')
% hold on
% plot(MVL(not_locked),peak_fft(not_locked),'*');
% xlabel('MVL');
% ylabel('Peak fft');
% set(gca,'fontsize',16)

figure
plot(MVL,1./peak_fft,'*')
hold on
plot(MVL(not_locked),1./peak_fft(not_locked),'*');
xlabel('MVL');
ylabel('Period (s)');
set(gca,'fontsize',16)

% figure
% plot(MVL,period_fft,'*')
% hold on
% plot(MVL(not_osc),period_fft(not_osc),'*');
% xlabel('MVL');
% ylabel('Period');
% set(gca,'fontsize',16)

figure
plot(FR_tot,1./peak_fft,'*')
hold on
plot(FR_tot(not_locked),1./peak_fft(not_locked),'*');
xlabel('Total # frames of 0.5 seconds with spikes');
ylabel('Period');
set(gca,'fontsize',16)


%% Rasterplot

per=find(period_fft(osc)<30);
loc=find(MVL(osc)>0.7);
int=intersect(per,loc);
int_pca=find(sorting_pca==osc(int(2)));

per2=find(period_fft(osc)>200);
loc2=find(MVL(osc)>0.7);
int2=intersect(per2,loc2);
int_pca2=find(sorting_pca==osc(int2(2)));


fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes,1)
    scatter((1:size(spikes,2))./8,i*spikes(sorting_pca(i),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf]);
%title([mouse,' Day',num2str(day)]);
xlabel('Time (s)');
ylabel('Neurons #');
set(gca,'fontsize',18);
yticks([100 400])
hold on
scatter((1:size(spikes,2))./8,int_pca*spikes(osc(int(1)),:),100,'b','filled')
scatter((1:size(spikes,2))./8,int_pca2*spikes(osc(int2(1)),:),100,'r','filled')

    

figure
subplot(2,1,1)
plot(val(osc(int(2)),:),'b')
xticks([1 maxlag maxlag*2])
xticklabels({-280 0 280});
box off
subplot(2,1,2)
[Powerspec2,Powerfreq2] = doPwelch(val(osc(int(2)),:),2,512*2);
plot(Powerfreq2(1:100),Powerspec2(1:100),'b');

figure
subplot(2,1,1)
plot(val(osc(int2(2)),:),'r')
xticks([1 maxlag maxlag*2])
xticklabels({-280 0 280});
box off
subplot(2,1,2)
[Powerspec2,Powerfreq2] = doPwelch(val(osc(int2(2)),:),2,512*2);
plot(Powerfreq2(1:100),Powerspec2(1:100),'r');




figure
plot(val(osc(int(1)),:),'b','linewidth',1.5)
box off
   
figure
plot(val(osc(int2(1)),:),'r','linewidth',1.5)
box off
