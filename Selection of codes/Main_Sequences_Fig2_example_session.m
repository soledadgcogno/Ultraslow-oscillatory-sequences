clear all
load('calcium_activity_matrix_60584_session17.mat')

%Preprocessing
dt=117;
spikes=full(spikes_d_s);
downsampling_factor=4;
fs_120 = 7.73;
fs_30 = 30.95;
fs = fs_120/downsampling_factor;
maxlag=floor(560*fs_120);
[N,T]=size(spikes);
FRp = spikes_downsample(spikes,N,downsampling_factor);

%% Oscillations in single cells

FR=spikes;
count=0;
for i=1:N
    [alfa(i,:),lags]=xcorr(FR(i,:),FR(i,:),maxlag,'coeff');
    [val(i,:),~]=zscore(alfa(i,:));
    signal=alfa(i,:);
    clear Powerspec2 Powerfreq2
    %Calculate spectogram using pwelch method
    [Powerspec2,Powerfreq2] = doPwelch(signal,fs_120,2*4096);
    pwelch_fft(i,:)=Powerspec2;
    pwelch_fft_section(i,:)=Powerspec2(5:100);

    %Check if the cell is oscillatory using the pwelch
    %method
    [peak_freq(i),peak_psd(i),quality(i)]=check_peak_quality_3c_f(Powerspec2,Powerfreq2,0.04);
    period_ps(i)=1/peak_freq(i);
    if quality(i)==1
        count=count+1;
        cells_osc(count)=i;
    else
    end
    clear  Powerspec2 P2 P1  Y
end

[val_sorting_ps2,sorting_ps2]=sort(peak_psd,'descend');
[~,sortind_descend,~] = get_sorting(spikes);

% Figures autocorrelograms for all cells
%Power sorting
figure
imagesc(val(sorting_ps2,:))
caxis([0 0.5])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 400])
axis square
colorbar

%Example cells
cell1=288;
cell2=323;

figure
subplot(1,2,1)
plot((alfa(cell1,:)),'k','linewidth',2)
xticks([1 maxlag maxlag*2])
xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Autocorrelation');
axis([-inf inf 0 0.3])
set(gca,'fontsize',20);
axis square
box off
subplot(1,2,2)
plot(Powerfreq2(1:65),(pwelch_fft(cell1,1:65)),'k','linewidth',2.5);
hold on
l=xline(peak_freq(cell1));
l.LineStyle='--';
l.Color=[170 170 170]/280;
l.LineWidth=1.5;
ylabel('PSD (mV^2/Hz)');
xlabel('Frequency (Hz)');
set(gca,'fontsize',20,'YColor','k','XColor','k');
box off
axis square

figure
subplot(1,2,1)
plot((alfa(cell2,:)),'k','linewidth',2)
xticks([1 maxlag maxlag*2])
xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Autocorrelation');
axis([-inf inf 0 0.3])
set(gca,'fontsize',20);
axis square
box off
subplot(1,2,2)
plot(Powerfreq2(1:65),(pwelch_fft(cell2,1:65)),'k','linewidth',2.5);
hold on
l=xline(peak_freq(cell2));
l.LineStyle='--';
l.Color=[170 170 170]/280;
l.LineWidth=1.5;
ylabel('PSD (mV^2/Hz)');
xlabel('Frequency (Hz)');
set(gca,'fontsize',20,'YColor','k','XColor','k');
box off
axis square



%% Sortings and rasterplots


spikes_d=full(spikes_d_s);
[N,T]=size(spikes_d);

%Preprocessing of data
downsampling_factor=4;
X = spikes_downsample(spikes_d,N,downsampling_factor);
FRp = spikes_downsample(spikes_d,N,downsampling_factor);
for i=1:N
    FR(i,:)=full(fire_rate(FRp(i,:),29,'g')); %Smoothing in about 10 seconds
end

% PCA
[coeff,score]=pca(spikes_d');
[sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
[~,sorting_w,~]=get_sorting_smoothed(spikes_d,117);
sorting_pca=sorting_descend;
ini=sorting_pca(1);
ini_w=sorting_w(1);

%Xcorr
maxlag=10;
downsampling_factor=4;
dt=117;
dt_sf=117/downsampling_factor;
sorting_corr=sorting_xcorr(spikes_d,maxlag,downsampling_factor,dt_sf);
aux=find(sorting_corr==ini);
sorting_corr_ini=circshift(sorting_corr,-(aux-1));
aux_w=find(sorting_corr==ini_w);
sorting_corr_ini_w=circshift(sorting_corr,-(aux_w-1));

addpath("C:\Users\xscogno\MATLAB\drtoolbox.tar");

%tSNE
P=30;
[Y] = tsne(FR,'NumDimensions',2,'NumPCAComponents',50,'Perplexity',P);
Y_tsne=Y;
angletsne=atan2(Y(:,2),Y(:,1));
[alfa,sorting_tsne]=sort(angletsne,'descend');
figure
spy(X(sorting_tsne,:),'k')
pbaspect([23,2,1])
aux=find(sorting_tsne==ini);
sorting_tsne_ini=circshift(sorting_tsne,-(aux-1));
aux_w=find(sorting_tsne==ini_w);
sorting_tsne_ini_w=circshift(sorting_tsne,-(aux_w-1));

%LEM
K=15;
[Y,~] = laplacian_eigen(FR,2,15,2);
anglelem=atan2(Y(:,2),Y(:,1));
[alfa,sorting_lem]=sort(anglelem,'descend');
aux=find(sorting_lem==ini);
sorting_lem_ini=circshift(sorting_lem,-(aux-1));
aux_w=find(sorting_tsne==ini_w);
sorting_lem_ini_w=circshift(sorting_lem,-(aux_w-1));

%ISOMAP
[Y_iso] = isomap(FR,2,15);
angleiso=atan2(Y_iso(:,2),Y_iso(:,1));
[alfa,sorting_iso]=sort(angleiso,'descend');
aux=find(sorting_iso==ini);
sorting_iso_ini=circshift(sorting_iso,-(aux-1));
aux_w=find(sorting_iso==ini_w);
sorting_iso_ini_w=circshift(sorting_iso,-(aux_w-1));

%UMAP
%                     downsampling_factor=8;
%                     X = spikes_downsample(spikes_d,N,downsampling_factor);

%30 and 0.1
for neigh=30%:20:70
    for  min_dis=0.3%:0.2:0.5
        [reduction]=run_umap(FR,'n_neighbors',neigh,'min_dist',min_dis,'metric','correlation');
        figure
        scatter(reduction(:,1),reduction(:,2),[],[ 0.2588    0.2588    0.2588],'filled')
        alpha 0.7
        set(gca,'fontsize',16)
        axis([-10 10 -10 10])
        ylabel('UMAP - Dim 2')
        xlabel('UMAP - Dim 1')
        xticks([-10 -5 0 5 10])
        yticks([-10 -5 0 5 10])
        set(gca,'fontsize',18)
        title(['Neigh:', num2str(neigh),' - Min Dist:',num2str(min_dis) ])
        com(:,1)=mean(reduction(:,1));
        com(:,2)=mean(reduction(:,2));
        %Now I move the points
        reduction_m(:,1)=reduction(:,1)-com(:,1);
        reduction_m(:,2)=reduction(:,2)-com(:,2);
        angleumap=atan2(reduction_m(:,2),reduction_m(:,1));
        [alfa,sorting_umap]=sort(angleumap,'descend');
        figure
        spy((X(sorting_umap,:)),'k');
        pbaspect([23,2,1])
        title(['Neigh:', num2str(neigh),' - Min Dist:',num2str(min_dis) ])
        %                             clear reduction
    end
end

com(:,1)=mean(reduction(:,1));
com(:,2)=mean(reduction(:,2));
%Now I move the points
reduction_m(:,1)=reduction(:,1)-com(:,1);
reduction_m(:,2)=reduction(:,2)-com(:,2);
angleumap=atan2(reduction_m(:,2),reduction_m(:,1));
[alfa,sorting_umap]=sort(angleumap,'descend');
aux=find(sorting_umap==ini);
sorting_umap_ini=circshift(sorting_umap,-(aux-1));
aux_w=find(sorting_umap==ini_w);
sorting_umap_ini_w=circshift(sorting_umap,-(aux_w-1));


rmpath("C:\Users\xscogno\MATLAB\drtoolbox.tar");



% Raster plots
close all

% PCA
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.3]);
sorting_descend=flip(sorting_descend);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_descend(i)),:),5,'k','filled')
    alpha 0.2
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 400])
ini=sorting_ascend(1);
set(gcf, 'Renderer', 'opengl');
title('PCA')

%Correlations
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.3]);
hold on
sorting_corr_ini=flip(sorting_corr_ini);
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_corr_ini(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 400])
set(gcf, 'Renderer', 'opengl');
title('Correlations')

%tsne
% sorting_tsne_ini=flip(sorting_tsne_ini);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.3]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_tsne_ini(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 400])
set(gcf, 'Renderer', 'opengl');
title('tsne')

%lem
sorting_lem_ini=flip(sorting_lem_ini);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.3]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_lem_ini(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 400])
set(gcf, 'Renderer', 'opengl');
title('LEM')

%isomap
sorting_iso_ini=flip(sorting_iso_ini);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.3]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_iso_ini(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 400])
set(gcf, 'Renderer', 'opengl');
title('Isomap')

%umap
sorting_umap_ini=flip(sorting_umap_ini);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.3]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_umap_ini(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 400])
set(gcf, 'Renderer', 'opengl');
title('UMAP')

