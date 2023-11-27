clear all
close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';


mice_number=6;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];
% mice=['L8M1';'L8M2';'L8M3';'L8M4';'L9M1';'L9M4'];


for m=6
        
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
    
    for day=17
        
        for s=1%:dates.sesnum(day)
            disp(s)
            munit=dates.ses{day}(s);
            file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
            
%             load([dpath ['WS_Osc_14_L08M4','.mat']]);
%             dt=WS_stat.dt(day,s);
            
            if (exist(file_name) == 2)
                disp(day)
                load(file_name,'-mat');
                spikes_d=full(spikes_d_s);
                [N,T]=size(spikes_d);
                
                if N>150
                    
                    %Preprocessing of data
                    downsampling_factor=4;
                    X = spikes_downsample(spikes_d,N,downsampling_factor);
                    FRp = spikes_downsample(spikes_d,N,downsampling_factor);
                    for i=1:N
                        FR(i,:)=full(fire_rate(FRp(i,:),29,'g')); %Smoothing in about 10 seconds
                    end
                    
                    
                    % Angle PCA
                    [coeff,score]=pca(spikes_d');
                    [sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
                    [~,sorting_w,~]=get_sorting_smoothed(spikes_d,117);
                    sorting_pca=sorting_descend;
                    figure
                    spy(spikes_d(sorting_pca,:),'k')
                    pbaspect([23,2,1])
                    ini=sorting_pca(1);
                    ini_w=sorting_w(1);
                    
                    %Xcorr
                    maxlag=10;
                    downsampling_factor=4;
                    dt=117;
                    dt_sf=117/downsampling_factor;
                    sorting_corr=sorting_xcorr(spikes_d,maxlag,downsampling_factor,dt_sf); %check the smoothing kernel
                    aux=find(sorting_corr==ini);
                    sorting_corr_ini=circshift(sorting_corr,-(aux-1));
                    aux_w=find(sorting_corr==ini_w);
                    sorting_corr_ini_w=circshift(sorting_corr,-(aux_w-1));
                    figure
                    spy(spikes_d(sorting_corr_ini,:),'k')
                    pbaspect([23,2,1])
                    figure
                    spy(spikes_d(sorting_corr_ini,:),'k')
                    pbaspect([23,2,1])
                    
                    


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
                    figure
                    spy(X(flip(sorting_tsne_ini),:),'k')
                    pbaspect([23,2,1])
                    
                    
                    %                     spk=spk_down(beta,:);
                    %
                    %                     figure
                    %                     set(gcf, 'PaperSize', [4 2]);
                    %                     spy(X(sorting_tsne,:),'k')
                    %                     pbaspect([25 3 1])
                    %                     title(['P ',num2str(P),' Mouse',num2str(mice(m,1:4)),' Day',num2str(day),' FoV',num2str(f)]);
                    %                     xlabel('Time [s]');
                    %                     ylabel('Cells');
                    %                     set(gca,'fontsize',18);
                    %                     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
                    
                    %LEM
                    %                     downsampling_factor=4;
                    %                     X = spikes_downsample(spikes_d,N,downsampling_factor);
                    %
                    %                     for i=1:N
                    %                         FR(i,:)=zscore(full(fire_rate(X(i,:),29,'g'))); %Smoothing in about 10 seconds
                    %                     end
                    %
                    %                     K=15;
                    %                     [Y] = laplacian_eigen(X,5,K,1);
                    %                     Y_lem=Y;
                    %                     anglelem=atan2(Y(:,2),Y(:,1));
                    %                     [alfa,sorting_lem]=sort(anglelem,'descend');
                    %
                    % %                     figure
                    % %                     scatter(Y_lem(:,1),Y_lem(:,2),[],[ 0.2588    0.2588    0.2588],'filled')
                    % %                     alpha 0.7
                    % %                     set(gca,'fontsize',16)
                    % %                     axis([-0.03 0.03 -0.03 0.03])
                    
                    K=15;
                    [Y,~] = laplacian_eigen(FR,2,15,2);
                    anglelem=atan2(Y(:,2),Y(:,1));
                    [alfa,sorting_lem]=sort(anglelem,'descend');
                    aux=find(sorting_lem==ini);
                    sorting_lem_ini=circshift(sorting_lem,-(aux-1));
                    aux_w=find(sorting_tsne==ini_w);
                    sorting_lem_ini_w=circshift(sorting_lem,-(aux_w-1));
                    figure
                    spy(X(flip(sorting_lem_ini),:),'k')
                    pbaspect([23,2,1])
                    
                    %ISOMAP
                   
                    [Y_iso] = isomap(FR,2,15);
                    angleiso=atan2(Y_iso(:,2),Y_iso(:,1));
                    [alfa,sorting_iso]=sort(angleiso,'descend');
                    aux=find(sorting_iso==ini);
                    sorting_iso_ini=circshift(sorting_iso,-(aux-1));
                    aux_w=find(sorting_iso==ini_w);
                    sorting_iso_ini_w=circshift(sorting_iso,-(aux_w-1));
                                        
                    figure
                    spy((X(sorting_iso,:)),'k');
                    pbaspect([23,2,1])
                    
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
                    
                    figure
                    spy((X(sorting_umap,:)),'k');
                    pbaspect([23,2,1])
%                     figure
%                     scatter(reduction_m(:,1),reduction_m(:,2),[],[ 0.2588    0.2588    0.2588],'filled')
%                     alpha 0.7
%                     set(gca,'fontsize',16)
%                     axis([-30 30 -30 30])
%                     ylabel('UMAP - Dim 2')
%                     xlabel('UMAP - Dim 1')
%                     xticks([-30 0 30])
%                     yticks([-30 0 30])
%                     set(gca,'fontsize',18)
                    
                    
                    %                     figure
                    %                     scatter(Y_iso(:,1),Y_iso(:,2),[],[ 0.2588    0.2588    0.2588],'filled')
                    %                     alpha 0.7
                    %                     set(gca,'fontsize',16)
                    %                     axis([-0.03 0.03 -0.03 0.03])
                    %
                    
                    % % % %                     %PCA on circularly shuffled data
                    % % % %                     spikes_sh=circular_shuffling(spikes_d);
                    % % % %                     [coeffsh]=pca(spikes_sh');
                    % % % %                     anglesh=atan2(coeffsh(:,2),coeffsh(:,1));
                    % % % %                     [~,sorting_sh_pca]=sort(anglesh,'ascend');
                    % % % %
                    % % % %
                    % % % %                     %PCA on regularly shuffled data
                    % % % %                     spikes_reg_sh=shuffle(spikes_d')';
                    % % % %                     mat_sh=shuffle(spikes_d')';
                    % % % %                     [coeffsh]=pca(spikes_reg_sh');
                    % % % %                     anglesh=atan2(coeffsh(:,2),coeffsh(:,1));
                    % % % %                     [a1,sorting_reg_sh_pca]=sort(anglesh,'ascend');
                    
                    rmpath("C:\Users\xscogno\MATLAB\drtoolbox.tar");

                    
                end
            end
        end
    end
end

%% Figures of raster plots

% PCA
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
sorting_pca=flip(sorting_pca);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_descend(i)),:),5,'k','filled')
    alpha 0.2
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
ini=sorting_ascend(1);
set(gcf, 'Renderer', 'opengl');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - PCA sorting.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - PCA sorting.svg');
% close all

% PCA - smoothed
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
sorting_w=flip(sorting_w);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_w(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
set(gcf, 'Renderer', 'opengl');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - PCA_w sorting.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - PCA_w sorting.svg');
% close all

%corr
% sorting_corr_ini=flip(sorting_corr_ini);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_corr_ini(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
set(gcf, 'Renderer', 'opengl');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - Corr sorting.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - Corr sorting.svg');
% close all


%corr - smoothed
sorting_corr_ini_w=flip(sorting_corr_ini_w);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_corr_ini_w(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 400])
set(gcf, 'Renderer', 'opengl');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - Corr_w sorting.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - Corr_w sorting.svg');
% close all


%tsne
% sorting_tsne_ini=flip(sorting_tsne_ini);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_tsne_ini(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
set(gcf, 'Renderer', 'opengl');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - tsne sorting.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - tsne sorting.svg');
close all


%tsne - smoothed
% sorting_tsne_ini_w=flip(sorting_tsne_ini_w);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_tsne_ini_w(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
set(gcf, 'Renderer', 'opengl');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - tsne_w sorting.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - tsne_w sorting.svg');
close all


%lem
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_lem_ini(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
set(gcf, 'Renderer', 'opengl');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - lem sorting.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - lem sorting.svg');
close all


%lem - smoothed
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_lem_ini_w(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
set(gcf, 'Renderer', 'opengl');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - lem_w sorting.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - lem_w sorting.svg');
close all


%isomap
% sorting_iso_ini=flip(sorting_iso_ini);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_iso_ini(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
set(gcf, 'Renderer', 'opengl');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - iso sorting.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - iso sorting.svg');
close all


%isomap - smoothed
% sorting_iso_ini_w=flip(sorting_iso_ini_w);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_iso_ini_w(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
set(gcf, 'Renderer', 'opengl');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - iso_w sorting.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - iso_w sorting.svg');
close all


%umap
% sorting_umap_ini=flip(sorting_umap_ini);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_umap_ini(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
set(gcf, 'Renderer', 'opengl');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - umap sorting.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - umap sorting.svg');
close all


%umap - smoothed
% sorting_umap_ini_w=flip(sorting_umap_ini_w);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_umap_ini_w(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 500])
set(gca,'fontsize',22)
set(gca,'fontsize',18)
yticks([100 200])
set(gcf, 'Renderer', 'opengl');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - umap_w sorting.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - umap_w sorting.svg');
close all

%% PCA maps the correlation structure - SORTING

[coeff,~,~] = pca(spikes_d');
angle=atan2(coeff(:,2),coeff(:,1));

[coeff_w,~,~] = pca(FR');
angle_w=atan2(coeff_w(:,2),coeff_w(:,1));
    
for i=1:N
    for j=1:N
        
        if i==j
            delta_ring_pca(i,j)=NaN;
            delta_ring_pca_w(i,j)=NaN;
        else            
            delta_ring_pca(i,j)=(angdiff(angle(i),angle(j)));
            delta_ring_pca_w(i,j)=(angdiff(angle_w(i),angle_w(j)));
        end
    end
end


for n=1:N
    for j=1:N
        %Distance for PCA
        s_n=find(sorting_pca==n);
        s_j=find(sorting_pca==j);
        dist=abs(s_n-s_j);
        if dist>N/2
            Dist_pca(n,j)=N-dist;
        else
            Dist_pca(n,j)=dist;
        end
        
        if n==j
            Dist_pca(n,j)=NaN;
        end

         %Distance for PCA smoothed
        s_n=find(sorting_w==n);
        s_j=find(sorting_w==j);
        dist=abs(s_n-s_j);
        if dist>N/2
            Dist_pca_w(n,j)=N-dist;
        else
            Dist_pca_w(n,j)=dist;
        end
        
        if n==j
            Dist_pca_w(n,j)=NaN;
        end
        
        %Distance for xcorr
        s_n=find(sorting_corr==n);
        s_j=find(sorting_corr==j);
        dist=abs(s_n-s_j);
        if dist>N/2
            Dist_corr(n,j)=N-dist;
        else
            Dist_corr(n,j)=dist;
        end
        
        if n==j
            Dist_corr(n,j)=NaN;
        end
    end
end

Dist_corr_vec=Dist_corr(:);
Dist_pca_vec=Dist_pca(:);
Dist_pca_vec_w=Dist_pca_w(:);
Dist_ring_pca_vec=delta_ring_pca(:);
Dist_ring_pca_w_vec=delta_ring_pca_w(:);

Delta_Dist=(Dist_pca_vec-Dist_corr_vec).*(Dist_pca_vec-Dist_corr_vec);
alfa_b=nanmean(Delta_Dist);
std_b=nanstd(Delta_Dist);
% figure
% hist(Delta_Dist)
drop=randperm(N*N,floor(0.3*N*N));

figure
scatter(Dist_corr_vec(drop),Dist_pca_vec(drop),3,'o','filled');
alpha 0.2
axis([0 floor(N/2) 0 floor(N/2)]);
xticks([0 floor(N/2)])
yticks([0 floor(N/2)])
axis square
set(gca,'fontsize',18)
ylabel({'Distance between cells';'PCA (cells #)'});
xlabel({'Distance between cells';'XCORR (cells #)'});

figure
scatter(Dist_corr_vec(drop),Dist_pca_vec_w(drop),3,'o','filled');
alpha 0.2
axis([0 floor(N/2) 0 floor(N/2)]);
xticks([0 floor(N/2)])
yticks([0 floor(N/2)])
axis square
set(gca,'fontsize',18)
ylabel({'Distance between cells';'PCA_W (cells #)'});
xlabel({'Distance between cells';'XCORR (cells #)'});

figure
scatter(Dist_corr_vec(drop),Dist_ring_pca_vec(drop),3,'o','filled');
alpha 0.2
axis([0 floor(N/2) -pi pi]);
xticks([0 floor(N/2)])
yticks([0 floor(N/2)])
axis square
set(gca,'fontsize',18)
ylabel({'Distance between cells';'Ring (cells #)'});
xlabel({'Distance between cells';'XCORR (cells #)'});

figure
scatter(Dist_corr_vec(drop),Dist_ring_pca_w_vec(drop),3,'o','filled');
alpha 0.2
axis([0 floor(N/2) -pi pi]);
xticks([0 floor(N/2)])
yticks([0 floor(N/2)])
axis square
set(gca,'fontsize',18)
ylabel({'Distance between cells';'Ring_W (cells #)'});
xlabel({'Distance between cells';'XCORR (cells #)'});

for sh=1:200
    sorting_sh=randperm(N);
    for n=1:N
        for j=1:N           
            s_n=find(sorting_sh==n);
            s_j=find(sorting_sh==j);
            dist=abs(s_n-s_j);
            if dist>N/2
                Dist_pca_sh(n,j,sh)=N-(abs(s_n-s_j));
            else
                Dist_pca_sh(n,j,sh)=dist;
            end
            
            if n==j
               Dist_pca_sh(n,j,sh)=NaN;
            end
        end
    end
end

for sh=1:200
    mat=Dist_pca_sh(:,:,sh);
    Dist_pca_sh_vec(:,sh)=mat(:);
    Delta_Dist_sh(:,sh)=(Dist_corr_sh_vec(:,sh)-Dist_pca_vec).*(Dist_corr_sh_vec(:,sh)-Dist_pca_vec);    
end


figure
scatter(Dist_corr_vec(drop),Dist_pca_sh_vec(drop,2),3,'o','filled');
alpha 0.2
axis([0 floor(N/2) 0 floor(N/2)]);
xticks([0 floor(N/2)])
yticks([0 floor(N/2)])
axis square
set(gca,'fontsize',18)
ylabel({'Distance between cells';'PCA_W - Shuffle (cells #)'});
xlabel({'Distance between cells';'XCORR (cells #)'});
xticks([0 240])

% figure
% scatter(Dist_corr_sh_vec(drop,2),Dist_pca_vec_w(drop),3,'o','filled');
% alpha 0.2
% axis([0 floor(N/2) 0 floor(N/2)]);
% xticks([0 floor(N/2)])
% yticks([0 floor(N/2)])
% axis square
% set(gca,'fontsize',18)
% ylabel({'Distance between cells';'PCA_W (cells #)'});
% xlabel({'Distance between cells';'XCORR - SH (cells #)'});

MSE_data=nansum((Dist_pca_vec_w-Dist_corr_vec).*(Dist_pca_vec_w-Dist_corr_vec))./(N*N - N);

for sh=1:200
    MSE_sh(sh)=nansum((Dist_pca_sh_vec(:,sh)-Dist_corr_vec).*(Dist_pca_sh_vec(:,sh)-Dist_corr_vec))./(N*N - N);
end
    

figure
y1=histogram(MSE_sh,20);
ylabel('Counts')
box off
set(gca,'fontsize',18)
alpha 0.3
hold on;
h=xline(MSE_data,'--k','linewidth',1.5);
axis([5000 10000 0 40]);
xticks([6000 9000])
yticks([0 20 40]);
xlabel('MSE')

figure
y1=histogram(MSE_sh,20);
ylabel('Counts')
box off
set(gca,'fontsize',18)
alpha 0.3
axis([9500 10000 0 40]);
xticks([9500 10000])
yticks([0 20 40]);
xlabel('MSE')

figure
hold on;
h=xline(MSE_data,'--k','linewidth',1.5);
axis([5000 6000 0 40]);
xticks([5000 6000])
yticks([0 20 40]);
xlabel('MSE')
ylabel('Counts')
set(gca,'fontsize',18)

%% PCA maps the correlation structure - SORTING

Pearson=corr(FR',FR');
Dist_ring_pca_vec=delta_ring_pca(:);
Dist_ring_pca_w_vec=delta_ring_pca_w(:);

figure
plot(Pearson(:),Dist_ring_pca_w_vec(:),'*')



%% PC loadings and correlations

X = spikes_downsample(spikes_d,N,downsampling_factor);
FRp = spikes_downsample(spikes_d,N,downsampling_factor);
for i=1:N
    FR(i,:)=full(fire_rate(FRp(i,:),29,'g')); %Smoothing in about 10 seconds
end


% Angle PCA
[coeff,score]=pca(spikes_d');
[coeff2]=pca(FR');


Pearson=corr(FR',FR');
for i=1:N
    for j=1:N
        prod_pca(i,j)=coeff2(i)*coeff2(j);
    end
end

figure
plot(Pearson(:),prod_pca(:),'*')
ylabel('LoadingPC1*LoadingPC2');
xlabel('Pearson Correlation');
set(gca,'fontsize',16);
xl=xline(0);   
yl=yline(0);    

%% Figures plane

aux1= find(coeff(:,1)<=-0.093601);
aux2= find(coeff(:,1)>-0.093602);
cell1=intersect(aux1,aux2);

aux1= find(coeff(:,1)>=0.05907);
aux2= find(coeff(:,1)<0.05908);
cell2=intersect(aux1,aux2);

figure
scatter(coeff(:,1),coeff(:,2),60,'MarkerEdgeColor',[ 85    85    85]./255,'MarkerFaceColor',[ 160    160    160]./255);
alpha 0.6
hold on
set(gca,'fontsize',16)
axis([-0.13 0.13 -0.13 0.13])
ylabel('PC2')
xlabel('PC1')
xticks([-0.1 0 0.1])
yticks([-0.1 0 0.1])
set(gca,'fontsize',18)
axis square


figure
scatter(coeff(:,1),coeff(:,2),60,'MarkerEdgeColor',[ 85    85    85]./255,'MarkerFaceColor',[ 160    160    160]./255);
alpha 0.6
hold on
scatter(coeff(cell1,1),coeff(cell1,2),70,'MarkerEdgeColor',[ 0   204    102]./255,'MarkerFaceColor',[21   255    153]./255);
scatter(coeff(cell2,1),coeff(cell2,2),70,'MarkerEdgeColor',[ 255   128    0]./255,'MarkerFaceColor',[255   178    102]./255);
set(gca,'fontsize',16)
axis([-0.13 0.13 -0.13 0.13])
ylabel('PC2')
xlabel('PC1')
xticks([-0.1 0 0.1])
yticks([-0.1 0 0.1])
set(gca,'fontsize',18)
axis square

figure
scatter(0.1*sin(angle2),0.1*cos(angle2),60,'MarkerEdgeColor',[ 85    85    85]./255,'MarkerFaceColor',[ 160    160    160]./255);
alpha 0.3
set(gca,'fontsize',16)
axis([-0.13 0.13 -0.13 0.13])
ylabel('PC2')
xlabel('PC1')
xticks([-0.1 0 0.1])
yticks([-0.1 0 0.1])
set(gca,'fontsize',18)
axis square

figure
scatter(0.1*sin(angle2),0.1*cos(angle2),60,'MarkerEdgeColor',[ 85    85    85]./255,'MarkerFaceColor',[ 160    160    160]./255);
alpha 0.3
hold on
scatter(-0.1*sin(angle2(cell1)),-0.1*cos(angle2(cell1)),70,'MarkerEdgeColor',[ 0   204    102]./255,'MarkerFaceColor',[21   255    153]./255);
scatter(0.1*sin(angle2(cell2)),0.1*cos(angle2(cell2)),70,'MarkerEdgeColor',[ 255   128    0]./255,'MarkerFaceColor',[255   178    102]./255);
set(gca,'fontsize',16)
axis([-0.13 0.13 -0.13 0.13])
ylabel('PC2')
xlabel('PC1')
xticks([-0.1 0 0.1])
yticks([-0.1 0 0.1])
set(gca,'fontsize',18)
axis square



figure
scatter(Y_tsne(:,1),Y_tsne(:,2),[],[ 0.2588    0.2588    0.2588],'filled')
alpha 0.7
set(gca,'fontsize',16)
axis([-100 100 -100 100])
ylabel('tSNE - Dim 2')
xlabel('tSNE - Dim 1')
xticks([-100 0 100])
yticks([-100 0 100])
set(gca,'fontsize',18)

figure
scatter(Y_lem(:,1),Y_lem(:,2),[],[ 0.2588    0.2588    0.2588],'filled')
alpha 0.7
set(gca,'fontsize',16)
axis([-0.05 0.02 -0.03 0.04])
ylabel('LEM - Dim 2')
xlabel('LEM - Dim 1')
xticks([-0.05 0 0.02])
yticks([-0.03 0 0.04])
set(gca,'fontsize',18)

figure
scatter(Y_iso(:,1),Y_iso(:,2),[],[ 0.2588    0.2588    0.2588],'filled')
alpha 0.7
set(gca,'fontsize',16)
axis([-30 30 -30 30])
ylabel('Isomap - Dim 2')
xlabel('Isomap - Dim 1')
xticks([-30 0 30])
yticks([-30 0 30])
set(gca,'fontsize',18)

figure
scatter(reduction_m(:,1),reduction_m(:,2),[],[ 0.2588    0.2588    0.2588],'filled')
alpha 0.7
set(gca,'fontsize',16)
axis([-5 5 -5 5])
ylabel('UMAP - Dim 2')
xlabel('UMAP - Dim 1')
xticks([-5 0 5])
yticks([-5 0 5])
set(gca,'fontsize',18)

%% Figures plane colouring according ensembles
cc=parula(10);
sorting=sorting_pca;
cells_per_ens=floor(N/10);
for e=1:10
    ens(e,:)=sorting((e-1)*cells_per_ens+1 : e*cells_per_ens);   
    ense_n(e,:)=e*ones(1,cells_per_ens);
end

Ens(:,1)=reshape(ens,cells_per_ens*10,1);
Ens(:,2)=reshape(ense_n,cells_per_ens*10,1);


figure
hold on
for e=1:10
    aux=find(Ens(:,2)==e);
    cells=Ens(aux,1);
scatter(coeff(cells,1),coeff(cells,2),[],cc(e,:),'filled');
alpha 0.7
set(gca,'fontsize',16)
axis([-0.2 0.2 -0.2 0.2])
ylabel('PCA - PC 2')
xlabel('PCA - PC 1')
xticks([-0.2 0 0.2])
yticks([-0.2 0 0.2])
set(gca,'fontsize',18)
end
% colormap(parula(10));    
% colorbar 

figure
hold on
for e=1:10
    aux=find(Ens(:,2)==e);
    cells=Ens(aux,1);
    scatter(Y_tsne(cells,1),Y_tsne(cells,2),[],cc(e,:),'filled')
    alpha 0.7
    set(gca,'fontsize',16)
    axis([-100 100 -100 100])
    ylabel('tSNE - Dim 2')
    xlabel('tSNE - Dim 1')
    xticks([-100 0 100])
    yticks([-100 0 100])
    set(gca,'fontsize',18)
end

figure
hold on
for e=1:10
    aux=find(Ens(:,2)==e);
    cells=Ens(aux,1);
    scatter(Y_lem(cells,1),Y_lem(cells,2),[],cc(e,:),'filled')
    alpha 0.7
    set(gca,'fontsize',16)
    axis([-0.05 0.02 -0.03 0.04])
    ylabel('LEM - Dim 2')
    xlabel('LEM - Dim 1')
    xticks([-0.05 0 0.02])
    yticks([-0.03 0 0.04])
    set(gca,'fontsize',18)
end

figure
hold on
for e=1:10
    aux=find(Ens(:,2)==e);
    cells=Ens(aux,1);
    scatter(Y_iso(cells,1),Y_iso(cells,2),[],cc(e,:),'filled')
    alpha 0.7
    set(gca,'fontsize',16)
    axis([-30 30 -30 30])
    ylabel('Isomap - Dim 2')
    xlabel('Isomap - Dim 1')
    xticks([-30 0 30])
    yticks([-30 0 30])
    set(gca,'fontsize',18)
end
colormap(parula(10))
colorbar

figure
hold on
for e=1:10
    aux=find(Ens(:,2)==e);
    cells=Ens(aux,1);
    scatter(reduction_m(cells,1),reduction_m(cells,2),[],cc(e,:),'filled')
    alpha 0.7
    set(gca,'fontsize',16)
    axis([-5 5 -5 5])
    ylabel('UMAP - Dim 2')
    xlabel('UMAP - Dim 1')
    xticks([-5 0 5])
    yticks([-5 0 5])
    set(gca,'fontsize',18)
end
colormap(parula(10))
colorbar