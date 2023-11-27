% Analysis of scatter plot of delta preferred phase VS pairwise anatomical
% distance
clear all
close all
rec_data_path='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
data=load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\locking_all_sessions.mat');
pixel_size_new=1.18185;
pixel_size_old=1.78211;
countsc=0;
individual_plots=0;
correlation_sh_pooled=[];
N_sh=100;
for w=1:length(data.waves)
    row_w=data.waves(w);
    disp(w)
    countsc=countsc+1;
    mouse=['L',num2str(data.big_table(row_w,1)),'M',num2str(data.big_table(row_w,2))];
    day=data.big_table(row_w,3);
    s=data.big_table(row_w,4);
    munit=data.big_table(row_w,5);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
%     file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_spk=[data.dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_anat=[data.dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

    load(file_name_anat,'-mat'); %Anatomical information
    load(file_name_spk,'-mat'); %Spike times
    spikes=full(spikes_d_s);
    [N,T]=size(spikes);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mean phases and locked cells
    not_locked=data.not_locked_all_sessions{w};
    locked=data.locked_all_sessions{w};
    mean_p=data.mean_p_all_sessions{w};
%     mean_p_locked=mean_p(locked);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distance in the tissue
    if w>10; pixel_size=pixel_size_old;
    else; pixel_size=pixel_size_new; end

    for i=1:N
        r_i(i,:)=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
        dist_origin(i)=norm(r_i(i,:));
    end
%     r_i(not_locked,:)=[];
    dist_origin_locked=dist_origin(locked);
    max_x=max(r_i(:,1));
    max_y=max(r_i(:,2));
    max_xy=max(max_x,max_y);
    edges_x=0:50:ceil(max_xy)+50;

    count=0;
    for i=1:N
        for j=i+1:N
            count=count+1;
            delta_tissue(i,j)=norm(r_i(j,:)-r_i(i,:));
            delta_tissue(j,i)= delta_tissue(i,j);
            delta_tissue_vec(count)=delta_tissue(i,j);
        end
    end
    delta_tissue(not_locked,:)=[];
    delta_tissue(:,not_locked)=[];
    dist_origin_locked=dist_origin(locked);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distance in preferred
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% phase
    count=0;
    for i=1:N
        for j=i+1:N
            count=count+1;
            delta_phase_mean(i,j)=angdiff(mean_p(i),mean_p(j));
            delta_phase_mean(j,i)=delta_phase_mean(i,j);
            delta_phase_vec(count)=delta_phase_mean(i,j);
        end
    end
    delta_phase_mean(not_locked,:)=[];
    delta_phase_mean(:,not_locked)=[];

    for sh=1:N_sh
        count=0;
        clear mean_p_sh delta_phase_vec_sh; 
        mean_p_sh=mean_p(randperm(N));
        for i=1:N
            for j=i+1:N
                count=count+1;
                delta_phase_vec_sh(count)=angdiff(mean_p_sh(i),mean_p_sh(j));
            end
        end
        [rho_sh(w,sh),pval_sh(w,sh)] = corr(delta_tissue_vec',(delta_phase_vec_sh)');
        [rho_sh_s(w,sh),pval_sh_s(w,sh)] = corr(delta_tissue_vec',(delta_phase_vec_sh)','type','Spearman');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save matrices
    matrix_distance_tissue{w}=delta_tissue_vec;
    matrix_distance_phase{w}=delta_phase_vec;
%     matrix_distance_tissue_sh{w}=repmat(delta_tissue_vec,N_sh,1);
%     matrix_distance_phase_sh{w}=delta_phase_vec_sh;
    correlation_sh{w}.rho=rho_sh;
    correlation_sh{w}.pval=pval_sh;
%     correlation_sh_pooled=[correlation_sh_pooled;[w*ones(N_sh,1),rho_sh',pval_sh']];

    [rho(w),pval(w)] = corr(delta_tissue_vec',(delta_phase_vec)');
    [rho_s(w),pval_s(w)] = corr(delta_tissue_vec',(delta_phase_vec)','type','Spearman');

    %     diff_phase_reshape=delta_phase_mean(:);
    %     diff_pos_reshape=delta_tissue(:);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation of traveling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wave - Distance in
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tissue and phase
%     [~,idx_p]=sort(mean_p,'ascend');
%     [~,idx_dist]=sort(dist_origin,'ascend');
%     [~,sorting_descend,~] = get_sorting(spikes);
%     for i=1:N; r_i_tw(sorting_descend(i),:)=r_i(idx_dist(i),:); end   
%     pp=discretize(mean_p,[-pi:2*pi/30:pi]);
%     cc=jet(30);
%     figure
%     subplot(1,2,1)
%     hold on
%     for n=1:N
%         scatter(r_i(n,1),r_i(n,2),[],cc(pp(n),:),"filled");
%     end    
%     subplot(1,2,2)
%     hold on
%     for n=1:N
%         scatter(r_i_tw(n,1),r_i_tw(n,2),[],cc(pp(n),:),"filled");
%     end 
%     count=0;
%     for i=1:N
%         for j=i+1:N
%             count=count+1;
%             delta_tissue_tw(i,j)=norm(r_i_tw(j,:)-r_i_tw(i,:));
%             delta_tissue_tw(j,i)= delta_tissue_tw(i,j);
%             delta_tissue_vec_tw(count)=delta_tissue_tw(i,j);
%         end
%     end
%     matrix_distance_tissue_tw{w}=delta_tissue_vec_tw;

    clear Anat cells_d count dates delta_phase_vec delta_phase_mean delta_tissue delta_tissue_vec dist_origin dist_origin_locked edges_x locked mean_p ...
        r_i spikes spikes_d_s i j not_locked delta_phase_vec_sh mean_p_sh  ...
       
end

%% Quantifications for phases after pooling all cycles
min_r=min(rho);
max_r=max(rho);
min_r_sh=min(rho_sh(:));
max_r_sh=max(rho_sh(:));

% Range of Pearson
for w=1:15
    min_r_sh_perw(w)=min(rho_sh(w,:));
    max_r_sh_perw(w)=max(rho_sh(w,:));
end

% Significance testing Pearson
rho_abs=abs(rho);
rho_sh_abs=abs(rho_sh);
for i=1:15
    cutoff(i)=prctile(rho_sh_abs(i,:),95);
end
aux=find(rho_abs>cutoff);

pout=myBinomTest(length(aux),15,0.05,'one');

% Significance testing Spearman
rho_abs=abs(rho_s);
rho_sh_abs=abs(rho_sh_s);
for i=1:15
    cutoff_s(i)=prctile(rho_sh_abs(i,:),95);
end
aux=find(rho_abs>cutoff_s);


%% FIGURES corresponding to cell 1
edges_for_boxplot=0:150:1200;

%One separate figure per session
for w=1:15   
    vect_d=matrix_distance_tissue{1,w};
    vect_p=matrix_distance_phase{1,w};
    vect_d_dis=discretize(vect_d,edges_for_boxplot);

    count_bp=zeros(1,length(edges_for_boxplot));
    mat=nan(length(vect_p),length(edges_for_boxplot));
    for i=1:length(vect_d_dis)
        col=vect_d_dis(i);
        count_bp(col)=count_bp(col)+1;
        mat(count_bp(col),col)=vect_p(i);
    end
    num_bins=find(count_bp==0,1);

    %Scatter plot
    figure
%     scatter(vect_d,vect_p,'k.');
      scatterhist(vect_d,vect_p, 'Marker','.','Color','k');
%     if mod(w,4)==1 %         ylabel('Data'); 
%         ylabel({'\Delta preferred';' phase (rad)'});   %Delta preferred phase
%     else
%         ylabel(' '); % dummy y-axis label
%     end
    ylabel({'Difference in';'preferred phase (rad)'});   %Delta preferred phase
    xlabel('Pairwise anatomical distance (um)');
    set(gca,'fontsize',20,'YColor','k','XColor','k');
    axis([0 inf -3.2 3.2]);
    yticks([-pi 0 pi]);
    yticklabels({'-\pi','0','\pi'});
    box off
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_ScatterHistPlot_session',num2str(w),'.png'));
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_ScatterHistPlot_session',num2str(w),'.fig'));
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_ScatterHistPlot_session',num2str(w),'.svg'));

    %Box plot
%     figure
%     boxplot(mat)
%     ylabel({'Difference in';'preferred phase (rad)'});   %Delta preferred phase
%     xlabel('Pairwise anatomical distance bin');
%     set(gca,'fontsize',20,'YColor','k','XColor','k');
%     axis([0 num_bins -3.2 3.2]);
%     xticks(1:num_bins-1);
%     yticks([-pi 0 pi]);
%     yticklabels({'-\pi','0','\pi'});
%     box off
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_BoxPlot_session',num2str(w),'.png'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_BoxPlot_session',num2str(w),'.fig'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_BoxPlot_session',num2str(w),'.svg'));    
    clear mat_d mat_p mat_d2 mat_p2 vect_d vect_p vect_d_dis count_bp mat
    close all
end

%Subplot with scatter plot
figure
for w=1:15   
    vect_d=matrix_distance_tissue{1,w};
    vect_p=matrix_distance_phase{1,w};
    vect_d_dis=discretize(vect_d,edges_for_boxplot);
    count_bp=zeros(1,length(edges_for_boxplot));
    mat=nan(length(vect_p),length(edges_for_boxplot));
    for i=1:length(vect_d_dis)
        col=vect_d_dis(i);
        count_bp(col)=count_bp(col)+1;
        mat(count_bp(col),col)=vect_p(i);
    end
    num_bins=find(count_bp==0,1);
    %Scatter plot
    subplot(4,4,w)
    scatter(vect_d,vect_p,'k.');
    if mod(w,4)==1 %         ylabel('Data'); 
        ylabel({'\Delta preferred';' phase (rad)'});   %Delta preferred phase
    else
        ylabel(' '); % dummy y-axis label
    end
    xlabel('Anatomical distance (\mum)');
    set(gca,'fontsize',16,'YColor','k','XColor','k');
    axis([0 1000 -3.2 3.2]);
    yticks([-pi 0 pi]);
    title(['Session ',num2str(w)]);
    yticklabels({'-\pi','0','\pi'});
    box off
%     [rho(w),pval(w)] = corr(vect_d',abs(vect_p)');   
    clear mat_d mat_p mat_d2 mat_p2 vect_d vect_p vect_d_dis count_bp mat
%     close all
end

% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_ScatterPlot_Allsessions.png'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_ScatterPlot_Allsessions.fig'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_ScatterPlot_Allsessions.svg'));


%Scatterhist for scatter plot
figure
for w=1:15   
    vect_d=matrix_distance_tissue{1,w};
    vect_p=matrix_distance_phase{1,w};
    vect_d_dis=discretize(vect_d,edges_for_boxplot);
    count_bp=zeros(1,length(edges_for_boxplot));
    mat=nan(length(vect_p),length(edges_for_boxplot));
    for i=1:length(vect_d_dis)
        col=vect_d_dis(i);
        count_bp(col)=count_bp(col)+1;
        mat(count_bp(col),col)=vect_p(i);
    end
    num_bins=find(count_bp==0,1);
    %Scatter plot
    subplot(4,4,w)
    scatterhist(vect_d,vect_p, 'Marker','.','Color','k');
    if mod(w,4)==1 %         ylabel('Data'); 
        ylabel({'\Delta preferred';' phase (rad)'});   %Delta preferred phase
    else
        ylabel(' '); % dummy y-axis label
    end
    xlabel('Anatomical distance (\mum)');
    set(gca,'fontsize',16,'YColor','k','XColor','k');
    axis([0 1000 -3.2 3.2]);
    yticks([-pi 0 pi]);
    title(['Session ',num2str(w)]);
    yticklabels({'-\pi','0','\pi'});
    box off
%     [rho(w),pval(w)] = corr(vect_d',abs(vect_p)');   
    clear mat_d mat_p mat_d2 mat_p2 vect_d vect_p vect_d_dis count_bp mat
%     close all
end

% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_ScatterPlot_Allsessions.png'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_ScatterPlot_Allsessions.fig'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_ScatterPlot_Allsessions.svg'));


%Figure of subplot with box plot
fig=figure;
% t = tiledlayout(4,4);
hold on
for w=1:15       
    vect_d=matrix_distance_tissue{1,w};
    vect_p=matrix_distance_phase{1,w};
    vect_d_dis=discretize(vect_d,edges_for_boxplot);
    count_bp=zeros(1,length(edges_for_boxplot));
    mat=nan(length(vect_p),length(edges_for_boxplot));
    for i=1:length(vect_d_dis)
        col=vect_d_dis(i);
        count_bp(col)=count_bp(col)+1;
        mat(count_bp(col),col)=vect_p(i);
    end
    num_bins=find(count_bp==0,1);

    %Box plot
    subplot(4,4,w)
    boxplot(mat)    
%     if mod(w,4)==1 %         ylabel('Data'); 
%         ylabel({'\Delta preferred';' phase (rad)'});   %Delta preferred phase
%     else
%         ylabel(' '); % dummy y-axis label
%     end    
%     xlabel('Anatomical distance (\mum)');
    set(gca,'fontsize',12,'YColor','k','XColor','k');
    axis([0 num_bins -3.2 3.2]);
    xticks(1:num_bins-1);
    yticks([-pi 0 pi]);
    yticklabels({'-\pi','0','\pi'});
    title(['Session ',num2str(w)]);
    box off
%     axis square
    clear mat_d mat_p mat_d2 mat_p2 vect_d vect_p vect_d_dis count_bp mat
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,{'\Delta preferred';' phase (rad)'});
xlabel(han,'Anatomical distance (bin #)');
set(gca,'fontsize',16,'YColor','k','XColor','k');

% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_BoxPlot_Allsessions.png'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_BoxPlot_Allsessions.fig'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_BoxPlot_Allsessions.svg'));
%     
% close all



%Heatmaps
edges_phase=-pi:2*pi/80:pi;
edges_dist=0:20:1000;
figure
for w=1:15
    vect_d=matrix_distance_tissue{1,w};
    vect_p=matrix_distance_phase{1,w};
    %     histogram2(X,'CdataMode','auto','Edges',{edges_dist edges_phase},'normalization','probability')
    
    subplot(4,4,w)
    histogram2(vect_d',vect_p',edges_dist,edges_phase,'normalization','probability','DisplayStyle','tile','ShowEmptyBins','on');
    ylabel({'\Delta preferred';' phase (rad)'});
    xlabel('Anatomical distance (\mum)');
%     colorbar
    view(2)
    axis([min(edges_dist) max(edges_dist) min(edges_phase) max(edges_phase)]);
    yticks([-pi 0 pi]);
    yticklabels({'-\pi','0','\pi'});
    xticks([0 500 1000]);
    axis square
    colormap pink
    box off
    set(gca,'fontsize',14)
    caxis([0 0.001])
    clear mat_d mat_p mat_d2 mat_p2 vect_d vect_p vect_d_dis count_bp mat
%     close all
end

saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_HeatMap_colorbar_Allsessions.png'));
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_HeatMap_colorbar_Allsessions.fig'));
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_HeatMap_colorbar_Allsessions.svg'));
    
close all


%% Analysis for each cycle of the oscillation separately - Phase
clear all
close all

dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
plot_all=1;
[big_table, waves] = get_big_table();
pixel_size_new=1.18185;
pixel_size_old=1.78211;
count_s=0;
correlation_sh_pooled=[];
correlation_abs_sh_pooled=[];
N_sh=100;

for w=1:length(waves)
    row_w=waves(w);
    disp(w);
    %     count=count+1;
    mouse=['L',num2str(big_table(row_w,1)),'M',num2str(big_table(row_w,2))];
    day=big_table(row_w,3);
    s=big_table(row_w,4);
    munit=big_table(row_w,5);
    clus=10;
    disc_phase=10;
    %     mean_p=data.mean_p_all_sessions{w};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
    load([dpath,strcat('recording_dates_',mouse,'.mat')]);
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_dff=[dpath ['DFF_120ms_Do_SNRH','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

    load(file_name_anat,'-mat'); %Anatomical information
    load(file_name_spk,'-mat'); %Spike times
    spikes_d=full(spikes_d_s);
    [N,~]=size(spikes_d);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Anatomical position of each cell
    if w>10; pixel_size=pixel_size_old; else; pixel_size=pixel_size_new; end
    for i=1:size(spikes_d,1)
        r_i(i,:)=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
        %         dist_origin(i)=norm(r_i(i,:));
    end
    delta_tissue=nan(N,N);
    count=0;
    for i=1:N
        for j=i+1:N
            count=count+1;
            delta_tissue(i,j)=norm(r_i(j,:)-r_i(i,:));
            delta_tissue(j,i)= delta_tissue(i,j);
            delta_tissue_vec(count)=delta_tissue(i,j);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Condition on having waves and phase per seq
    dt=big_table(row_w,8); if ~isinteger(dt); dt=floor(dt); end
    clear FRp;
    for i=1:N
        FRp(i,:)=full(fire_rate(spikes_d(i,:),1*dt,'g')); %smooth using as kernel the dt chosen for each session
    end
    [coefft,scoret,~] = pca(FRp');
    phase_f=(atan2(smooth(scoret(:,2),floor(1*dt)),smooth(scoret(:,1),floor(1*dt))));
    radius_f=sqrt(coefft(:,2).*coefft(:,2)+coefft(:,1).*coefft(:,1));
    num_clus_discr=10;
    make_fig=0;
    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes_d);
    spikes_r=[]; %Reduced spike matrix; only contains wave epochs
    phase_r=[]; %Reduced phase; only contains wave epochs
    for wa=1:size(table_u,1)
        spikes_r=[spikes_r,spikes_d(:,table_u(wa,1):table_u(wa,2))];
        phase_r=[phase_r;phase_f(table_u(wa,1):table_u(wa,2))];
    end
    spikes=spikes_d;
    phase=phase_r;
    spikes_d=[];
    spikes_d=spikes_r; %spikes only during sequences

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Here I calculate the phase per sequence
    for s=1:size(table_u,1)
        clear spikes_seq; spikes_seq=spikes(:,table_u(s,1):table_u(s,2));
        clear phase_seq; phase_seq=phase_f(table_u(s,1):table_u(s,2));

        for i=1:N
            clear p; p=phase_seq(find(spikes_seq(i,:)));
            if length(p)>5; mean_p_seq(i,s)=circ_mean(p); else mean_p_seq(i,s)=nan; end
        end
    end
    %     mean_p_all_sessions{w}=mean_p;
    %     std_p_all_sessions{w}=std_p;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loop on sequences
    tissue_vec{w}=[];
    phase_vec{w}=[];
    for s=1:size(table_u,1)
        count_s=count_s+1;
        mean_p=mean_p_seq(:,s);
        delta_phase_mean=nan(N,N);
        count=0;
        for i=1:N
            for j=i+1:N
                count=count+1;
                delta_phase_mean(i,j)=angdiff(mean_p(i),mean_p(j));
                delta_phase_mean(j,i)=delta_phase_mean(i,j);
                delta_phase_vec(count)=delta_phase_mean(i,j);
            end
        end

        for sh=1:N_sh %Differences in preferred phase in shuffled data
            count=0;
            clear mean_p_sh delta_phase_vec_sh;
            mean_p_sh=mean_p(randperm(N));
            for i=1:N
                for j=i+1:N
                    count=count+1;
                    delta_phase_vec_sh(count)=angdiff(mean_p_sh(i),mean_p_sh(j));
                end
            end
            clear aux; aux=~isnan(delta_phase_vec_sh);
            [rho_sh(count_s,sh),pval_sh(count_s,sh)] = corr(delta_tissue_vec(aux)',(delta_phase_vec_sh(aux))');
            [rho_sh_abs(count_s,sh),pval_sh_abs(count_s,sh)] = corr(delta_tissue_vec(aux)',abs(delta_phase_vec_sh(aux))');

        end

        clear aux; aux=~isnan(delta_phase_vec);
        if (w==8 && s==19)
%             scatterhist(delta_tissue_vec(aux)',(delta_phase_vec(aux))', 'Marker','.','Color','k');
%             ylabel({'\Delta preferred';' phase (rad)'});   %Delta preferred phase
%             xlabel('Anatomical distance (\mum)');
%             set(gca,'fontsize',16,'YColor','k','XColor','k');
%             axis([0 1000 -3.2 3.2]);
%             yticks([-pi 0 pi]);
%             title(['Session ',num2str(w),' - Sequence # ',num2str(s)]);
%             yticklabels({'-\pi','0','\pi'});
%             box off
            [rho_seq19,pval_seq19] = corr(delta_tissue_vec(aux)',(delta_phase_vec(aux))');

            edges_phase=-pi:2*pi/80:pi;
            edges_dist=0:20:1000;
            figure
            histogram2(delta_tissue_vec(aux)',(delta_phase_vec(aux))',edges_dist,edges_phase,'normalization','probability','DisplayStyle','tile','ShowEmptyBins','on');
            ylabel({'\Delta preferred';' phase (rad)'});
            xlabel('Anatomical distance (\mum)');
            %     colorbar
            view(2)
            axis([min(edges_dist) max(edges_dist) min(edges_phase) max(edges_phase)]);
            yticks([-pi 0 pi]);
            yticklabels({'-\pi','0','\pi'});
            xticks([0 500 1000]);
            axis square
            colormap pink
            box off
            set(gca,'fontsize',14)
            caxis([0 0.001])

            
        end

%         correlation_sh_pooled=[correlation_sh_pooled;[count_s*ones(N_sh,1),rho_sh',pval_sh']];
%         correlation_abs_sh_pooled=[correlation_abs_sh_pooled;[count_s*ones(N_sh,1),rho_sh_abs',pval_sh_abs']];

        [rho(count_s),pval(count_s)] = corr(delta_tissue_vec(aux)',(delta_phase_vec(aux))');
        [rho_abs(count_s),pval_abs(count_s)] = corr(delta_tissue_vec(aux)',abs(delta_phase_vec(aux))');
        tissue_vec{w}=[tissue_vec{w};delta_tissue_vec(aux)'];
        phase_vec{w}=[phase_vec{w};delta_phase_vec(aux)'];

        clear mean_p mean_p_sh delta_phase_mean delta_phase_sh  delta_phase_vec aux
    end

    clear phase_f phase_r phase_seq radius_f spikes_seq spikes_r spikes_d_s spikes_d spikes scoret template ...
        phase mean_p_seq FRp r_i r_i_sh table_u sh_ord cells_d delta_tissue delta_phase_mean delta_phase_sh ...
        dist_origin coefft Anat phases_sh delta_tissue_vec aux delta_phase_vec_sh mean_p_sh pval_sh  
end

%PREFERED PHASES
% Smalles and largest correlation
min_r=min(rho);
max_r=max(rho);
min_r_sh=min(rho_sh(:));
max_r_sh=max(rho_sh(:));

% Range of Pearson
for w=1:421
    min_r_sh_perw(w)=min(rho_sh(w,:));
    max_r_sh_perw(w)=max(rho_sh(w,:));
end

% Significance testing Pearson
rho_abs2=abs(rho);
rho_sh_abs2=abs(rho_sh);
for i=1:421
    cutoff(i)=prctile(rho_sh_abs2(i,:),95);
end
aux=find(rho_abs2>cutoff);

pout=myBinomTest(length(aux),421,0.05,'one');

%Figures of distributions
[bar_hexp,bins_exp]=histcounts(rho(:),[-0.1:1/100:0.1],'Normalization','probability');
[bar_hshuffle,bins_shuffle]=histcounts(rho_sh(:),[-0.1:1/200:0.1],'Normalization','probability');

figure
yyaxis left
bar(bins_exp(1:end-1)+0.005,bar_hexp)
yticks([0 0.15 0.3]);
axis([-0.2 0.2 0 0.35])
ylabel('Normalized frequency - Data');
hold on
yyaxis right
plot(bins_shuffle(1:end-1) + (bins_shuffle(end) - bins_shuffle(end-1))/2,bar_hshuffle,'-.','linewidth',3)
axis([-0.1 0.1 0 0.18])
yticks([0 0.15]);
ylabel('Normalized frequency - Shuffle');
xlabel('Correlation values');
set(gca,'fontsize',16);

%ABSOLUTE VALUE OT PREFERRED PHASES
min_r=min(rho_abs);
max_r=max(rho_abs);
min_r_sh=min(rho_sh_abs(:));
max_r_sh=max(rho_sh_abs(:));

% Range of Pearson
for w=1:421
    min_r_sh_perw(w)=min(rho_sh(w,:));
    max_r_sh_perw(w)=max(rho_sh(w,:));
end

% Significance testing Pearson
rho_abs_reshape=abs(rho_abs);
rho_sh_abs_reshape=abs(rho_sh_abs);
for i=1:421
    cutoff(i)=prctile(rho_sh_abs_reshape(i,:),95);
end
aux=find(abs(rho_abs)>cutoff);

pout=myBinomTest(length(aux),421,0.05,'one');






%% Analysis of scatter plot of delta PARTICIPATION INDEX VS pairwise anatomical distance

clear all
close all
rec_data_path='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
data=load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\locking_all_sessions.mat');
pixel_size_new=1.18185;
pixel_size_old=1.78211;
countsc=0;
individual_plots=0;
correlation_sh_pooled=[];
N_sh=100;

for w=1:length(data.waves)
    row_w=data.waves(w);
    disp(w)
    countsc=countsc+1;
    mouse=['L',num2str(data.big_table(row_w,1)),'M',num2str(data.big_table(row_w,2))];
    day=data.big_table(row_w,3);
    s=data.big_table(row_w,4);
    munit=data.big_table(row_w,5);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
%     file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_spk=[data.dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_anat=[data.dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

    load(file_name_anat,'-mat'); %Anatomical information
    load(file_name_spk,'-mat'); %Spike times
    spikes=full(spikes_d_s);
    [N,T]=size(spikes);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mean phases and locked cells
    PI=data.PR_all_sessions{1,w}(:,1);
    not_locked=data.not_locked_all_sessions{w};
    locked=data.locked_all_sessions{w};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distance in the tissue
    if w>10; pixel_size=pixel_size_old;
    else; pixel_size=pixel_size_new; end

    for i=1:N
        r_i(i,:)=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
        dist_origin(i)=norm(r_i(i,:));
    end
    dist_origin_locked=dist_origin(locked);

    count=0;
    for i=1:N
        for j=i+1:N
            count=count+1;
            delta_tissue(i,j)=norm(r_i(j,:)-r_i(i,:));
            delta_tissue(j,i)= delta_tissue(i,j);
            delta_tissue_vec(count)=delta_tissue(i,j);
        end
    end
    delta_tissue(not_locked,:)=[];
    delta_tissue(:,not_locked)=[];
    dist_origin_locked=dist_origin(locked);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distance in PI
    count=0;
    for i=1:N
        for j=i+1:N
            count=count+1;
            delta_PI(i,j)=abs(PI(i)-PI(j));
            delta_PI(j,i)=delta_PI(i,j);
            delta_PI_vec(count)=delta_PI(i,j);
        end
    end
    delta_PI(not_locked,:)=[];
    delta_PI(:,not_locked)=[];

    for sh=1:N_sh
        count=0;
        clear PI_sh delta_PI_vec_sh;
        PI_sh=PI(randperm(N));
        for i=1:N
            for j=i+1:N
                count=count+1;
                delta_PI_vec_sh(count)=abs(PI_sh(i)-PI_sh(j));
            end
        end
        [rho_sh(w,sh),pval_sh(w,sh)] = corr(delta_tissue_vec',delta_PI_vec_sh');
%         [rho_sh_s(w,sh),pval_sh_s(w,sh)] = corr(delta_tissue_vec',delta_PI_vec_sh','type','Spearman');
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save matrices
    matrix_distance_tissue{w}=delta_tissue_vec;
    matrix_distance_PI{w}=delta_PI_vec;
    [rho(w),pval(w)] = corr(delta_tissue_vec',delta_PI_vec');
    [rho_s(w),pval_s(w)] = corr(delta_tissue_vec',delta_PI_vec','type','Spearman');
   

%     correlation_sh_pooled=[correlation_sh_pooled;[w*ones(N_sh,1),rho_sh',pval_sh']];
    %     diff_phase_reshape=delta_phase_mean(:);
    %     diff_pos_reshape=delta_tissue(:);
    clear Anat cells_d count dates delta_PI_vec delta_PI delta_tissue delta_tissue_vec dist_origin dist_origin_locked edges_x locked mean_p ...
        r_i spikes spikes_d_s i j not_locked PI PI_sh   
end


min_r=min(rho);
max_r=max(rho);
min_r_sh=min(rho_sh(:));
max_r_sh=max(rho_sh(:));

% Range of Pearson
for w=1:15
    min_r_sh_perw(w)=min(rho_sh(w,:));
    max_r_sh_perw(w)=max(rho_sh(w,:));
end

% Significance testing Pearson
rho_abs=abs(rho);
rho_sh_abs=abs(rho_sh);
for i=1:15
    cutoff(i)=prctile(rho_sh_abs(i,:),95);
end
aux=find(rho_abs>cutoff);

pout=myBinomTest(length(aux),15,0.05,'one');



% Significance testing
rho_abs=abs(rho);
rho_sh_abs=abs(rho_sh);
for i=1:15
    cutoff(i)=prctile(rho_sh_abs(i,:),95);
end
aux=find(rho_abs>cutoff);


rho_abs=abs(rho_s);
rho_sh_abs=abs(rho_sh_s);
for i=1:15
    cutoff_s(i)=prctile(rho_sh_abs(i,:),95);
end
aux=find(rho_abs>cutoff_s);


%Figures Participation Index
edges_for_boxplot=0:50:1200;

% Individual plots
for w=1:15   
    vect_d=matrix_distance_tissue{1,w};
    vect_p=matrix_distance_PI{1,w};
    vect_d_dis=discretize(vect_d,edges_for_boxplot);

    count_bp=zeros(1,length(edges_for_boxplot));
    mat=nan(length(vect_p),length(edges_for_boxplot));
    for i=1:length(vect_d_dis)
        col=vect_d_dis(i);
        count_bp(col)=count_bp(col)+1;
        mat(count_bp(col),col)=vect_p(i);
    end
    num_bins=find(count_bp==0,1);

    %Scatter plot
    figure
    scatterhist(vect_d,vect_p, 'Marker','.','Color','k');
%     scatter(vect_d,vect_p,'k.');
    ylabel({'Difference in PI'});   %Delta preferred phase
    xlabel('Pairwise anatomical distance (um)');
    set(gca,'fontsize',20,'YColor','k','XColor','k');
    axis([0 inf 0 1]);
    yticks([0 0.5 1]);
    box off
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\PI\DeltaPIVSDist_ScatterPlot_session',num2str(w),'.png'));
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\PI\DeltaPIVSDist_ScatterPlot_session',num2str(w),'.fig'));
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\PI\DeltaPIVSDist_ScatterPlot_session',num2str(w),'.svg'));

    %Box plot
    figure
    boxplot(mat)
    ylabel({'Difference in PI'});   %Delta preferred phase
    xlabel('Pairwise anatomical distance bin');
    set(gca,'fontsize',20,'YColor','k','XColor','k');
    axis([0 num_bins 0 1]);
    xticks(1:num_bins-1);
    yticks([0 0.5 1]);
    box off
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\PI\DeltaPIVSDist_BoxPlot_session',num2str(w),'.png'));
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\PI\DeltaPIVSDist_BoxPlot_session',num2str(w),'.fig'));
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\PI\DeltaPIVSDist_BoxPlot_session',num2str(w),'.svg'));
    
    clear mat_d mat_p mat_d2 mat_p2 vect_d vect_p vect_d_dis count_bp mat
    close all
end



%Heatmaps
edges_phase=0:0.01:1;
edges_dist=0:20:1000;
figure
for w=1:15
    vect_d=matrix_distance_tissue{1,w};
    vect_p=matrix_distance_PI{1,w};
    %     histogram2(X,'CdataMode','auto','Edges',{edges_dist edges_phase},'normalization','probability')
    
    subplot(4,4,w)
    histogram2(vect_d',vect_p',edges_dist,edges_phase,'normalization','probability','DisplayStyle','tile','ShowEmptyBins','on');
    ylabel({'Difference in PI'});   %Delta preferred phase
    xlabel('Anatomical distance (\mum)');
    colorbar
    view(2)
    axis([min(edges_dist) max(edges_dist) min(edges_phase) max(edges_phase)]);
%     xticks([0 500 1000]);
    axis square
    colormap pink
    box off
    set(gca,'fontsize',14,'YColor','k','XColor','k');
    axis([0 1000 0 1]);
    yticks([0 0.5 1]);    
    caxis([0 0.001])
    clear mat_d mat_p mat_d2 mat_p2 vect_d vect_p vect_d_dis count_bp mat
%     close all
end

saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPIVSDist_HeatMap_colorbar_Allsessions.png'));
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPIVSDist_HeatMap_colorbar_Allsessions.fig'));
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPIVSDist_HeatMap_colorbar_Allsessions.svg'));
    
close all


% % 
% % %subplots - scatter plots 
% % figure
% % for w=1:15   
% %     vect_d=matrix_distance_tissue{1,w};
% %     vect_p=matrix_distance_PI{1,w};
% %     vect_d_dis=discretize(vect_d,edges_for_boxplot);
% % 
% %     count_bp=zeros(1,length(edges_for_boxplot));
% %     mat=nan(length(vect_p),length(edges_for_boxplot));
% %     for i=1:length(vect_d_dis)
% %         col=vect_d_dis(i);
% %         count_bp(col)=count_bp(col)+1;
% %         mat(count_bp(col),col)=vect_p(i);
% %     end
% %     num_bins=find(count_bp==0,1);
% % 
% %     %Scatter plot
% %     subplot(4,4,w)
% %     scatter(vect_d,vect_p,'k.');
% %     ylabel({'| \Delta PI |'});   %Delta preferred phase
% %     xlabel('Anatomical distance (um)');
% %     set(gca,'fontsize',16,'YColor','k','XColor','k');
% %     title(['Session ',num2str(w)]);
% %     axis([0 inf -inf 1]);
% %     yticks([0 0.5 1]);
% %     box off
% %    
% %     [rho(w),pval(w)] = corr(vect_d',(vect_p)');
% % 
% %     clear mat_d mat_p mat_d2 mat_p2 vect_d vect_p vect_d_dis count_bp mat
% % end
% % % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\PI\DeltaPIVSDist_ScatterPlot_Allsessions.png'));
% % % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\PI\DeltaPIVSDist_ScatterPlot_Allsessions.fig'));
% % % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\PI\DeltaPIVSDist_ScatterPlot_Allsessions.svg'));
% % % % 
% % 
% % % subplots - box plots 
% % fig=figure;
% % hold on
% % for w=1:15   
% %     vect_d=matrix_distance_tissue{1,w};
% %     vect_p=matrix_distance_PI{1,w};
% %     vect_d_dis=discretize(vect_d,edges_for_boxplot);
% % 
% %     count_bp=zeros(1,length(edges_for_boxplot));
% %     mat=nan(length(vect_p),length(edges_for_boxplot));
% %     for i=1:length(vect_d_dis)
% %         col=vect_d_dis(i);
% %         count_bp(col)=count_bp(col)+1;
% %         mat(count_bp(col),col)=vect_p(i);
% %     end
% %     num_bins=find(count_bp==0,1);
% % 
% %     %Box plot
% % %     subplot(4,4,w)
% % %     boxplot(mat)
% % % %     ylabel({'| \Delta PI |'});   %Delta preferred phase
% % % %     xlabel('Anatomical distance bin');
% % %     set(gca,'fontsize',12,'YColor','k','XColor','k');
% % %     axis([0 num_bins 0 1]);
% % %     xticks(1:num_bins-1);
% % %     yticks([0 0.5 1]);
% % %     title(['Session ',num2str(w)]);
% %     
% %     [p(w),tbl,stats] = anova1(mat);
% %     clear mat_d mat_p mat_d2 mat_p2 vect_d vect_p vect_d_dis count_bp mat
% % end
% % han=axes(fig,'visible','off');
% % han.Title.Visible='on';
% % han.XLabel.Visible='on';
% % han.YLabel.Visible='on';
% % ylabel(han,{'| \Delta PI |'});
% % xlabel(han,'Anatomical distance (bin #)');
% % set(gca,'fontsize',16,'YColor','k','XColor','k');
% % 
% % 
% % 
% % 


% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_BoxPlot_Allsessions.png'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_BoxPlot_Allsessions.fig'));
% saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\DeltaPhase VS DeltaDist\Mean phase\DeltaPPVSDist_BoxPlot_Allsessions.svg'));
%     
% close all
