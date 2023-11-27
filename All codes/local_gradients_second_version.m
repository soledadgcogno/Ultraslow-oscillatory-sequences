%% Analysis for each cycle of the oscillation separately - Phase
clear all
close all

dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
plot_all=1;
[big_table, waves] = get_big_table();

%Params
pixel_size_new=1.18185;
pixel_size_old=1.78211;
N_sh=200;
count_s=0;
% data=load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\locking_all_sessions.mat');

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
        dist_origin(i)=norm(r_i(i,:));
    end
    % Shuffled anatomical position of each cell
%     N=size(spikes_d,1);
%     for sh=1:N_sh
%         sh_ord=randperm(N);
%         r_i_sh(:,:,sh)=r_i(sh_ord,:);
%     end

    delta_tissue=nan(N,N);
    for i=1:N
        for j=i+1:N
            delta_tissue(i,j)=norm(r_i(j,:)-r_i(i,:));
            delta_tissue(j,i)= delta_tissue(i,j);
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

    sequences_per_session(w)=size(table_u,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Here I calculate the phase per sequence
    for s=1:size(table_u,1)
        clear spikes_seq; spikes_seq=spikes(:,table_u(s,1):table_u(s,2));
        clear phase_seq; phase_seq=phase_f(table_u(s,1):table_u(s,2));

        for i=1:N
            clear p; p=phase_seq(find(spikes_seq(i,:)));
            if length(p)>5; mean_p_seq(i,s)=circ_mean(p); else mean_p_seq(i,s)=nan; end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of gradient per sequence
    for s=1:size(table_u,1)
        count_s=count_s+1;
        mean_p=mean_p_seq(:,s);
        delta_phase_mean=nan(N,N);
        for i=1:N
            for j=i+1:N
                delta_phase_mean(i,j)=angdiff(mean_p(i),mean_p(j));
                delta_phase_mean(j,i)=delta_phase_mean(i,j);
            end
        end

        delta_phase_sh=nan(N,N,N_sh);
        for sh=1:N_sh
            mean_p_sh=mean_p(randperm(N));
            for i=1:N
                for j=i+1:N
                    delta_phase_sh(i,j,sh)=angdiff(mean_p_sh(i),mean_p_sh(j));
                    delta_phase_sh(j,i,sh)=delta_phase_sh(i,j,sh);
                end
            end
        end

        delta_phase_50=[];delta_phase_100=[];delta_phase_200=[];
        delta_phase_sh_50=[];delta_phase_sh_100=[];delta_phase_sh_200=[];
        clear template; template=1:N;
        tic
        for i=1:N
            clear neighbor_cells_50; neighbor_cells_50=find(delta_tissue(i,:)<50 & delta_tissue(i,:)>0); delta_phase_50=[delta_phase_50,delta_phase_mean(i,neighbor_cells_50)];
            clear neighbor_cells_100; neighbor_cells_100=find(delta_tissue(i,:)<100); delta_phase_100=[delta_phase_100,delta_phase_mean(i,neighbor_cells_100)];
            clear neighbor_cells_200; neighbor_cells_200=find(delta_tissue(i,:)<200); delta_phase_200=[delta_phase_200,delta_phase_mean(i,neighbor_cells_200)];

            clear phases_sh; phases_sh=reshape(delta_phase_sh(i,neighbor_cells_50,:),[length(neighbor_cells_50),N_sh,1]);
%             clear temp50; temp50=repmat(ismember(template,neighbor_cells_50),1,N_sh);
            delta_phase_sh_50=[delta_phase_sh_50;phases_sh];

            clear phases_sh; phases_sh=reshape(delta_phase_sh(i,neighbor_cells_100,:),[length(neighbor_cells_100),N_sh,1]);
            delta_phase_sh_100=[delta_phase_sh_100;phases_sh];

            clear phases_sh; phases_sh=reshape(delta_phase_sh(i,neighbor_cells_200,:),[length(neighbor_cells_200),N_sh,1]);
            delta_phase_sh_200=[delta_phase_sh_200;phases_sh];
        end
        toc

        %Medians and means
        mean_data50(count_s)= mean(abs(delta_phase_50(~isnan(delta_phase_50))));
        median_data50(count_s)= median(abs(delta_phase_50(~isnan(delta_phase_50))));
        var_data50(count_s)= var(abs(delta_phase_50(~isnan(delta_phase_50))));

        mean_data100(count_s)= mean(abs(delta_phase_100(~isnan(delta_phase_100))));
        median_data100(count_s)= median(abs(delta_phase_100(~isnan(delta_phase_100))));
        var_data100(count_s)= var(abs(delta_phase_50(~isnan(delta_phase_50))));

        mean_data200(count_s)= mean(abs(delta_phase_200(~isnan(delta_phase_200))));
        median_data200(count_s)= median(abs(delta_phase_200(~isnan(delta_phase_200))));
        var_data200(count_s)= var(abs(delta_phase_50(~isnan(delta_phase_50))));

        tic
        for sh=1:N_sh
            clear aux aux_wonan ; aux=delta_phase_sh_50(:,sh); aux_wonan=aux(~isnan(aux));
            mean_data50_sh(count_s,sh)= mean(abs(aux_wonan));
            median_data50_sh(count_s,sh)= median(abs(aux_wonan));
            var_data50_sh(count_s,sh)= var(abs(aux_wonan));

            clear aux aux_wonan ; aux=delta_phase_sh_100(:,sh); aux_wonan=aux(~isnan(aux));
            mean_data100_sh(count_s,sh)= mean(abs(aux_wonan));
            median_data100_sh(count_s,sh)= median(abs(aux_wonan));
            var_data100_sh(count_s,sh)= var(abs(aux_wonan));

            clear aux aux_wonan ; aux=delta_phase_sh_200(:,sh); aux_wonan=aux(~isnan(aux));
            mean_data200_sh(count_s,sh)= mean(abs(aux_wonan));
            median_data200_sh(count_s,sh)= median(abs(aux_wonan));
            var_data200_sh(count_s,sh)= var(abs(aux_wonan));
        end
        toc

        mean_50_5(count_s)=prctile(mean_data50_sh(count_s,:),5);
        mean_50_95(count_s)=prctile(mean_data50_sh(count_s,:),95);
        mean_100_5(count_s)=prctile(mean_data100_sh(count_s,:),5);
        mean_100_95(count_s)=prctile(mean_data100_sh(count_s,:),95);
        mean_200_5(count_s)=prctile(mean_data200_sh(count_s,:),5);
        mean_200_95(count_s)=prctile(mean_data200_sh(count_s,:),95);

        median_50_5(count_s)=prctile(median_data50_sh(count_s,:),5);
        median_50_95(count_s)=prctile(median_data50_sh(count_s,:),95);
        median_100_5(count_s)=prctile(median_data100_sh(count_s,:),5);
        median_100_95(count_s)=prctile(median_data100_sh(count_s,:),95);
        median_200_5(count_s)=prctile(median_data200_sh(count_s,:),5);
        median_200_95(count_s)=prctile(median_data200_sh(count_s,:),95);

        var_50_5(count_s)=prctile(var_data50_sh(count_s,:),5);
        var_50_95(count_s)=prctile(var_data50_sh(count_s,:),95);
        var_100_5(count_s)=prctile(var_data100_sh(count_s,:),5);
        var_100_95(count_s)=prctile(var_data100_sh(count_s,:),95);
        var_200_5(count_s)=prctile(var_data200_sh(count_s,:),5);
        var_200_95(count_s)=prctile(var_data200_sh(count_s,:),95);

        if (w==8 && s==19)
            figure
            subplot(3,2,1)
            histogram(mean_data50_sh(count_s,:),[1.4:0.01:1.7]);
            hold on
            xline(mean_data50(count_s),'--r','LineWidth',2);
            title('Mean - 50 \mum');
            subplot(3,2,2)
            histogram(median_data50_sh(count_s,:),[1.4:0.01:1.7]);
            hold on
            xline(median_data50(count_s),'--r','LineWidth',2);
            title('Median - 50 \mum');
%             subplot(3,3,3)
%             histogram(var_data50_sh(count_s,:),[0.7:0.01:1]);
%             hold on
%             xline(var_data50(count_s),'--r','LineWidth',2);
%             title('Var - 50 \mum');
            
            subplot(3,2,3)
            histogram(mean_data100_sh(count_s,:),[1.4:0.01:1.7]);
            hold on
            xline(mean_data100(count_s),'--r','LineWidth',2);
            title('Mean - 100 \mum');
            subplot(3,2,4)
            histogram(median_data100_sh(count_s,:),[1.4:0.01:1.7]);
            hold on
            xline(median_data100(count_s),'--r','LineWidth',2);
            title('Median - 100 \mum');
%             subplot(3,3,6)
%             histogram(var_data100_sh(count_s,:),[0.7:0.01:1]);
%             hold on
%             xline(var_data100(count_s),'--r','LineWidth',2);
%             title('Var - 100 \mum');

            subplot(3,2,5)
            histogram(mean_data200_sh(count_s,:),[1.4:0.01:1.7]);
            hold on
            xline(mean_data200(count_s),'--r','LineWidth',2);
            title('Mean - 200 \mum');
            subplot(3,2,6)
            histogram(median_data200_sh(count_s,:),[1.4:0.01:1.7]);
            hold on
            xline(median_data200(count_s),'--r','LineWidth',2);
            title('Median - 200 \mum');
%             subplot(3,3,9)
%             histogram(var_data200_sh(count_s,:),[0.7:0.01:1]);
%             hold on
%             xline(var_data200(count_s),'--r','LineWidth',2);
%             title('Var - 200 \mum');
        end

        %Tests
%         [h50(count_s),p50(count_s)] = kstest2(delta_phase_50(~isnan(delta_phase_50)),delta_phase_sh_50(~isnan(delta_phase_sh_50)));
%         [h100(count_s),p100(count_s)] = kstest2(delta_phase_100,delta_phase_sh_100);
%         [h200(count_s),p200(count_s)] = kstest2(delta_phase_200(~isnan(delta_phase_200)),delta_phase_sh_200(~isnan(delta_phase_sh_200)));

        delta_phase_sh_50_vec=abs(delta_phase_sh_50(:));
        delta_phase_sh_100_vec=abs(delta_phase_sh_100(:));
        delta_phase_sh_200_vec=abs(delta_phase_sh_200(:));

        [p50_RS(count_s),h50_RS(count_s),stats{count_s}] = ranksum(abs(delta_phase_50(~isnan(delta_phase_50))),delta_phase_sh_50_vec(~isnan(delta_phase_sh_50_vec)),'tail','left');
        [p100_RS(count_s),h100_RS(count_s),stats{count_s}] = ranksum(abs(delta_phase_100(~isnan(delta_phase_100))),delta_phase_sh_100_vec(~isnan(delta_phase_sh_100_vec)),'tail','left');
        [p200_RS(count_s),h200_RS(count_s),stats{count_s}] = ranksum(abs(delta_phase_200(~isnan(delta_phase_200))),delta_phase_sh_200_vec(~isnan(delta_phase_sh_200_vec)),'tail','left');

        clear mean_p mean_p_sh delta_phase_mean delta_phase_sh delta_phase_sh_100_vec delta_phase_sh_200_vec delta_phase_sh_50_vec
    end
    clear phase_f phase_r phase_seq radius_f spikes_seq spikes_r spikes_d_s spikes_d spikes scoret template ...
        phase mean_p_seq FRp r_i r_i_sh table_u sh_ord cells_d delta_tissue delta_phase_mean delta_phase_sh ...
        dist_origin coefft Anat phases_sh

end

%% Quantifications

res_mean50=length(find(mean_data50<mean_50_5));
res_median50=length(find(median_data50<median_50_5));
res_var50=length(find(var_data50<var_50_5));

res_mean100=length(find(mean_data100<mean_100_5));
res_median100=length(find(median_data100<median_100_5));
res_var100=length(find(var_data100<var_100_5));

res_mean200=length(find(mean_data200<mean_200_5));
res_median200=length(find(median_data200<median_200_5));
res_var200=length(find(var_data200<var_200_5));

pout=myBinomTest(15,421,0.05,'one');


l1=find(p50_RS<0.05);
l2=find(p100_RS<0.05);
l3=find(p200_RS<0.05);

myBinomTest(1,421,0.05,'one')

%% Figures for all sessions

%1. Figures with percentile
ed=[0:0.01:pi];

[mean_50_shuffle_dist,edges]=histcounts(mean_50_5,ed,'Normalization','probability');
binwidth=edges(2)-edges(1);
[mean_50_dist,edges]=histcounts(mean_data50,ed,'Normalization','probability');
[median_50_shuffle_dist,edges]=histcounts(median_50_5,ed,'Normalization','probability');
[median_50_dist,edges]=histcounts(median_data50,ed,'Normalization','probability');

[mean_100_shuffle_dist,edges]=histcounts(mean_100_5,ed,'Normalization','probability');
binwidth=edges(2)-edges(1);
[mean_100_dist,edges]=histcounts(mean_data100,ed,'Normalization','probability');
[median_100_shuffle_dist,edges]=histcounts(median_100_5,ed,'Normalization','probability');
[median_100_dist,edges]=histcounts(median_data100,ed,'Normalization','probability');

[mean_200_shuffle_dist,edges]=histcounts(mean_200_5,ed,'Normalization','probability');
binwidth=edges(2)-edges(1);
[mean_200_dist,edges]=histcounts(mean_data200,ed,'Normalization','probability');
[median_200_shuffle_dist,edges]=histcounts(median_200_5,ed,'Normalization','probability');
[median_200_dist,edges]=histcounts(median_data200,ed,'Normalization','probability');


fig=figure;
subplot(3,2,1)
plot(edges(1:end-1)+binwidth/2,cumsum(mean_50_dist),'linewidth',2.5);
axis([0.4 1.7 0 inf]);
title('Mean - 50 \mum');
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(mean_50_shuffle_dist),'linewidth',2.5);
box off; axis square;
legend({'Data','Shuffle'});
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
subplot(3,2,2)
plot(edges(1:end-1)+binwidth/2,cumsum(median_50_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(median_50_shuffle_dist),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Median - 50 \mum');
subplot(3,2,3)
plot(edges(1:end-1)+binwidth/2,cumsum(mean_100_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(mean_100_shuffle_dist),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Mean - 100 \mum');
subplot(3,2,4)
plot(edges(1:end-1)+binwidth/2,cumsum(median_100_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(median_100_shuffle_dist),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Median - 100 \mum');
subplot(3,2,5)
plot(edges(1:end-1)+binwidth/2,cumsum(mean_200_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(mean_200_shuffle_dist),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Mean - 200 \mum');
subplot(3,2,6)
plot(edges(1:end-1)+binwidth/2,cumsum(median_200_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(median_200_shuffle_dist),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Median - 200 \mum');
han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Cumulative normalized frequency');
xlabel(han,'Difference in preferred phase (rad)');
set(gca,'fontsize',16);


%1. Figures with mean or median


for seq=1:421
    mean_mean_sh_50(seq)=mean(mean_data50_sh(seq,:));
    median_mean_sh_50(seq)=median(mean_data50_sh(seq,:));
    mean_median_sh_50(seq)=mean(median_data50_sh(seq,:));
    median_median_sh_50(seq)=median(median_data50_sh(seq,:));

    mean_mean_sh_100(seq)=mean(mean_data100_sh(seq,:));
    median_mean_sh_100(seq)=median(mean_data100_sh(seq,:));
    mean_median_sh_100(seq)=mean(median_data100_sh(seq,:));
    median_median_sh_100(seq)=median(median_data100_sh(seq,:));

    mean_mean_sh_200(seq)=mean(mean_data200_sh(seq,:));
    median_mean_sh_200(seq)=median(mean_data200_sh(seq,:));
    mean_median_sh_200(seq)=mean(median_data200_sh(seq,:));
    median_median_sh_200(seq)=median(median_data200_sh(seq,:));
end

[mean_mean_50_shuffle,edges]=histcounts(mean_mean_sh_50,ed,'Normalization','probability');
binwidth=edges(2)-edges(1);
[median_mean_50_shuffle,~]=histcounts(mean_mean_sh_50,ed,'Normalization','probability');

[mean_median_50_shuffle,edges]=histcounts(mean_median_sh_50,ed,'Normalization','probability');
binwidth=edges(2)-edges(1);
[median_median_50_shuffle,~]=histcounts(mean_median_sh_50,ed,'Normalization','probability');

[mean_mean_100_shuffle,edges]=histcounts(mean_mean_sh_100,ed,'Normalization','probability');
binwidth=edges(2)-edges(1);
[median_mean_100_shuffle,~]=histcounts(mean_mean_sh_100,ed,'Normalization','probability');

[mean_median_100_shuffle,edges]=histcounts(mean_median_sh_100,ed,'Normalization','probability');
binwidth=edges(2)-edges(1);
[median_median_100_shuffle,~]=histcounts(mean_median_sh_100,ed,'Normalization','probability');

[mean_mean_200_shuffle,edges]=histcounts(mean_mean_sh_200,ed,'Normalization','probability');
binwidth=edges(2)-edges(1);
[median_mean_200_shuffle,~]=histcounts(mean_mean_sh_200,ed,'Normalization','probability');

[mean_median_200_shuffle,edges]=histcounts(mean_median_sh_200,ed,'Normalization','probability');
binwidth=edges(2)-edges(1);
[median_median_200_shuffle,~]=histcounts(mean_median_sh_200,ed,'Normalization','probability');

%Plotting Mean of shuffled distribution
fig=figure;
subplot(3,2,1)
plot(edges(1:end-1)+binwidth/2,cumsum(mean_50_dist),'linewidth',2.5);
axis([0.4 1.7 0 inf]);
title('Mean - 50 \mum');
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(mean_mean_50_shuffle),'linewidth',2.5);
box off; axis square;
legend({'Data','Shuffle'});
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
subplot(3,2,2)
plot(edges(1:end-1)+binwidth/2,cumsum(median_50_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(mean_median_50_shuffle),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Median - 50 \mum');
subplot(3,2,3)
plot(edges(1:end-1)+binwidth/2,cumsum(mean_100_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(mean_mean_100_shuffle),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Mean - 100 \mum');
subplot(3,2,4)
plot(edges(1:end-1)+binwidth/2,cumsum(median_100_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(mean_median_100_shuffle),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Median - 100 \mum');
subplot(3,2,5)
plot(edges(1:end-1)+binwidth/2,cumsum(mean_200_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(mean_mean_200_shuffle),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Mean - 200 \mum');
subplot(3,2,6)
plot(edges(1:end-1)+binwidth/2,cumsum(median_200_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(mean_median_200_shuffle),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Median - 200 \mum');
han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Cumulative normalized frequency');
xlabel(han,'Difference in preferred phase (rad)');
set(gca,'fontsize',16);



%Plotting Median of shuffled distribution
fig=figure;
subplot(3,2,1)
plot(edges(1:end-1)+binwidth/2,cumsum(mean_50_dist),'linewidth',2.5);
axis([0.4 1.7 0 inf]);
title('Mean - 50 \mum');
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(median_mean_50_shuffle),'linewidth',2.5);
box off; axis square;
legend({'Data','Shuffle'});
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
subplot(3,2,2)
plot(edges(1:end-1)+binwidth/2,cumsum(median_50_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(median_median_50_shuffle),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Median - 50 \mum');
subplot(3,2,3)
plot(edges(1:end-1)+binwidth/2,cumsum(mean_100_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(median_mean_100_shuffle),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Mean - 100 \mum');
subplot(3,2,4)
plot(edges(1:end-1)+binwidth/2,cumsum(median_100_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(median_median_100_shuffle),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Median - 100 \mum');
subplot(3,2,5)
plot(edges(1:end-1)+binwidth/2,cumsum(mean_200_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(median_mean_200_shuffle),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Mean - 200 \mum');
subplot(3,2,6)
plot(edges(1:end-1)+binwidth/2,cumsum(median_200_dist),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(median_median_200_shuffle),'linewidth',2.5);
box off; axis square;
set(gca,'fontsize',16);
axis([0.4 1.7 0 inf]);
title('Median - 200 \mum');
han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Cumulative normalized frequency');
xlabel(han,'Difference in preferred phase (rad)');
set(gca,'fontsize',16);


%% Figure for one session (I run everything again for that one session and then I make the figure)
clear all
close all
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
plot_all=1;
[big_table, waves] = get_big_table();
%Params
pixel_size_new=1.18185;
pixel_size_old=1.78211;
N_sh=200;
count_s=0;
% data=load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\locking_all_sessions.mat');
for w=8%:length(waves)
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
        dist_origin(i)=norm(r_i(i,:));
    end
    % Shuffled anatomical position of each cell
%     N=size(spikes_d,1);
%     for sh=1:N_sh
%         sh_ord=randperm(N);
%         r_i_sh(:,:,sh)=r_i(sh_ord,:);
%     end

    delta_tissue=nan(N,N);
    for i=1:N
        for j=i+1:N
            delta_tissue(i,j)=norm(r_i(j,:)-r_i(i,:));
            delta_tissue(j,i)= delta_tissue(i,j);
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

    sequences_per_session(w)=size(table_u,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Here I calculate the phase per sequence
    for s=1:size(table_u,1)
        clear spikes_seq; spikes_seq=spikes(:,table_u(s,1):table_u(s,2));
        clear phase_seq; phase_seq=phase_f(table_u(s,1):table_u(s,2));

        for i=1:N
            clear p; p=phase_seq(find(spikes_seq(i,:)));
            if length(p)>5; mean_p_seq(i,s)=circ_mean(p); else mean_p_seq(i,s)=nan; end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of gradient per sequence
    for s=19%1:size(table_u,1)
        count_s=count_s+1;
        mean_p=mean_p_seq(:,s);
        delta_phase_mean=nan(N,N);
        for i=1:N
            for j=i+1:N
                delta_phase_mean(i,j)=angdiff(mean_p(i),mean_p(j));
                delta_phase_mean(j,i)=delta_phase_mean(i,j);
            end
        end

        delta_phase_sh=nan(N,N,N_sh);
        for sh=1:N_sh
            mean_p_sh=mean_p(randperm(N));
            for i=1:N
                for j=i+1:N
                    delta_phase_sh(i,j,sh)=angdiff(mean_p_sh(i),mean_p_sh(j));
                    delta_phase_sh(j,i,sh)=delta_phase_sh(i,j,sh);
                end
            end
        end

        delta_phase_50=[];delta_phase_100=[];delta_phase_200=[];
        delta_phase_sh_50=[];delta_phase_sh_100=[];delta_phase_sh_200=[];
        clear template; template=1:N;
        tic
        for i=1:N
            clear neighbor_cells_50; neighbor_cells_50=find(delta_tissue(i,:)<50 & delta_tissue(i,:)>0); delta_phase_50=[delta_phase_50,delta_phase_mean(i,neighbor_cells_50)];
            clear neighbor_cells_100; neighbor_cells_100=find(delta_tissue(i,:)<100); delta_phase_100=[delta_phase_100,delta_phase_mean(i,neighbor_cells_100)];
            clear neighbor_cells_200; neighbor_cells_200=find(delta_tissue(i,:)<200); delta_phase_200=[delta_phase_200,delta_phase_mean(i,neighbor_cells_200)];

            clear phases_sh; phases_sh=reshape(delta_phase_sh(i,neighbor_cells_50,:),[length(neighbor_cells_50),N_sh,1]);
%             clear temp50; temp50=repmat(ismember(template,neighbor_cells_50),1,N_sh);
            delta_phase_sh_50=[delta_phase_sh_50;phases_sh];

            clear phases_sh; phases_sh=reshape(delta_phase_sh(i,neighbor_cells_100,:),[length(neighbor_cells_100),N_sh,1]);
            delta_phase_sh_100=[delta_phase_sh_100;phases_sh];

            clear phases_sh; phases_sh=reshape(delta_phase_sh(i,neighbor_cells_200,:),[length(neighbor_cells_200),N_sh,1]);
            delta_phase_sh_200=[delta_phase_sh_200;phases_sh];
        end
        toc

        %Medians and means
        mean_data50(count_s)= mean(abs(delta_phase_50(~isnan(delta_phase_50))));
        median_data50(count_s)= median(abs(delta_phase_50(~isnan(delta_phase_50))));
        var_data50(count_s)= var(abs(delta_phase_50(~isnan(delta_phase_50))));

        mean_data100(count_s)= mean(abs(delta_phase_100(~isnan(delta_phase_100))));
        median_data100(count_s)= median(abs(delta_phase_100(~isnan(delta_phase_100))));
        var_data100(count_s)= var(abs(delta_phase_50(~isnan(delta_phase_50))));

        mean_data200(count_s)= mean(abs(delta_phase_200(~isnan(delta_phase_200))));
        median_data200(count_s)= median(abs(delta_phase_200(~isnan(delta_phase_200))));
        var_data200(count_s)= var(abs(delta_phase_50(~isnan(delta_phase_50))));

        tic
        for sh=1:N_sh
            clear aux aux_wonan ; aux=delta_phase_sh_50(:,sh); aux_wonan=aux(~isnan(aux));
            mean_data50_sh(count_s,sh)= mean(abs(aux_wonan));
            median_data50_sh(count_s,sh)= median(abs(aux_wonan));
            var_data50_sh(count_s,sh)= var(abs(aux_wonan));

            clear aux aux_wonan ; aux=delta_phase_sh_100(:,sh); aux_wonan=aux(~isnan(aux));
            mean_data100_sh(count_s,sh)= mean(abs(aux_wonan));
            median_data100_sh(count_s,sh)= median(abs(aux_wonan));
            var_data100_sh(count_s,sh)= var(abs(aux_wonan));

            clear aux aux_wonan ; aux=delta_phase_sh_200(:,sh); aux_wonan=aux(~isnan(aux));
            mean_data200_sh(count_s,sh)= mean(abs(aux_wonan));
            median_data200_sh(count_s,sh)= median(abs(aux_wonan));
            var_data200_sh(count_s,sh)= var(abs(aux_wonan));
        end
        toc

        mean_50_5(count_s)=prctile(mean_data50_sh(count_s,:),5);
        mean_50_95(count_s)=prctile(mean_data50_sh(count_s,:),95);
        mean_100_5(count_s)=prctile(mean_data100_sh(count_s,:),5);
        mean_100_95(count_s)=prctile(mean_data100_sh(count_s,:),95);
        mean_200_5(count_s)=prctile(mean_data200_sh(count_s,:),5);
        mean_200_95(count_s)=prctile(mean_data200_sh(count_s,:),95);

        median_50_5(count_s)=prctile(median_data50_sh(count_s,:),5);
        median_50_95(count_s)=prctile(median_data50_sh(count_s,:),95);
        median_100_5(count_s)=prctile(median_data100_sh(count_s,:),5);
        median_100_95(count_s)=prctile(median_data100_sh(count_s,:),95);
        median_200_5(count_s)=prctile(median_data200_sh(count_s,:),5);
        median_200_95(count_s)=prctile(median_data200_sh(count_s,:),95);

        var_50_5(count_s)=prctile(var_data50_sh(count_s,:),5);
        var_50_95(count_s)=prctile(var_data50_sh(count_s,:),95);
        var_100_5(count_s)=prctile(var_data100_sh(count_s,:),5);
        var_100_95(count_s)=prctile(var_data100_sh(count_s,:),95);
        var_200_5(count_s)=prctile(var_data200_sh(count_s,:),5);
        var_200_95(count_s)=prctile(var_data200_sh(count_s,:),95);

        delta_phase_sh_50_vec=abs(delta_phase_sh_50(:));
        delta_phase_sh_100_vec=abs(delta_phase_sh_100(:));
        delta_phase_sh_200_vec=abs(delta_phase_sh_200(:));

        [p50_RS(count_s),h50_RS(count_s),stats{count_s}] = ranksum(abs(delta_phase_50(~isnan(delta_phase_50))),delta_phase_sh_50_vec(~isnan(delta_phase_sh_50_vec)),'tail','left');
        [p100_RS(count_s),h100_RS(count_s),stats{count_s}] = ranksum(abs(delta_phase_100(~isnan(delta_phase_100))),delta_phase_sh_100_vec(~isnan(delta_phase_sh_100_vec)),'tail','left');
        [p200_RS(count_s),h200_RS(count_s),stats{count_s}] = ranksum(abs(delta_phase_200(~isnan(delta_phase_200))),delta_phase_sh_200_vec(~isnan(delta_phase_sh_200_vec)),'tail','left');

%         clear mean_p mean_p_sh delta_phase_mean delta_phase_sh delta_phase_sh_100_vec delta_phase_sh_200_vec delta_phase_sh_50_vec
    end
%     clear phase_f phase_r phase_seq radius_f spikes_seq spikes_r spikes_d_s spikes_d spikes scoret template ...
%         phase mean_p_seq FRp r_i r_i_sh table_u sh_ord cells_d delta_tissue delta_phase_mean delta_phase_sh ...
%         dist_origin coefft Anat phases_sh
end

ed=[0:0.01:pi];

ed=[1.4:0.01:1.7];
fig=figure;
subplot(3,2,1)
histogram(mean_data50_sh(count_s,:),ed);
hold on
xline( mean_50_5(count_s),'b-','linewidth',5);
xline( mean_data50(count_s),'r-','linewidth',5);
axis([1.45 1.7 0 inf])
box off
axis square; set(gca,'fontsize',16);
title('Mean - 50 \mum');
subplot(3,2,2)
histogram(median_data50_sh(count_s,:),ed);
hold on
xline( median_50_5(count_s),'b-','linewidth',5);
xline( median_data50(count_s),'r-','linewidth',5);
axis([1.45 1.7 0 inf])
box off
axis square; set(gca,'fontsize',16);
title('Median - 50 \mum');
subplot(3,2,3)
histogram(mean_data100_sh(count_s,:),ed);
hold on
xline( mean_100_5(count_s),'b-','linewidth',5);
xline( mean_data100(count_s),'r-','linewidth',5);
axis([1.45 1.7 0 inf])
box off
axis square; set(gca,'fontsize',16);
title('Mean - 100 \mum');
subplot(3,2,4)
histogram(median_data100_sh(count_s,:),ed);
hold on
xline( median_100_5(count_s),'b-','linewidth',5);
xline( median_data100(count_s),'r-','linewidth',5);
axis([1.45 1.7 0 inf])
box off
axis square; set(gca,'fontsize',16);
title('Median - 100 \mum');
subplot(3,2,5)
histogram(mean_data200_sh(count_s,:),ed);
hold on
xline( mean_200_5(count_s),'b-','linewidth',5);
xline( mean_data200(count_s),'r-','linewidth',5);
axis([1.45 1.7 0 inf])
box off
axis square; set(gca,'fontsize',16);
title('Mean - 200 \mum');
subplot(3,2,6)
histogram(median_data200_sh(count_s,:),ed);
hold on
xline( median_200_5(count_s),'b-','linewidth',5);
xline( median_data200(count_s),'r-','linewidth',5);
axis([1.45 1.7 0 inf])
box off
axis square; set(gca,'fontsize',16);
title('Median - 200 \mum');
han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Cumulative normalized frequency');
xlabel(han,'Difference in preferred phase (rad)');
set(gca,'fontsize',16);

