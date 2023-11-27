%% Analysis for each cycle of the oscillation separately - Phase
%Computes the MVL in spatial bins

clear all
close all

dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
plot_all=1;
[big_table, waves] = get_big_table();

%Params
pixel_size_new=1.18185;
pixel_size_old=1.78211;
N_sh=1000;
count_s=0;
min_cutoff_n_neurons=10; %threshold for neurons in spatial bins
% data=load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\locking_all_sessions.mat');

%Params for spatial bins
spatial_bin_size=200;
bin_size=spatial_bin_size;      %in um
max_length=1000;     %size of FOV in x direction
n_bins=max_length/bin_size;

%Arrays initialization
circp_pos_sh=nan(n_bins,n_bins,N_sh);
meanp_pos_sh=nan(n_bins,n_bins,N_sh);
occupancy_mat=zeros(n_bins,n_bins,500);
occupancy_mat_b=zeros(n_bins,n_bins,500);
meanp_pos=nan(n_bins,n_bins,500);
circp_pos=nan(n_bins,n_bins,500);
% prctile_1=nan(n_bins,n_bins,500);
% prctile_5=nan(n_bins,n_bins,500);
% prctile_10=nan(n_bins,n_bins,500);
prctile_90=nan(n_bins,n_bins,500);
prctile_95=nan(n_bins,n_bins,500);
prctile_99=nan(n_bins,n_bins,500);

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Here I calculate the phase per sequence
    for s=1:size(table_u,1)
        clear spikes_seq; spikes_seq=spikes(:,table_u(s,1):table_u(s,2));
        clear phase_seq; phase_seq=phase_f(table_u(s,1):table_u(s,2));

        for i=1:N
            clear p; p=phase_seq(find(spikes_seq(i,:)));
            if length(p)>5; mean_p_seq(i,s)=circ_mean(p); else mean_p_seq(i,s)=nan; end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binning of cells' position
    [C, idx] = bin_cells_position(r_i,bin_size,max_length); %C is the table with the bins, and the data is in idx
    %     [~, idx_tw] = bin_cells_position(r_i_tw,bin_size,max_length);
    %     for sh=1:N_sh
    %         [~, idx_sh(:,sh)] = bin_cells_position(r_i_sh(:,:,sh),bin_size,max_length);
    %     end
    %Reaccomodate the C matrix for plotting it
    C_s=[];
    for m=1:n_bins
        clear aux; aux=1:n_bins;
        C_s=[C_s;[aux',(n_bins-m+1)*ones(n_bins,1)]];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of local phases per sequence
    for s=1:size(table_u,1)
        count_s=count_s+1;
        clear mean_p; mean_p=mean_p_seq(:,s);
        meanp_pos_sh=nan(n_bins,n_bins,N_sh);
        circp_pos_sh=nan(n_bins,n_bins,N_sh);
%         clear meanp_pos_sh circp_pos_sh;
        for sh=1:N_sh
            clear mean_p_sh; mean_p_sh=mean_p(randperm(N));
            for m=1:size(C,1) % Loop on all spatial bins
                if sum(idx==m)>0
                    clear aux aux2; aux=(mean_p_sh(idx==m)); aux2=aux(~isnan(aux));
                    if length(aux2)>min_cutoff_n_neurons
                        meanp_pos_sh(C(m,1),C(m,2),sh)=circ_mean(aux2);
                        circp_pos_sh(C(m,1),C(m,2),sh)=circ_r(aux2);
                    end
                end
            end
        end

        for m=1:size(C,1) % Loop on all spatial bins
            if sum(idx==m)>0
                occupancy_mat(C(m,1),C(m,2),count_s)=sum(idx==m);
                clear aux aux2; aux=(mean_p(idx==m)); aux2=aux(~isnan(aux));
                if length(aux2)>min_cutoff_n_neurons
                    occupancy_mat_b(C(m,1),C(m,2),count_s)=length(aux2);
                    meanp_pos(C(m,1),C(m,2),count_s)=circ_mean(aux2);
                    circp_pos(C(m,1),C(m,2),count_s)=circ_r(aux2);

                    prctile_90(C(m,1),C(m,2),count_s)=prctile(circp_pos_sh(C(m,1),C(m,2),:),90);
                    prctile_95(C(m,1),C(m,2),count_s)=prctile(circp_pos_sh(C(m,1),C(m,2),:),95);
                    prctile_99(C(m,1),C(m,2),count_s)=prctile(circp_pos_sh(C(m,1),C(m,2),:),99);
                end
            end
        end
    end

    clear phase_f phase_r phase_seq radius_f spikes_seq spikes_r spikes_d_s spikes_d spikes scoret template ...
        phase mean_p_seq FRp r_i r_i_sh table_u sh_ord cells_d delta_tissue delta_phase_mean delta_phase_sh ...
        dist_origin coefft Anat phases_sh C idx mean_p mean_p_sh aux aux2 C_s 
end

% figure
% scatter(circp_pos(:),prctile_90(:),'k','filled');
% ylabel({'90^t^h percentile' ; 'Shuffled data'});
% set(gca,'fontsize',16);
% xlabel('Experimental data');
% % title([num2str(n_bins),'x',num2str(n_bins),' bins']);
% title(['200\mum x 200\mum bins']);
% hl=refline(1,0);



%New quantification
tot_largerthanten=find(occupancy_mat_b(:)>10); %spatial bins with more than 10 cells
circp_pos_b=circp_pos(:);
circp_pos_b=circp_pos_b(tot_largerthanten);
prctile_95_b=prctile_95(:);
prctile_95_b=prctile_95_b(tot_largerthanten);

num=find(circp_pos_b>prctile_95_b);
tot=length(num)*100/length(circp_pos_b);
disp(['Prob = ',num2str(tot)])

% Figure
figure
scatter(circp_pos_b(:),prctile_95_b(:),'k','filled');
ylabel({'95^t^h percentile' ; 'Shuffled data'});
set(gca,'fontsize',16);
xlabel('Experimental data');
% title([num2str(n_bins),'x',num2str(n_bins),' bins']);
title(['200\mum x 200\mum bins']);
hl=refline(1,0);
hl.LineWidth=5;
hl.LineStyle='--';

% %Old quantification
% tot=sum(~isnan(circp_pos(:)));
% aux=find(circp_pos(:)>prctile_95(:));
% pout=myBinomTest(length(aux),tot,0.05,'one');

% figure
% subplot(1,3,1)
% scatter(circp_pos(:),prctile_90(:))
% ylabel('Shuffled data');
% set(gca,'fontsize',16);
% xlabel('Experimental data');
% title('90^t^h percentile');
% hl=refline(1,0);
% hl.Color = 'r';
% subplot(1,3,2)
% scatter(circp_pos(:),prctile_95(:))
% ylabel('Shuffled data');
% title('95^t^h percentile');
% ylabel('Shuffled data');
% xlabel('Data');
% set(gca,'fontsize',16);
% xlabel('Experimental data');
% hl=refline(1,0);
% hl.Color = 'r';
% subplot(1,3,3)
% scatter(circp_pos(:),prctile_99(:))
% ylabel('99th percentile - Shuffle');
% ylabel('Shuffled data');
% xlabel('Experimental data');
% set(gca,'fontsize',16);
% xlabel('Experimental data');
% title('99^t^h percentile');
% hl=refline(1,0);
% hl.Color = 'r';

% %Quantification
% tot=sum(~isnan(circp_pos(:)));
% ratio_99=length(find(circp_pos(:)>prctile_99(:)))/sum(~isnan(circp_pos(:)));
% ratio_95=length(find(circp_pos(:)>prctile_95(:)))/sum(~isnan(circp_pos(:)));
% ratio_90=length(find(circp_pos(:)>prctile_90(:)))/sum(~isnan(circp_pos(:)));

%% Example session and sequences
% Analysis for each cycle of the oscillation separately - Phase
clear all
% close all

dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
plot_all=1;
[big_table, waves] = get_big_table();

%Params
pixel_size_new=1.18185;
pixel_size_old=1.78211;
N_sh=1000;
count_s=0;
% data=load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\locking_all_sessions.mat');

%Params for spatial bins
spatial_bin_size=100;
bin_size=spatial_bin_size;      %in um
max_length=600;     %size of FOV in x direction
n_bins=max_length/bin_size;
min_cutoff_n_neurons=10; %threshold for neurons in spatial bins

%Arrays initialization
circp_pos_sh=nan(n_bins,n_bins,N_sh);
meanp_pos_sh=nan(n_bins,n_bins,N_sh);
occupancy_mat=zeros(n_bins,n_bins,1);
occupancy_mat_b=zeros(n_bins,n_bins,1);
meanp_pos=nan(n_bins,n_bins,1);
circp_pos=nan(n_bins,n_bins,1);
prctile_90=nan(n_bins,n_bins,1);
prctile_95=nan(n_bins,n_bins,1);
prctile_99=nan(n_bins,n_bins,1);

w=8;%1:length(waves)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binning of cells' position
[C, idx] = bin_cells_position(r_i,bin_size,max_length); %C is the table with the bins, and the data is in idx
%Reaccomodate the C matrix for plotting it
C_s=[];
for m=1:n_bins
    clear aux; aux=1:n_bins;
    C_s=[C_s;[aux',(n_bins-m+1)*ones(n_bins,1)]];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of local phases per sequence
count_s=0;
for s=19%17:19
    count_s=count_s+1;
    clear mean_p; mean_p=mean_p_seq(:,s);
    meanp_pos_sh=nan(n_bins,n_bins,N_sh);
    circp_pos_sh=nan(n_bins,n_bins,N_sh);
    %         clear meanp_pos_sh circp_pos_sh;
    for sh=1:N_sh
        clear mean_p_sh; mean_p_sh=mean_p(randperm(N));
        for m=1:size(C,1) % Loop on all spatial bins
            if sum(idx==m)>0
                clear aux aux2; aux=(mean_p_sh(idx==m)); aux2=aux(~isnan(aux));
                if length(aux2)>min_cutoff_n_neurons
                    meanp_pos_sh(C(m,1),C(m,2),sh)=circ_mean(aux2);
                    circp_pos_sh(C(m,1),C(m,2),sh)=circ_r(aux2);
                end
            end
        end
    end

    for m=1:size(C,1) % Loop on all spatial bins
        if sum(idx==m)>0
            clear aux aux2; aux=(mean_p(idx==m)); aux2=aux(~isnan(aux));
            occupancy_mat(C(m,1),C(m,2),count_s)=length(aux2);
            if length(aux2)>min_cutoff_n_neurons
                occupancy_mat_b(C(m,1),C(m,2),count_s)=length(aux2);
                meanp_pos(C(m,1),C(m,2),count_s)=circ_mean(aux2);
                circp_pos(C(m,1),C(m,2),count_s)=circ_r(aux2);

                prctile_90(C(m,1),C(m,2),count_s)=prctile(circp_pos_sh(C(m,1),C(m,2),:),90);
                prctile_95(C(m,1),C(m,2),count_s)=prctile(circp_pos_sh(C(m,1),C(m,2),:),95);
                prctile_99(C(m,1),C(m,2),count_s)=prctile(circp_pos_sh(C(m,1),C(m,2),:),99);
            end
        end
    end

    %Reaccomodate matrix for occupancy and plot
    occupancy_mat_b_plot=occupancy_mat';
    %     figure
    %     imagesc(flip(occupancy_mat_b_plot));
    %     colorbar
    %     set(gca,'fontsize',10);
    %     ylabel('Y bin #');
    %     xlabel('X bin #');
    %     set(gca,'fontsize',16);

    figure
    yvalues=n_bins:-1:1;
    xvalues=1:n_bins;
    heatmap(xvalues,yvalues,flip(occupancy_mat_b_plot));
    %     c=colorbar;
    set(gca,'fontsize',10);
    ylabel('Y bin #');
    xlabel('X bin #');
    set(gca,'fontsize',16);
    %     title(c,'Counts');


    figure
    for m=1:size(C,1) % Loop on all spatial bins
        if isnan(circp_pos(C_s(m,1),C_s(m,2),count_s))
        else
            subplot(n_bins,n_bins,m)
            histogram(circp_pos_sh(C_s(m,1),C_s(m,2),:),[0:1/40:1]);
            title({['Bin = (',num2str(C_s(m,1)),',',num2str(C_s(m,2)),')']});
            axis([0 1 -inf inf])
            hold on
            xline( prctile_95(C_s(m,1),C_s(m,2),count_s),'b-','linewidth',5);
            xline( circp_pos(C_s(m,1),C_s(m,2),count_s),'r-','linewidth',5);
            ylabel('Counts');
            xlabel('MVL');
            set(gca,'fontsize',16);
        end
    end

    %       figure
    %     for m=1:size(C,1) % Loop on all spatial bins
    %         if isnan(circp_pos(C(m,1),C(m,2),count_s))
    %         else
    %             subplot(n_bins,n_bins,m)
    %             histogram(circp_pos_sh(C(m,1),C(m,2),:),[0:1/40:1]);
    %             title({['Bin = (',num2str(C(m,1)),',',num2str(C(m,2)),')']});
    %             axis([0 1 -inf inf])
    %             hold on
    %             xline( prctile_90(C(m,1),C(m,2),count_s),'b--','linewidth',3);
    %             xline( circp_pos(C(m,1),C(m,2),count_s),'k--','linewidth',3);
    %             ylabel('Counts');
    %             xlabel('MVL');
    %             set(gca,'fontsize',10);
    %         end
    %     end
end

%New quantification
tot_largerthanten=find(occupancy_mat_b(:)>10); %spatial bins with more than 10 cells
circp_pos_b=circp_pos(:);
circp_pos_b=circp_pos_b(tot_largerthanten);
prctile_95_b=prctile_95(:);
prctile_95_b=prctile_95_b(tot_largerthanten);

num=find(circp_pos_b>prctile_95_b);
tot=length(num)*100/length(circp_pos_b);
disp(['Prob = ',num2str(tot)])
