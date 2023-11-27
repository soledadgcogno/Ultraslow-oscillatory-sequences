%% One session - PI

N_sh_dist=500;
w=8;

row_w=waves(w);

count=count+1;
mouse=['L',num2str(big_table(row_w,1)),'M',num2str(big_table(row_w,2))];
day=big_table(row_w,3);
s=big_table(row_w,4);
munit=big_table(row_w,5);

clus=10;
disc_phase=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
file_name_dff=[dpath ['DFF_120ms_Do_SNRH','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

load(file_name_anat,'-mat'); %Anatomical information
load(file_name_spk,'-mat'); %Spike times
spikes_d=full(spikes_d_s);
[N,T]=size(spikes_d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PI and locked cells

not_locked=not_locked_all_sessions{w};
locked=locked_all_sessions{w};
PI=PR_all_sessions{w};
PI_locked=PI(locked);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distance in the tissue

pixel_size_new=1.18185;
pixel_size_old=1.78211;

if w>10
    pixel_size=pixel_size_old;
else
    pixel_size=pixel_size_new;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    r_i(i,:)=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
    dist_origin(i)=norm(r_i(i,:));
end
r_i(not_locked,:)=[];
dist_origin_locked=dist_origin(locked);
max_x=max(r_i(:,1));
max_y=max(r_i(:,2));
max_xy=max(max_x,max_y);

edges_x=0:50:ceil(max_xy)+50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=0;
for i=1:N
    r_i=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
    for j=i+1:N
        count=count+1;
        r_j=[Anat.med{1,j}(1)*pixel_size,Anat.med{1,j}(2)*pixel_size];
        delta_tissue(i,j)=norm(r_j-r_i);
        delta_tissue(j,i)=norm(r_j-r_i);
        delta_tissue_vec(count)=norm(r_j-r_i);
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
        delta_PI_mean(i,j)=PI(i)-PI(j);
        delta_PI_mean(j,i)=PI(j)-PI(i);
    end
end
delta_PI_mean(not_locked,:)=[];
delta_PI_mean(:,not_locked)=[];

matrix_distance_tissue{w}=delta_tissue;
matrix_distance_phase{w}=delta_PI_mean;

diff_PI_reshape=delta_PI_mean(:);
diff_pos_reshape=delta_tissue(:);

clear dist_diff p_sim p_diff dist_sim
pairs_sim=find(abs(diff_PI_reshape)<0.001);
pairs_diff=find(abs(diff_PI_reshape)>0.25);

dist_sim=diff_pos_reshape(pairs_sim);
dist_diff=diff_pos_reshape(pairs_diff);

p_sim=histcounts(dist_sim,edges_x,'Normalization','Probability');
p_diff=histcounts(dist_diff,edges_x,'Normalization','Probability');

centroid=edges_x-25;
centroid(1)=[];
figure
plot(centroid,cumsum(p_sim),'linewidth',2.5);
hold on
plot(centroid,cumsum(p_diff),'linewidth',2.5);
ylabel('Cumulative probability');
xlabel('Pairwise distance [um]');
set(gca,'fontsize',16);
legend('Similar preferred phase','Differente preferred phase');
title('Data - PI');
legend boxoff 

% [h_small(w),p_small(w)] = kstest2(dist_sim,dist_diff,'Tail','smaller') ;
% [h_large(w),p_large(w)] = kstest2(dist_sim,dist_diff,'Tail','larger') ;
% [h_unequal(w),p_unequal(w)] = kstest2(dist_sim,dist_diff,'Tail','unequal') ;
% [h_default(w),p_default(w)] = kstest2(dist_sim,dist_diff) ;

KL1_one=kldiv(p_sim,p_diff);
KL2_one=kldiv(p_diff,p_sim);

%shuffle

for sh=1:N_sh_dist
   
    PI_sh=PI_locked(randperm(length(PI_locked)));
    for i=1:length(locked)
        for j=i+1:length(locked)
            delta_PI_mean_sh(i,j)=PI_sh(i)-PI_sh(j);
            delta_PI_mean_sh(j,i)=PI_sh(j)-PI_sh(i);
        end
    end

    diff_PI_reshape_sh=delta_PI_mean_sh(:);

    pairs_sim_sh=find(abs(diff_PI_reshape_sh)<0.001);
    pairs_diff_sh=find(abs(diff_PI_reshape_sh)>0.25);

    dist_sim_sh=diff_pos_reshape(pairs_sim_sh);
    dist_diff_sh=diff_pos_reshape(pairs_diff_sh);

    p_sim_sh=histcounts(dist_sim_sh,edges_x,'Normalization','Probability');
    p_diff_sh=histcounts(dist_diff_sh,edges_x,'Normalization','Probability');

    if sh==1
        PI_shuffled=PI_sh;
        figure
        plot(centroid,cumsum(p_sim_sh),'linewidth',2.5);
        hold on
        plot(centroid,cumsum(p_diff_sh),'linewidth',2.5);
        ylabel('Cumulative probability');
        xlabel('Pairwise distance [um]');
        set(gca,'fontsize',16);
        legend('Similar preferred phase','Differente preferred phase');
        title('Shuffle - PI')
        legend boxoff 
    end
    
    [h_sh_one(sh),p_sh_one(sh)] = kstest2(dist_sim_sh,dist_diff_sh,'Tail','smaller') ;
    KL1_sh_one(sh)=kldiv(p_sim_sh,p_diff_sh);
    KL2_sh_one(sh)=kldiv(p_diff_sh,p_sim_sh);
    
    clear p_sim_sh p_diff_sh pairs_sim_sh pairs_diff_sh dist_sim_sh dist_diff_sh delta_phase_mean_sh ...
            diff_phase_reshape_sh
end
% % % %%%%%%% Order
% % % 
% % % for i=1:N
% % %     r_i=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
% % %     dist_origin(i)=norm(r_i);
% % % end
% % % dist_origin_locked=dist_origin(locked);
% % % clear a b c d mean_phase_locked_ordered
% % % [a,b]=sort(dist_origin_locked,'ascend');
% % % [c,d]=sort(mean_p_locked,'ascend');
% % % 
% % % mean_phase_locked_ordered=nan(1,length(locked));
% % % for i=1:length(locked)
% % %     mean_phase_locked_ordered(b(i))=c(i);
% % % end
% % % 
% % % for i=1:length(locked)
% % %     for j=i+1:length(locked)
% % %         delta_phase_mean_or(i,j)=angdiff(mean_phase_locked_ordered(i),mean_phase_locked_ordered(j));
% % %         delta_phase_mean_or(j,i)=angdiff(mean_phase_locked_ordered(j),mean_phase_locked_ordered(i));
% % %     end
% % % end
% % % 
% % % diff_phase_reshape_or=delta_phase_mean_or(:);
% % % 
% % % %     clear dist_diff p_sim p_diff dist_sim
% % % pairs_sim_or=find(abs(diff_phase_reshape_or)<delta_phase);
% % % pairs_diff_or=find(abs(diff_phase_reshape_or)>pi-delta_phase);
% % % 
% % % dist_sim_or=diff_pos_reshape(pairs_sim_or);
% % % dist_diff_or=diff_pos_reshape(pairs_diff_or);
% % % 
% % % p_sim_or=histcounts(dist_sim_or,edges_x,'Normalization','Probability');
% % % p_diff_or=histcounts(dist_diff_or,edges_x,'Normalization','Probability');
% % % 
% % % figure
% % % plot(centroid,cumsum(p_sim_or),'linewidth',2.5);
% % % hold on
% % % plot(centroid,cumsum(p_diff_or),'linewidth',2.5);
% % % ylabel('Cumulative probability');
% % % xlabel('Distance between pairs of cells [um]');
% % % set(gca,'fontsize',16);
% % % legend('Similar preferred phase','Differente preferred phase');
% % % title('Ordered');
% % % 
% % % [h_small_or_one,p_small_or_one] = kstest2(dist_sim_or,dist_diff_or,'Tail','smaller') ;
% % % KL2_or_one=kldiv(p_diff_or,p_sim_or);
% % % KL1_or_one=kldiv(p_sim_or,p_diff_or);

% Figures of maps

figure
count_map=0;
[b,E]=discretize(PI_locked,40);    % [~,b2]=sort(b,'ascend');
cc=winter(40);
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_anat,'-mat');
figure
hold on
pixel_size_new=1.18185;
for i=1:N
    x=Anat.pos{1,i}(:,1);
    y=Anat.pos{1,i}(:,2);
    ov=Anat.overlap{1,i};
    
    if ismember(i,locked)==1
        count_map=count_map+1;
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(b(count_map),:),'filled');
    else
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
        
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
colormap winter(40)
co=colorbar('XTick',0:1,'XTickLabel',{'0.24','0.7'});
co.Label.String = 'PI';

figure
count_map=0;
[b,E]=discretize(PI_shuffled,40);    % [~,b2]=sort(b,'ascend');
cc=winter(40);
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_anat,'-mat');
figure
hold on
pixel_size_new=1.18185;
for i=1:N
    x=Anat.pos{1,i}(:,1);
    y=Anat.pos{1,i}(:,2);
    ov=Anat.overlap{1,i};
    
    if ismember(i,locked)==1
        count_map=count_map+1;
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(b(count_map),:),'filled');
    else
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
        
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
colormap winter(40)
co=colorbar('XTick',0:1,'XTickLabel',{'-\pi','\pi'});
co.Label.String = 'PI';


[Y,edges_kl]=histcounts(KL2_sh_one)
centrod_kl=edges_kl-(edges_kl(2)-edges_kl(1))/2;
centrod_kl(1)=[];

figure
bar(centrod_kl,Y)
hold on
ylabel('Counts');
xlabel('KL divergence');
set(gca,'fontsize',16)
hold on
l=xline(KL2_one);
l.LineStyle='--';
l.LineWidth=2.5;
l.Color='k';
box off
l=xline(prctile(KL2_sh_one,95));
l.LineStyle='-.';
l.LineWidth=2.5;
l.Color='k';
legend
legend boxoff

figure
bar(centrod_kl,Y)
hold on
ylabel('Counts');
xlabel('KL divergence');
set(gca,'fontsize',16)
hold on
l=xline(KL2_one);
l.LineStyle='--';
l.LineWidth=2.5;
l.Color='k';
box off
l=xline(prctile(KL2_sh_one,99));
l.LineStyle='-.';
l.LineWidth=2.5;
l.Color='k';
legend
legend boxoff


clear diff_phase_reshape diff_pos_reshape delta_phase_mean delta_tissue spikes_d not_locked locked mean_p mean_p_locked prob_join dist_origin_locked dist_origin
clear mean_phase_locked_ordered a b c d diff_phase_ordered_reshape  delta_phase_mean_ordered x y ov Anat ...
    p_sim p_diff pairs_sim pairs_diff dist_sim dist_diff delta_phase_mean diff_phase_reshape diff_pos_reshape delta_tissue r_i ...
    locked not_locked mean_p mean_p_locked dist_origin...
    p_sim_or p_diff_or dist_sim_or dist_diff_or pairs_sim_or pairs_diff_or diff_phase_reshape_or mean_phase_locked_ordered ...
    a b c d max_x max_y max_xy mean_p_sh N spikes_d delta_tissue edges_x dist_sim_or dist_diff_or pairs_sim_or pairs_diff_or  ...
    delta_phase_mean_or r_i dist_origin_locked dist_origin p_sim_or p_diff_or


%% All sessions - PI

N_sh_dist=500;
for w=1:length(waves)
    
    row_w=waves(w);
    disp(w)
    
    count=count+1;
    mouse=['L',num2str(big_table(row_w,1)),'M',num2str(big_table(row_w,2))];
    day=big_table(row_w,3);
    s=big_table(row_w,4);
    munit=big_table(row_w,5);
    
    clus=10;
    disc_phase=10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_dff=[dpath ['DFF_120ms_Do_SNRH','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    
    load(file_name_anat,'-mat'); %Anatomical information
    load(file_name_spk,'-mat'); %Spike times
    spikes_d=full(spikes_d_s);
    [N,T]=size(spikes_d);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Phases and locked cells
    
    not_locked=not_locked_all_sessions{w};
    locked=locked_all_sessions{w};
    PI=PR_all_sessions{w};
    PI_locked=PI(locked);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distance in the tissue
    
    pixel_size_new=1.18185;
    pixel_size_old=1.78211;
    
    if w>10
        pixel_size=pixel_size_old;
    else
        pixel_size=pixel_size_new;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:N
        r_i(i,:)=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
        dist_origin(i)=norm(r_i);
    end
    r_i(not_locked,:)=[];
    dist_origin_locked=dist_origin(locked);   
    max_x=max(r_i(:,1));
    max_y=max(r_i(:,2));
    max_xy=max(max_x,max_y);    
 
    edges_x=0:50:ceil(max_xy)+50;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    count=0;
    for i=1:N
        r_i=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
        dist_origin(i)=norm(r_i);
        for j=i+1:N
            count=count+1;
            r_j=[Anat.med{1,j}(1)*pixel_size,Anat.med{1,j}(2)*pixel_size];
            delta_tissue(i,j)=norm(r_j-r_i);
            delta_tissue(j,i)=norm(r_j-r_i);
            delta_tissue_vec(count)=norm(r_j-r_i);
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
            delta_PI_mean(i,j)=PI(i)-PI(j);
            delta_PI_mean(j,i)=PI(j)-PI(i);
        end
    end
    delta_PI_mean(not_locked,:)=[];
    delta_PI_mean(:,not_locked)=[];
    
    matrix_distance_tissue{w}=delta_tissue;
    matrix_distance_PI{w}=delta_PI_mean;
    
    diff_PI_reshape=delta_PI_mean(:);
    diff_pos_reshape=delta_tissue(:);
    
    
    clear dist_diff p_sim p_diff dist_sim
    pairs_sim=find(abs(diff_PI_reshape)<0.001);
    pairs_diff=find(abs(diff_PI_reshape)>0.25);
    
    dist_sim=diff_pos_reshape(pairs_sim);
    dist_diff=diff_pos_reshape(pairs_diff);
    
    p_sim=histcounts(dist_sim,edges_x,'Normalization','Probability');
    p_diff=histcounts(dist_diff,edges_x,'Normalization','Probability');
    
    p_sim_w{w}=p_sim;
    p_sim_w{w}=p_diff;
    
%     figure
%     plot(cumsum(p_sim))
%     hold on
%     plot(cumsum(p_diff))
    
    [p_rank(w),h_rank(w)] = ranksum(dist_sim,dist_diff,'Tail','left') ;

    [h_small(w),p_small(w)] = kstest2(dist_sim,dist_diff,'Tail','smaller') ;
    [h_large(w),p_large(w)] = kstest2(dist_sim,dist_diff,'Tail','larger') ;
    [h_unequal(w),p_unequal(w)] = kstest2(dist_sim,dist_diff,'Tail','unequal') ;
    [h_default(w),p_default(w)] = kstest2(dist_sim,dist_diff) ;

    KL1(w)=kldiv(p_sim,p_diff);
    KL2(w)=kldiv(p_diff,p_sim);

    %shuffle
    for sh=1:N_sh_dist
        PI_sh=PI_locked(randperm(length(PI_locked)));
        for i=1:length(locked)
            for j=i+1:length(locked)
                delta_PI_mean_sh(i,j)=PI_sh(i)-PI_sh(j);
                delta_PI_mean_sh(j,i)=PI_sh(j)-PI_sh(i);
            end
        end
       
        diff_PI_reshape_sh=delta_PI_mean_sh(:);
                
        pairs_sim_sh=find(abs(diff_PI_reshape_sh)<0.001);
        pairs_diff_sh=find(abs(diff_PI_reshape_sh)>0.25);
        
        dist_sim_sh=diff_pos_reshape(pairs_sim_sh);
        dist_diff_sh=diff_pos_reshape(pairs_diff_sh);
        
        p_sim_sh=histcounts(dist_sim_sh,edges_x,'Normalization','Probability');
        p_diff_sh=histcounts(dist_diff_sh,edges_x,'Normalization','Probability');
        
%         figure
%         plot(cumsum(p_sim_sh))
%         hold on
%         plot(cumsum(p_diff_sh))
        
%         [h_sh(w,sh),p_sh(w,sh)] = kstest2(dist_sim_sh,dist_diff_sh,'Tail','smaller') ;
        KL1_sh(w,sh)=kldiv(p_sim_sh,p_diff_sh);
        KL2_sh(w,sh)=kldiv(p_diff_sh,p_sim_sh);

        if sh==1
            p_sim_w_sh{w}=p_sim_sh;
            p_sim_w_sh{w}=p_diff_sh;
        end
        clear p_sim_sh p_diff_sh pairs_sim_sh pairs_diff_sh dist_sim_sh dist_diff_sh delta_PI_mean_sh ...
            diff_PI_reshape_sh
    end
    
    %%%%%%% Order
%     
%     [a,b]=sort(dist_origin,'ascend');
%     [c,d]=sort(mean_p,'ascend');
%     mean_PI_locked_ordered=nan(1,length(locked));
%     for i=1:length(locked)
%         mean_PI_locked_ordered(b(i))=c(i);
%     end  
%     mean_PI_locked_ordered_full{w}=mean_PI_locked_ordered;
%     
%     clear a b c d mean_PI_locked_ordered
%     [a,b]=sort(dist_origin(locked),'ascend');
%     [c,d]=sort(mean_p_locked,'ascend');
%     
%     mean_PI_locked_ordered=nan(1,length(locked));
%     for i=1:length(locked)
%         mean_PI_locked_ordered(b(i))=c(i);
%     end    
%     
%     for i=1:length(locked)
%         for j=i+1:length(locked)
%             delta_PI_mean_or(i,j)=angdiff(mean_PI_locked_ordered(i),mean_PI_locked_ordered(j));
%             delta_PI_mean_or(j,i)=angdiff(mean_PI_locked_ordered(j),mean_PI_locked_ordered(i));
%         end
%     end
%      
%     diff_PI_reshape_or=delta_PI_mean_or(:);
%        
% %     clear dist_diff p_sim p_diff dist_sim
%     pairs_sim_or=find(abs(diff_PI_reshape_or)<delta_PI);
%     pairs_diff_or=find(abs(diff_PI_reshape_or)>pi-delta_PI);
%     
%     dist_sim_or=diff_pos_reshape(pairs_sim_or);
%     dist_diff_or=diff_pos_reshape(pairs_diff_or);
%     
%     p_sim_or=histcounts(dist_sim_or,edges_x,'Normalization','Probability');
%     p_diff_or=histcounts(dist_diff_or,edges_x,'Normalization','Probability');
%     
% %     figure
% %     plot(cumsum(p_sim_or))
% %     hold on
% %     plot(cumsum(p_diff_or))
%     
%     [h_small_or(w),p_small_or(w)] = kstest2(dist_sim_or,dist_diff_or,'Tail','smaller') ;
%     KL2_or(w)=kldiv(p_diff_or,p_sim_or);
%     KL1_or(w)=kldiv(p_sim_or,p_diff_or);
%     
%     p_sim_w_or{w}=p_sim_or;
%     p_diff_w_or{w}=p_diff_or;
% 
%    
    clear diff_PI_reshape diff_pos_reshape delta_PI_mean delta_tissue spikes_d not_locked locked mean_p mean_p_locked prob_join dist_origin_locked dist_origin
    clear mean_PI_locked_ordered a b c d diff_PI_ordered_reshape  delta_PI_mean_ordered x y ov Anat ...
        p_sim p_diff pairs_sim pairs_diff dist_sim dist_diff delta_PI_mean diff_PI_reshape diff_pos_reshape delta_tissue r_i ...
        locked not_locked mean_p mean_p_locked dist_origin...
        p_sim_or p_diff_or dist_sim_or dist_diff_or pairs_sim_or pairs_diff_or diff_PI_reshape_or mean_PI_locked_ordered ...
        a b c d max_x max_y max_xy mean_p_sh N spikes_d delta_tissue edges_x dist_sim_or dist_diff_or pairs_sim_or pairs_diff_or  ...
        delta_PI_mean_or r_i dist_origin_locked dist_origin p_sim_or p_diff_or PI PI_locked
end

for w=1:15
    
thr_999(w)= prctile(KL2_sh(w,~isinf(KL2_sh(w,:))),99.9);
thr_99(w)= prctile(KL2_sh(w,~isinf(KL2_sh(w,:))),99);
thr_95(w)= prctile(KL2_sh(w,~isinf(KL2_sh(w,:))),95);
end
find(KL2>thr_999)
find(KL2>thr_99)
find(KL2>thr_95)

figure
scatter(KL2,thr_95,40,'o','filled')
h=refline(1,0);
h.LineStyle='--';
h.LineWidth=2.5;
h.Color='k';
ylabel('95th percentile of null distribution');
xlabel('KL divergence - data');
set(gca,'fontsize',16)
title('PI')

figure
scatter(KL2,thr_99,40,'o','filled')
h=refline(1,0);
h.LineStyle='--';
h.LineWidth=2.5;
h.Color='k';
ylabel('99th percentile of null distribution');
xlabel('KL divergence - data');
set(gca,'fontsize',16)
title('PI')

[h,c]=histcounts(p_rank,0:0.05:1,'Normalization','cdf');
