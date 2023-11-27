%% Load data

clear all
close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

mice_number=12;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];

ml=[2.64,3.7,2.7,2.64,3.2,3.3,2.52,3.3,3.5,-10,-10,-10]; %ML coordinates

count=0;
for m=1:mice_number
    
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
    load([save_data_path ['WS_Osc_15_sf7p73II_',mice(m,:),'.mat']]);
    ws_ent=load([save_data_path ['WS_Entropy_',mice(m,:),'.mat']]);
 
    for day=1:dates.daysnum
        for s=1:dates.sesnum(day)
            disp(s)
            munit=dates.ses{day}(s);
            
            file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
            
            if (exist(file_name_spk) == 2)
                disp(day)
                spk=load(file_name_spk);
                
                if s<100
                    count=count+1;
                    
                    if isfield(dates,'actual_day')  == 1
                        big_table(count,3)=dates.actual_day(day);   %Day # on the wheel
                    else
                        big_table(count,3)=day; %Day # on the wheel
                    end
                    
                    if mouse(1) == '9' || mouse(1) == '6'
                        big_table(count,1)=1;%str2num(mice(m,2:3));
                        big_table(count,10)=-10; %ML position
                        if mouse=='92227'
                            big_table(count,2)=1;
                        elseif mouse == '92229'
                            big_table(count,2)=2;
                        elseif mouse == '60961'
                            big_table(count,2)=3;
                        end
                    else
                        big_table(count,1)=str2num(mice(m,2:3)); %Litter
                        big_table(count,2)=str2num(mice(m,end)); %Mouse number
                        big_table(count,10)=ml(m); %ML position
                    end
                    
                    big_table(count,4)=s; %Session
                    big_table(count,5)=munit; %munit
                    big_table(count,6)=WS_stat.WS(day,s); %WS - Osc
                    big_table(count,11)=ws_ent.WS_stat.wave_score_ent(day,s); %WS - Ent
                    big_table(count,12)=size(spk.spikes_d_s,1);
                    
                    threshold_kl=1;
                    if ( big_table(count,6)>=threshold_kl )
                        big_table(count,7)=1;
                    else
                        big_table(count,7)=0;
                    end
                    
                    if (isnan(big_table(count,6)))
                        big_table(count,7)=NaN;
                        big_table(count,8)=NaN; %dt
                        big_table(count,9)=NaN; %fft
                    else
                        big_table(count,8)=WS_stat.dt(day,s); %dt
                        big_table(count,9)=WS_stat.FFT(day,s); %fft
                    end
                    
                end
               clear spk 
            end
        end
        
    end
    clear WS_stat 
end

% Wave sessions

% figpath='C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Raster Plots all MEC sessions\';
MEC_ses=find(big_table(:,10)>3); %Sessions inthe in MEC
adults=find(big_table(:,3)>15); %Sessions inthe of adults
sessions=intersect(MEC_ses,adults);
count=0;

ws=big_table(sessions,6);
index_nan=isnan(ws);
sessions(index_nan)=[];
session_number=big_table(sessions,4);
days=big_table(sessions,3);
repeated_days=find(diff(days)==0);
clear index_to_delete

count=0;
for w=1:length(repeated_days)
    ws1=big_table(sessions(repeated_days(w)),6);
    ws2=big_table(sessions(repeated_days(w)+1),6);
    
    if (ws1==1 && ws2==0)
        count=count+1;
        index_to_delete(count)=repeated_days(w)+1;
    elseif  (ws1==0 && ws2==1)
        count=count+1;
        index_to_delete(count)=repeated_days(w);
    elseif (ws1==0 && ws2==0)
        count=count+1;
        index_to_delete(count)=repeated_days(w);
    else
        print('Other case');
    end   
end

sessions(index_to_delete)=[];
big_table_2=big_table(sessions,:);
mec_sessions=sessions;
waves=mec_sessions(find(big_table(mec_sessions,6)==1));




% PI one session
load('locking_all_sessions.mat')
w=8;
N_sh_dist=1;
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

[a,b]=sort(abs(diff_PI_reshape),'ascend');

diff_PI_reshape_s=abs(diff_PI_reshape);
diff_PI_reshape_s=diff_PI_reshape_s(b);
diff_pos_reshape_s=diff_pos_reshape(b);

% thr=median(abs(diff_PI_reshape));
% pairs_sim=find(abs(diff_PI_reshape)<prctile(abs(diff_PI_reshape),25));
% pairs_diff=find(abs(diff_PI_reshape)>prctile(abs(diff_PI_reshape),75));
% 
% dist_sim=diff_pos_reshape(pairs_sim);
% dist_diff=diff_pos_reshape(pairs_diff);

thr=floor(0.2*length(diff_pos_reshape_s));
dist_sim=diff_pos_reshape_s(1:floor(thr));
dist_diff=diff_pos_reshape_s(end-floor(thr):end);

% dist_sim=diff_pos_reshape(pairs_sim);
% dist_diff=diff_pos_reshape(pairs_diff);

% p_sim=histcounts(dist_sim,edges_x,'Normalization','Probability');
% p_diff=histcounts(dist_diff,edges_x,'Normalization','Probability');

% centroid=edges_x-25;
% centroid(1)=[];
% figure
% plot(centroid,cumsum(p_sim),'linewidth',2.5);
% hold on
% plot(centroid,cumsum(p_diff),'linewidth',2.5);
% ylabel('Cumulative probability');
% xlabel('Pairwise distance [um]');
% set(gca,'fontsize',16);
% legend('Similar preferred phase','Differente preferred phase');
% title('Data - PI');
% legend boxoff 
% 
% % [h_small(w),p_small(w)] = kstest2(dist_sim,dist_diff,'Tail','smaller') ;
% % [h_large(w),p_large(w)] = kstest2(dist_sim,dist_diff,'Tail','larger') ;
% % [h_unequal(w),p_unequal(w)] = kstest2(dist_sim,dist_diff,'Tail','unequal') ;
% % [h_default(w),p_default(w)] = kstest2(dist_sim,dist_diff) ;
% 
% KL1_one=kldiv(p_sim,p_diff);
% KL2_one=kldiv(p_diff,p_sim);
% 
% %shuffle
% 
% for sh=1:N_sh_dist
%    
%     PI_sh=PI_locked(randperm(length(PI_locked)));
%     for i=1:length(locked)
%         for j=i+1:length(locked)
%             delta_PI_mean_sh(i,j)=PI_sh(i)-PI_sh(j);
%             delta_PI_mean_sh(j,i)=PI_sh(j)-PI_sh(i);
%         end
%     end
% 
%     diff_PI_reshape_sh=delta_PI_mean_sh(:);
% 
%     pairs_sim_sh=find(abs(diff_PI_reshape_sh)<0.001);
%     pairs_diff_sh=find(abs(diff_PI_reshape_sh)>0.3);
% 
%     dist_sim_sh=diff_pos_reshape(pairs_sim_sh);
%     dist_diff_sh=diff_pos_reshape(pairs_diff_sh);
% 
%     p_sim_sh=histcounts(dist_sim_sh,edges_x,'Normalization','Probability');
%     p_diff_sh=histcounts(dist_diff_sh,edges_x,'Normalization','Probability');
% 
%     if sh==1
%         PI_shuffled=PI_sh;
%         figure
%         plot(centroid,cumsum(p_sim_sh),'linewidth',2.5);
%         hold on
%         plot(centroid,cumsum(p_diff_sh),'linewidth',2.5);
%         ylabel('Cumulative probability');
%         xlabel('Pairwise distance [um]');
%         set(gca,'fontsize',16);
%         legend('Similar preferred phase','Differente preferred phase');
%         title('Shuffle - PI')
%         legend boxoff 
%     end
%     
%     [h_sh_one(sh),p_sh_one(sh)] = kstest2(dist_sim_sh,dist_diff_sh,'Tail','smaller') ;
%     KL1_sh_one(sh)=kldiv(p_sim_sh,p_diff_sh);
%     KL2_sh_one(sh)=kldiv(p_diff_sh,p_sim_sh);
%     
%     clear p_sim_sh p_diff_sh pairs_sim_sh pairs_diff_sh dist_sim_sh dist_diff_sh delta_phase_mean_sh ...
%             diff_phase_reshape_sh
% end
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

% % Figures of maps
% 
% figure
% count_map=0;
% [b,E]=discretize(PI_locked,40);    % [~,b2]=sort(b,'ascend');
% cc=winter(40);
% file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
% load(file_name_anat,'-mat');
% figure
% hold on
% pixel_size_new=1.18185;
% for i=1:N
%     x=Anat.pos{1,i}(:,1);
%     y=Anat.pos{1,i}(:,2);
%     ov=Anat.overlap{1,i};
%     
%     if ismember(i,locked)==1
%         count_map=count_map+1;
%         scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(b(count_map),:),'filled');
%     else
%         scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
%         
%     end
%     clear x y ov
% end
% box on
% xlabel('X [\mum]');
% ylabel('Y [\mum]');
% axis([0 580 20 600])
% yticks([120 420]);
% yticklabels([100 400]);
% xticks([100 400]);
% set(gca,'fontsize',16)
% axis square
% colormap winter(40)
% co=colorbar('XTick',0:1,'XTickLabel',{'0.24','0.7'});
% co.Label.String = 'PI';
% 
% figure
% count_map=0;
% [b,E]=discretize(PI_shuffled,40);    % [~,b2]=sort(b,'ascend');
% cc=winter(40);
% file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
% load(file_name_anat,'-mat');
% figure
% hold on
% pixel_size_new=1.18185;
% for i=1:N
%     x=Anat.pos{1,i}(:,1);
%     y=Anat.pos{1,i}(:,2);
%     ov=Anat.overlap{1,i};
%     
%     if ismember(i,locked)==1
%         count_map=count_map+1;
%         scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(b(count_map),:),'filled');
%     else
%         scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
%         
%     end
%     clear x y ov
% end
% box on
% xlabel('X [\mum]');
% ylabel('Y [\mum]');
% axis([0 580 20 600])
% yticks([120 420]);
% yticklabels([100 400]);
% xticks([100 400]);
% set(gca,'fontsize',16)
% axis square
% colormap winter(40)
% co=colorbar('XTick',0:1,'XTickLabel',{'-\pi','\pi'});
% co.Label.String = 'PI';
% 
% 
% [Y,edges_kl]=histcounts(KL2_sh_one)
% centrod_kl=edges_kl-(edges_kl(2)-edges_kl(1))/2;
% centrod_kl(1)=[];
% 
% figure
% bar(centrod_kl,Y)
% hold on
% ylabel('Counts');
% xlabel('KL divergence');
% set(gca,'fontsize',16)
% hold on
% l=xline(KL2_one);
% l.LineStyle='--';
% l.LineWidth=2.5;
% l.Color='k';
% box off
% l=xline(prctile(KL2_sh_one,95));
% l.LineStyle='-.';
% l.LineWidth=2.5;
% l.Color='k';
% legend
% legend boxoff
% 
% figure
% bar(centrod_kl,Y)
% hold on
% ylabel('Counts');
% xlabel('KL divergence');
% set(gca,'fontsize',16)
% hold on
% l=xline(KL2_one);
% l.LineStyle='--';
% l.LineWidth=2.5;
% l.Color='k';
% box off
% l=xline(prctile(KL2_sh_one,99));
% l.LineStyle='-.';
% l.LineWidth=2.5;
% l.Color='k';
% legend
% legend boxoff

mat_PI_distance=nan(max([length(dist_sim),length(dist_diff)]),2);
mat_PI_distance(1:length(dist_sim),1)=dist_sim;
mat_PI_distance(1:length(dist_diff),2)=dist_diff;



[p,h,stat]=ranksum(dist_sim,dist_diff);
figure
boxplot(mat_PI_distance)
xticklabels({'Similar PI','Different PI'});
ylabel('Pairwise anatomical distance (um)');
set(gca,'fontsize',16,'YColor','k','XColor','k');
ylim([-5 1000]);
box off
title(num2str(p))
% 
% clear diff_phase_reshape diff_pos_reshape delta_phase_mean delta_tissue spikes_d not_locked locked mean_p mean_p_locked prob_join dist_origin_locked dist_origin
% clear mean_phase_locked_ordered a b c d diff_phase_ordered_reshape  delta_phase_mean_ordered x y ov Anat ...
%     p_sim p_diff pairs_sim pairs_diff dist_sim dist_diff delta_phase_mean diff_phase_reshape diff_pos_reshape delta_tissue r_i ...
%     locked not_locked mean_p mean_p_locked dist_origin...
%     p_sim_or p_diff_or dist_sim_or dist_diff_or pairs_sim_or pairs_diff_or diff_phase_reshape_or mean_phase_locked_ordered ...
%     a b c d max_x max_y max_xy mean_p_sh N spikes_d delta_tissue edges_x dist_sim_or dist_diff_or pairs_sim_or pairs_diff_or  ...
%     delta_phase_mean_or r_i dist_origin_locked dist_origin p_sim_or p_diff_or


%% All sessions - PI


delta_PI(:,1)=[0.001,0.001,0.001,0.05,0.05,0.05,0.1,0.1,0.1,0.15,0.15,0.15];
delta_PI(:,2)=[0.3,0.25,0.2,0.3,0.25,0.2,0.3,0.25,0.2,0.3,0.25,0.2];

load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\locking_all_sessions.mat');

N_sh_dist=5;
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
    matrix_distance_PI{w}=delta_PI_mean;

    diff_PI_reshape=delta_PI_mean(:);
    diff_pos_reshape=delta_tissue(:);

    clear dist_diff p_sim p_diff dist_sim

    count=0;
    for count=1:size(delta_PI,1)

        pairs_sim=find(abs(diff_PI_reshape)<delta_PI(count,1));
        pairs_diff=find(abs(diff_PI_reshape)>delta_PI(count,2));

        num_sim(w,count)=length(pairs_sim);
        num_diff(w,count)=length(pairs_diff);

        dist_sim=diff_pos_reshape(pairs_sim);
        dist_diff=diff_pos_reshape(pairs_diff);

        p_sim=histcounts(dist_sim,edges_x,'Normalization','Probability');
        p_diff=histcounts(dist_diff,edges_x,'Normalization','Probability');

        p_sim_w{w,count}=p_sim;
        p_diff_w{w,count}=p_diff;

        mean_p_sim_w(w,count)=mean(dist_sim);
        mean_p_diff_w(w,count)=mean(dist_diff);

        [p_rank(w),h_rank(w)] = ranksum(dist_sim,dist_diff,'Tail','left') ;

        [h_small(w),p_small(w)] = kstest2(dist_sim,dist_diff,'Tail','smaller') ;
        [h_large(w),p_large(w)] = kstest2(dist_sim,dist_diff,'Tail','larger') ;
        [h_unequal(w),p_unequal(w)] = kstest2(dist_sim,dist_diff,'Tail','unequal') ;
        [h_default(w),p_default(w)] = kstest2(dist_sim,dist_diff) ;

        KL1(w,count)=kldiv(p_sim,p_diff);
        KL2(w,count)=kldiv(p_diff,p_sim);

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

            pairs_sim_sh=find(abs(diff_PI_reshape_sh)<delta_PI(count,1));
            pairs_diff_sh=find(abs(diff_PI_reshape_sh)>delta_PI(count,2));

            dist_sim_sh=diff_pos_reshape(pairs_sim_sh);
            dist_diff_sh=diff_pos_reshape(pairs_diff_sh);

            p_sim_sh=histcounts(dist_sim_sh,edges_x,'Normalization','Probability');
            p_diff_sh=histcounts(dist_diff_sh,edges_x,'Normalization','Probability');

            %         figure
            %         plot(cumsum(p_sim_sh))
            %         hold on
            %         plot(cumsum(p_diff_sh))

            %         [h_sh(w,sh),p_sh(w,sh)] = kstest2(dist_sim_sh,dist_diff_sh,'Tail','smaller') ;
            KL1_sh(w,sh,count)=kldiv(p_sim_sh,p_diff_sh);
            KL2_sh(w,sh,count)=kldiv(p_diff_sh,p_sim_sh);

            if sh==1
                p_sim_w_sh{w}=p_sim_sh;
                p_sim_w_sh{w}=p_diff_sh;
            end
            clear p_sim_sh p_diff_sh pairs_sim_sh pairs_diff_sh dist_sim_sh dist_diff_sh delta_PI_mean_sh ...
                diff_PI_reshape_sh
        end
        clear pairs_sim pairs_diff p_sim p_diff dist_sim dist_diff
    end

    clear diff_PI_reshape diff_pos_reshape delta_PI_mean delta_tissue spikes_d not_locked locked mean_p mean_p_locked prob_join dist_origin_locked dist_origin
    clear mean_PI_locked_ordered a b c d diff_PI_ordered_reshape  delta_PI_mean_ordered x y ov Anat ...
        p_sim p_diff pairs_sim pairs_diff dist_sim dist_diff delta_PI_mean diff_PI_reshape diff_pos_reshape delta_tissue r_i ...
        locked not_locked mean_p mean_p_locked dist_origin...
        p_sim_or p_diff_or dist_sim_or dist_diff_or pairs_sim_or pairs_diff_or diff_PI_reshape_or mean_PI_locked_ordered ...
        a b c d max_x max_y max_xy mean_p_sh N spikes_d delta_tissue edges_x dist_sim_or dist_diff_or pairs_sim_or pairs_diff_or  ...
        delta_PI_mean_or r_i dist_origin_locked dist_origin p_sim_or p_diff_or PI PI_locked
end

% for i=1:count
%     for w=1:15
%         thr_999(w,i)= prctile(KL2_sh(w,~isinf(KL2_sh(w,:,i)),i),99.9);
%         thr_99(w,i)= prctile(KL2_sh(w,~isinf(KL2_sh(w,:,i)),i),99);
%         thr_95(w,i)= prctile(KL2_sh(w,~isinf(KL2_sh(w,:,i)),i),95);
%     end
% end

figure
boxplot([mean_p_sim_w(:,2),mean_p_diff_w(:,2)]);
ylim([0 1000])
xticklabels({'Similar PI','Different PI'});
ylabel('Pairwise anatomical distance (um)');
title('All sessions');
set(gca,'FontSize',16,'XColor','k','Ycolor','k');
box off

[p,h,stat]=ranksum(mean_p_sim_w(:,2),mean_p_diff_w(:,2));

for i=1:size(mean_p_sim_w,2)
    [p(i),h(i)]=ranksum(mean_p_sim_w(:,i),mean_p_diff_w(:,i));
end

figure
plot(1:length(delta_PI),p,'k-*','linewidth',2);
hold on
yline(0.05,'--','LineWidth',2,'Color',[47,79,79]/255);
ylabel('p value');
xlabel('Threshold (\delta/\mu) ');
xticks([1:12])
xlim([0 13])
xticklabels({'(0.001/0.3)','(0.001/0.25)','(0.001/0.2)','(0.05/0.3)','(0.05/0.25)','(0.05/0.2)'...
    '(0.1/0.3)','(0.1/0.25)','(0.1/0.2)','(0.15/0.3)','(0.15/0.25)','(0.15/0.2)'});
set(gca,'fontsize',16,'Ycolor','k','Xcolor','k');
box off

delta_PI(:,1)=[0.001,0.001,0.001,0.05,0.05,0.05,0.1,0.1,0.1,0.15,0.15,0.15];
delta_PI(:,2)=[0.3,0.25,0.2,0.3,0.25,0.2,0.3,0.25,0.2,0.3,0.25,0.2];
% 
% mean_p_sim_w{w,count}=mean(p_sim);
% mean_p_sim_w{w,count}=mean(p_diff);
% 
% figure
% for i=[2,3,5,6,8,9,11,12]
%     hold on
%     if i==2
%         scatter(KL2(:,i),thr_99(:,i),50,'ko','filled')
%     else
%         scatter(KL2(:,i),thr_99(:,i),40,'o','filled')
%     end
% end
% h=refline(1,0);
% h.LineStyle='--';
% h.LineWidth=2.5;
% h.Color='k';
% ylabel('99th percentile null distribution');
% xlabel('KL divergence - data');
% set(gca,'fontsize',16);
% legend('\DeltaPI = 0.001/0.25','\DeltaPI = 0.001/0.2','\DeltaPI = 0.05/0.25','\DeltaPI = 0.05/0.2','\DeltaPI = 0.1/0.25',...
%     '\DeltaPI = 0.1/0.2','\DeltaPI = 0.15/0.25','\DeltaPI = 0.15/0.2')
% 
% % diffPI=delta_PI(:,2)-delta_PI(:,1);
% % diffPI([2,3,5,6,8,9,11,12])
% 
% samples_diff=sum(num_diff);
% samples_sim=sum(num_diff);
% samples_tot=samples_diff+samples_sim;
% samples_subset=samples_tot([2,3,5,6,8,9,11,12]);
% [~,index]=sort(samples_subset,'descend');
% 
% %color code
% cc=cool(length(index));
% 
% figure
% count=0;
% for i=[2,3,5,6,8,9,11,12]
%     count=count+1;
%     hold on
%     %     if i==2
%     %         scatter(KL2(:,i),thr_99(:,i),50,'ko','filled')
%     %     else
%     scatter(KL2(:,i),thr_99(:,i),40,'o','filled','MarkerFaceColor',cc(index(count),:));
%     %     end
% end
% h=refline(1,0);
% h.LineStyle='--';
% h.LineWidth=2.5;
% h.Color='k';
% ylabel('99th percentile null distribution');
% xlabel('KL divergence - data');
% set(gca,'fontsize',16);
% legend('\DeltaPI = 0.001/0.25','\DeltaPI = 0.001/0.2','\DeltaPI = 0.05/0.25','\DeltaPI = 0.05/0.2','\DeltaPI = 0.1/0.25',...
%     '\DeltaPI = 0.1/0.2','\DeltaPI = 0.15/0.25','\DeltaPI = 0.15/0.2')
% 
% 
% 
% 
% 
% find(KL2>thr_999)
% find(KL2>thr_99)
% find(KL2>thr_95)
% 
% figure
% scatter(KL2,thr_95,40,'o','filled')
% h=refline(1,0);
% h.LineStyle='--';
% h.LineWidth=2.5;
% h.Color='k';
% ylabel('95th percentile of null distribution');
% xlabel('KL divergence - data');
% set(gca,'fontsize',16)
% title('PI')
% 
% figure
% scatter(KL2,thr_99,40,'o','filled')
% h=refline(1,0);
% h.LineStyle='--';
% h.LineWidth=2.5;
% h.Color='k';
% ylabel('99th percentile of null distribution');
% xlabel('KL divergence - data');
% set(gca,'fontsize',16)
% title('PI')
% 
% [h,c]=histcounts(p_rank,0:0.05:1,'Normalization','cdf');
