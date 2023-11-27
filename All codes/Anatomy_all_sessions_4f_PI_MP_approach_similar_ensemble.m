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




%% PI one session
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
file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation


[a,b]=sort(PI_locked,'ascend');
PI_sorted=PI_locked(b);

Nl=length(PI_sorted);
fraction=0.1;
clear group1 group2
group1=b(1:floor(fraction*Nl));
group2=b(end-floor(fraction*Nl)+1:end);

%%%%%%%%%%%%%%%%%

clear distg1 distg2 distg12
countg1=0;
for i=1:length(group1)
    for j=i+1:length(group1)
        countg1=countg1+1;
        distg1(countg1)=delta_tissue(group1(i),group1(j));
    end
end

countg2=0;
for i=1:length(group2)
    for j=i+1:length(group2)
        countg2=countg2+1;
        distg2(countg2)=delta_tissue(group2(i),group2(j));
    end
end

countg12=0;
for i=1:length(group1)
    for j=1:length(group2)
        countg12=countg12+1;
        distg12(countg12)=delta_tissue(group1(i),group2(j));
    end
end

[p,h,stat]=ranksum([distg1],distg12);

clear mat_PP_distance
mat_PP_distance=nan(length(distg12),2);
mat_PP_distance(1:length([distg1]),1)=[distg1];
% mat_PP_distance(1:length([distg2]),2)=[distg2];
mat_PP_distance(1:length(distg12),2)=distg12;

figure
boxplot(mat_PP_distance)
xticklabels({'Similar PI','Different PI'});
ylabel('Pairwise anatomical distance (um)');
ax=set(gca,'fontsize',16,'YColor','k','XColor','k');
ylim([0 800]);
box off
title('one session')


% clear mat_PP_distance
mat_PI=nan(length(group1),2);
mat_PI(1:length([group1]),1)=PI_locked(group1);
mat_PI(1:length(group2),2)=PI_locked(group2);

figure
boxplot(mat_PI)
xticklabels({'Group1','Group2'});
ylabel('PI');
ax=set(gca,'fontsize',16,'YColor','k','XColor','k');
ylim([0 1]);
box off
% title(num2str(p))


figure
histogram(PI_locked(group1),[0:1/30:1]);
xlabel('Participation index');
ylabel('Number of cells');
% title('Group 1');
% set(gca,'fontsize',16,'XColor','k','YColor','k');
hold on
% figure
histogram(PI_locked(group2),[0:1/30:1]);
% xlabel('Preferred phase (rad)');
% ylabel('Counts');
% title('Group 2');
set(gca,'fontsize',16,'XColor','k','YColor','k');
legend('Group 1', 'Group2');
xticks([0 0.25 0.5 0.75 1]);
% xticklabels({'0','0','\pi'});
box off
xlim([0 1])

%% All sessions - PI

load('locking_all_sessions.mat')
load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\locking_all_sessions.mat');

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
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation


    [a,b]=sort(PI_locked,'ascend');
    PI_sorted=PI_locked(b);

    Nl=length(PI_sorted);

    count=0;
    for fraction=[0.05,0.1,0.2,0.3,0.4,0.5]
        count=count+1;
        clear group1 group2
        group1=b(1:floor(fraction*Nl));
        group2=b(end-floor(fraction*Nl)+1:end);

        %%%%%%%%%%%%%%%%%

        clear distg1 distg2 distg12
        countg1=0;
        for i=1:length(group1)
            for j=i+1:length(group1)
                countg1=countg1+1;
                distg1(countg1)=delta_tissue(group1(i),group1(j));
            end
        end

        countg2=0;
        for i=1:length(group2)
            for j=i+1:length(group2)
                countg2=countg2+1;
                distg2(countg2)=delta_tissue(group2(i),group2(j));
            end
        end

        countg12=0;
        for i=1:length(group1)
            for j=1:length(group2)
                countg12=countg12+1;
                distg12(countg12)=delta_tissue(group1(i),group2(j));
            end
        end
    
         mean_dist1(w,count)=mean(distg1);
        mean_dist12(w,count)=mean(distg12);

        distg1_s{w,count}=distg1;
        distg12_s{w,count}=distg12;

        clear pairs_sim pairs_diff p_sim p_diff dist_sim dist_diff distg12 distg1 group1 group2 distg2
    end

    clear diff_phase_reshape diff_pos_reshape delta_phase_mean delta_tissue spikes_d not_locked locked mean_p mean_p_locked prob_join dist_origin_locked dist_origin
    clear mean_phase_locked_ordered a b c d diff_phase_ordered_reshape  delta_phase_mean_ordered x y ov Anat ...
        p_sim p_diff pairs_sim pairs_diff dist_sim dist_diff delta_phase_mean diff_phase_reshape diff_pos_reshape delta_tissue r_i ...
        locked not_locked mean_p mean_p_locked dist_origin...
        p_sim_or p_diff_or dist_sim_or dist_diff_or pairs_sim_or pairs_diff_or diff_phase_reshape_or mean_phase_locked_ordered ...
        a b c d max_x max_y max_xy mean_p_sh N spikes_d delta_tissue edges_x dist_sim_or dist_diff_or pairs_sim_or pairs_diff_or  ...
        delta_phase_mean_or r_i dist_origin_locked dist_origin p_sim_or p_diff_or b
end


figure
boxplot([mean_dist1(:,2),mean_dist12(:,2)]);
ylim([0 800])
xticklabels({'Similar PI','Different PI'});
ylabel('Pairwise anatomical distance (um)');
title('All sessions');
set(gca,'FontSize',16,'XColor','k','Ycolor','k');
box off

[p,h,stat]=ranksum(mean_dist1(:,2),mean_dist12(:,2));

for i=1:size(mean_dist1,2)
    [p(i),h(i)]=ranksum(mean_dist1(:,i),mean_dist12(:,i));
end

fraction=[0.05,0.1,0.2,0.3,0.4, 0.5];
figure
plot(fraction,p,'k-*','linewidth',2);
xlim([0 0.55])
hold on
yline(0.05,'--','LineWidth',2,'Color',[47,79,79]/255);
ylabel('p-value');
xlabel({'Fraction of neurons used';'to define the groups with similar PP'});
set(gca,'fontsize',16,'Ycolor','k','Xcolor','k');
box off

