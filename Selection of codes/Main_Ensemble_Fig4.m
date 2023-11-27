%Fig 4b,c
%xtended data Figure 10o,p

%% clear all
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
    ws_prob=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_',mice(m,:),'.mat']]); 
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
                    big_table(count,11)=ws_prob.WS_stat.wave_score_prob(day,s); %WS - Ent
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


%% Wave sessions

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

%% Calculates the transition matrices
N_sh=5; %in the paper we used 100
clus=10;
count=0;

for w=1:length(waves)
    row_w=waves(w);
    disp(w);    
    count=count+1;
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);
    munit=big_table(waves(w),5);    
    dt=floor(big_table(waves(w),8));
    if isinteger(dt)
    else
        dt=floor(dt);
    end
    
%     dt_Ent=floor(big_table(waves(w),11));

    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_spk,'-mat');
    spikes=full(spikes_d_s);
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);

    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);        
    
    spikes_w=[];
    for i=1:size(table_u,1)
        spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
    end
        
    %Calculates the probability of sequence length
    [prob_data(w,:)] = comp_wave_prob_f(spikes_w,floor(dt),clus);
    for i=1:N_sh
        mat_sh=shuffle(spikes_w')';
        [prob_data_sh(i,:)]=comp_wave_prob_f(mat_sh,floor(dt),clus);
        clear mat_sh
    end
    prob_data_sh_allsessions{count}=prob_data_sh;
    
    %TM
    [TM(:,:,count),adj_sh] = comp_TM_quick(spikes_w,floor(dt),clus,N_sh);
    TM_sh_allsessions{count}=adj_sh;
    
    
    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2
    
end


%% Threshold transition matrices

prctile_th=95;

for count=1:length(waves)
    
    adj_sh=TM_sh_allsessions{count};

    %Calculating the threshold
    subset_size=ceil(N_sh/2);
    subset2_sh=adj_sh(:,:,subset_size:N_sh);
    thr_mat=prctile(subset2_sh,prctile_th,3);
    
    aux_mat=zeros(clus,clus); %TM thresholded
    aux_mat_2=TM(:,:,count);  %TM without threshold
    aux_mat(aux_mat_2>thr_mat)=aux_mat_2(aux_mat_2>thr_mat);
    TM_th(:,:,count)=aux_mat;
    
    for sh=1:subset_size
        aux_mat_sh=zeros(clus,clus);
        aux2=adj_sh(:,:,sh);
        aux_mat_sh(aux2>thr_mat)=aux2(aux2>thr_mat);
        adj_sh_th(:,:,sh)=aux_mat_sh;
        clear aux_mat_sh aux2
    end
    
    TM_sh_th_allsessions{count}=adj_sh_th;
    
end
    
%% Plot prob wave length

% Pooling sessions of recorded data
% for i=1:length(waves)
mean_data=mean(prob_data,1);
%     mean_sh=mean(mean_prob_sh,1);
sem_data=std(prob_data,[],1)./sqrt(length(waves));
std_data=std(prob_data,[],1);
%     std_sh=std(mean_prob_sh,[],1);
% end

prob_sh=[];
for i=1:length(waves)
    prob_sh=[prob_sh;prob_data_sh_allsessions{i}];
    mean_prob_sh(i,:)=mean(prob_data_sh_allsessions{i});

end

mean_sh=mean(prob_sh,1);
sem_sh=std(prob_sh,[],1)./sqrt(size(prob_sh,1));
std_sh=std(prob_sh,[],1);

mean_meansh=mean(mean_prob_sh);
sd_meansh=std(mean_prob_sh);

% Figure with S.D
figure
errorbar(2:clus,mean_sh(2:end),std_sh(2:end),'linewidth',2,'Color',[0.8500, 0.3250, 0.0980]	);
hold on
errorbar(2:clus,mean_data(2:end),std_data(2:end),'linewidth',2,'Color',[0, 0.4470, 0.7410]);
axis([2 10 0 1]);
axis square
box off
xlabel('Number of ensembles');
ylabel({'Probability of';'sequential activation'});
set(gca,'fontsize',20, 'YColor','k', 'XColor','k')
legend('Data','Shuffle')
legend boxoff
axis()

for i=1:10
    [p_probr(i),h_probr(i),stats_probr{i}]=ranksum(prob_data(:,i)',prob_sh(:,i)','tail','right');
end

for i=1:10
    [p_probr(i),h_probr(i),stats_probr{i}]=ranksum(prob_data(:,i)',prob_sh(:,i)');
end

%% Quantification of transitions 10 ->1

for w=1:15
    for i=1:9
        trans(w,i)=TM(i,i+1,w);
    end
    trans(w,10)=TM(10,1,w);
end

figure
boxplot(trans)
ylabel('Probability');
xlabel('Transition');
xticks([1,2,3,4,5,6,7,8,9,10]);
xticklabels({'1->2','2->3','3->4','4->5','5->6','6->7','7->8','8->9','9->10','10->1'});
box off
set(gca,'fontsize',16,'xcolor','k','ycolor','k');
yticks([0 0.06 0.12])


[p,tbl,stat] = friedman(trans);

count=0;
for i=1:10
    for j=i+1:10
        count=count+1;
        [p(count)]=ranksum(trans(:,i),trans(:,j));
    end
end
nonsig=length(find(p>(0.05/45)));

%% Plot TM and graphs
% 
% count=length(waves);
% 
% figure
% for i=1:count
%     subplot(4,4,i)
%     imagesc(TM(:,:,i));
%     colormap hot
% %     colorbar
%     caxis([0 0.1])
%     set(gca,'fontsize',20)
%     xticks([]);
%     yticks([]);
%     axis square
% end
% 
% 
% 
% figure
% for i=1:count
%     subplot(4,4,i)
%     imagesc(TM_th(:,:,i));
%     colormap hot
% %     colorbar   
%     caxis([0 0.1])
%   set(gca,'fontsize',20)
%     xticks([]);
%     yticks([]);
%     axis square
% end
% close all
% % c=colorbar;
% % set(c,'XTick',[0 0.08])
% % set(gca,'fontsize',20)
% % ylabel('Ensemble #');
% % xlabel('Ensemble #');
% 
% figure
% for i=1:count
%     subplot(4,4,i)
%     imagesc(TM_sh_allsessions{1,i}(:,:,2));
%     colormap hot
% %     colorbar
%     caxis([0 0.1])
%     axis square
%     xticks([]);
%     yticks([]);
% end
% close all
% 
% figure
% for i=1:count
%     subplot(4,4,i)
%     imagesc(TM_sh_th_allsessions{1,i}(:,:,2));
%     colormap hot
% %     colorbar
%     caxis([0 0.08])
%     axis square
% end
% close all
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Graphs
% 
% figure
% for i=1:count
%     G{i}=digraph(TM_th(:,:,i));
%     subplot(4,4,i)
%     LWidths = 5* G{i}.Edges.Weight/max( G{i}.Edges.Weight);
%     p=plot(G{i},'LineWidth',LWidths,'Layout','circle');
%     p.EdgeColor= 'k';
%     p.MarkerSize = 12;
%     cc=parula(10);
%     for n=1:10
%         highlight(p,n,'NodeColor',cc(n,:))
%     end
%     axis off
%     axis square
%     box off
%     labelnode(p,1:10,'');
%     colormap(parula(10));    
% end
% close all
% 
% figure
% for i=1:count   
%     G{i}=digraph(TM(:,:,i));
%     subplot(4,4,i)
%     LWidths = 5* G{i}.Edges.Weight/max( G{i}.Edges.Weight);
%     p=plot(G{i},'LineWidth',LWidths,'Layout','circle');
%     p.EdgeColor= 'k';
%     p.MarkerSize = 12;
%     cc=parula(10);
%     for n=1:10
%         highlight(p,n,'NodeColor',cc(n,:))
%     end
%     axis off
%     axis square
%     box off
%     labelnode(p,1:10,'');
%     colormap(parula(10));
%     
% end
% close all
% 
% figure
% for i=1:count   
%     G{i}=digraph(TM_sh_allsessions{1,i}(:,:,2));
%     subplot(4,4,i)
%     LWidths = 5* G{i}.Edges.Weight/max( G{i}.Edges.Weight);
%     p=plot(G{i},'LineWidth',LWidths,'Layout','circle');
%     p.EdgeColor= 'k';
%     p.MarkerSize = 12;
%     cc=parula(10);
%     for n=1:10
%         highlight(p,n,'NodeColor',cc(n,:))
%     end
%     axis off
%     axis square
%     box off
%     labelnode(p,1:10,'');
%     colormap(parula(10));
%     
% end
% close all
% 
% figure
% for i=1:count   
%     G{i}=digraph(TM_sh_th_allsessions{1,i}(:,:,2));
%     subplot(4,4,i)
%     LWidths = 5* G{i}.Edges.Weight/max( G{i}.Edges.Weight);
%     p=plot(G{i},'LineWidth',LWidths,'Layout','circle');
%     p.EdgeColor= 'k';
%     p.MarkerSize = 12;
%     cc=parula(10);
%     for n=1:10
%         highlight(p,n,'NodeColor',cc(n,:))
%     end
%     axis off
%     axis square
%     box off
%     labelnode(p,1:10,'');
%     colormap(parula(10));
%     
% end
% close all


%% Figures for Example session

i=8;
figure
imagesc(TM(:,:,i));
colormap hot
%     colorbar
caxis([0 0.1])
set(gca,'fontsize',20)
xticks([]);
yticks([]);
axis square

figure
imagesc(TM_th(:,:,i));
colormap hot
%     colorbar
caxis([0 0.1])
set(gca,'fontsize',20)
xticks([]);
yticks([]);
axis square

figure
imagesc(TM_sh_allsessions{1,i}(:,:,2));
colormap hot
%     colorbar
caxis([0 0.1])
axis square
xticks([]);
yticks([]);

figure
imagesc(TM_sh_th_allsessions{1,i}(:,:,2));
colormap hot
%     colorbar
caxis([0 0.08])
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Graphs

figure
G{i}=digraph(TM_th(:,:,i));
LWidths = 5* G{i}.Edges.Weight/max( G{i}.Edges.Weight);
p=plot(G{i},'LineWidth',LWidths,'Layout','circle');
p.EdgeColor= 'k';
p.MarkerSize = 12;
cc=parula(10);
for n=1:10
    highlight(p,n,'NodeColor',cc(n,:))
end
axis off
axis square
box off
labelnode(p,1:10,'');
colormap(parula(10));

figure
G{i}=digraph(TM(:,:,i));
LWidths = 5* G{i}.Edges.Weight/max( G{i}.Edges.Weight);
p=plot(G{i},'LineWidth',LWidths,'Layout','circle');
p.EdgeColor= 'k';
p.MarkerSize = 12;
cc=parula(10);
for n=1:10
    highlight(p,n,'NodeColor',cc(n,:))
end
axis off
axis square
box off
labelnode(p,1:10,'');
colormap(parula(10));

figure
G{i}=digraph(TM_sh_allsessions{1,i}(:,:,2));
LWidths = 5* G{i}.Edges.Weight/max( G{i}.Edges.Weight);
p=plot(G{i},'LineWidth',LWidths,'Layout','circle');
p.EdgeColor= 'k';
p.MarkerSize = 12;
cc=parula(10);
for n=1:10
    highlight(p,n,'NodeColor',cc(n,:))
end
axis off
axis square
box off
labelnode(p,1:10,'');
colormap(parula(10));

figure
G{i}=digraph(TM_sh_th_allsessions{1,i}(:,:,2));
LWidths = 5* G{i}.Edges.Weight/max( G{i}.Edges.Weight);
p=plot(G{i},'LineWidth',LWidths,'Layout','circle');
p.EdgeColor= 'k';
p.MarkerSize = 12;
cc=parula(10);
for n=1:10
    highlight(p,n,'NodeColor',cc(n,:))
end
axis off
axis square
box off
labelnode(p,1:10,'');
colormap(parula(10));
