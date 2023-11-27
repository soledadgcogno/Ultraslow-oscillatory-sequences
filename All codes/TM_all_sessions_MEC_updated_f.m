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
N_sh=100;
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
        
    %Calculates the probability of wave length
    
%     if dt>96
%         dt=96;
%     end
            
    [prob_data(w,:)] = comp_wave_prob_version2_f(spikes_w,floor(dt),clus);
    for i=1:N_sh
        mat_sh=shuffle(spikes_w')';
        [prob_data_sh(i,:)]=comp_wave_prob_version2_f(mat_sh,floor(dt),clus);
        clear mat_sh
    end
    prob_data_sh_allsessions{count}=prob_data_sh;
    
    %TM
%     [TM(:,:,count),adj_sh] = comp_TM_quick(spikes_w,floor(dt),clus,N_sh);
%     TM_sh_allsessions{count}=adj_sh;
    
    
    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2
    
end

%Recalculating threshold
% for w=1:length(waves)
%     prctile_th=70;
%     adj_sh=TM_sh_allsessions{w};
%     subset_size=ceil(N_sh/2);
%     subset2_sh=adj_sh(:,:,subset_size:N_sh);
%     thr_mat=prctile(subset2_sh,prctile_th,3);
%     
%     aux_mat=zeros(clus,clus);
%     aux_mat_2=TM(:,:,count);
%     aux_mat(aux_mat_2>thr_mat)=aux_mat_2(aux_mat_2>thr_mat);
%     TM_th(:,:,count)=aux_mat;
%     
%     for sh=1:subset_size
%         aux_mat_sh=zeros(clus,clus);
%         aux2=adj_sh(:,:,sh);
%         aux_mat_sh(aux2>thr_mat)=aux2(aux2>thr_mat);
%         adj_sh_th(:,:,sh)=aux_mat_sh;
%         clear aux_mat_sh
%     end
%     TM_sh_th_allsessions{w}=adj_sh_th;
% 
% end

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

for i=1:10
    [p_probr(i),h_probr(i),stats_probr{i}]=ranksum(prob_data(:,i)',prob_sh(:,i)','tail','right');
end

for i=1:10
    [p_probr(i),h_probr(i),stats_probr{i}]=ranksum(prob_data(:,i)',prob_sh(:,i)');
end


% for c=1:clus
%     variance(c,1)=var(prob_data(:,c));
%     variance(c,2)=var(prob_sh(:,c));
%     [h_prob_wave(c),p_prob_wave(c),ci,stat_prob_wave(c)]=ttest2(prob_data(:,c),prob_sh(:,c));
%     [h_prob_wave_uv(c),p_prob_wave_uv(c),ci_uv,stat_prob_wave_uv(c)]=ttest2(prob_data(:,c),prob_sh(:,c),'Vartype','unequal');
% end
% 
% %Mean distribution per session
% sig_99=zeros(length(waves),clus);
% sig_999=zeros(length(waves),clus);
% sig_95=zeros(length(waves),clus);
% 
% for i=1:length(waves)
% mean_prob_sh(i,:)=mean(prob_data_sh_allsessions{i},1);
% sem_prob_sh(i,:)=std(prob_data_sh_allsessions{i},[],1)./sqrt(N_sh);
% sd_prob_sh(i,:)=std(prob_data_sh_allsessions{i},[],1);
% 
% 
%     for c=1:clus
%         th99(i,c)=prctile(prob_data_sh_allsessions{i}(:,c),99);
%         th999(i,c)=prctile(prob_data_sh_allsessions{i}(:,c),99.9);
%         th95(i,c)=prctile(prob_data_sh_allsessions{i}(:,c),95);
%         
%         if (prob_data(i,c)>th99(i,c)) sig_99(i,c)=1; end
%         if (prob_data(i,c)>th999(i,c)) sig_999(i,c)=1; end
%         if (prob_data(i,c)>th95(i,c)) sig_95(i,c)=1; end
%     end   
% end

% % % % % Figure with S.E.M
% % % % figure
% % % % errorbar(2:clus,mean_sh(2:end),sem_sh(2:end),'linewidth',2);
% % % % hold on
% % % % errorbar(2:clus,mean_data(2:end),sem_data(2:end),'linewidth',2);
% % % % axis([2 10 0 1]);
% % % % axis square
% % % % box off
% % % % xlabel('Number of ensembles');
% % % % ylabel({'Probability of';'sequential activation'});
% % % % set(gca,'fontsize',20)
% % % % legend('Data','Shuffle')
% % % % legend boxoff
% % % % 





% % % % % % % % Figure with S.D with n=15 both for data and shuffle
% % % % % % % 
% % % % % % % meansh=mean(mean_prob_sh);
% % % % % % % stdsh=std(mean_prob_sh);
% % % % % % % 
% % % % % % % figure
% % % % % % % errorbar(2:clus,meansh(2:end),stdsh(2:end),'linewidth',2,'Color',[0.8500, 0.3250, 0.0980]	);
% % % % % % % hold on
% % % % % % % errorbar(2:clus,mean_data(2:end),std_data(2:end),'linewidth',2,'Color',[0, 0.4470, 0.7410]);
% % % % % % % axis([2 10 0 1]);
% % % % % % % axis square
% % % % % % % box off
% % % % % % % xlabel('Number of ensembles');
% % % % % % % ylabel({'Probability';'of sequential activation'});
% % % % % % % set(gca,'fontsize',20, 'YColor','k', 'XColor','k')
% % % % % % % legend('Data','Shuffle')
% % % % % % % legend boxoff
% % % % % % % 
% % % % % % % for i=1:10
% % % % % % %     [p_prob(i),h_prob(i),stats_prob{i}]=ranksum(prob_data(:,i),mean_prob_sh(:,i));
% % % % % % % end

% for i=1:10
%     [p_prob(i),h_prob(i),stats_prob{i}]=ranksum(prob_data(:,i)',mean_prob_sh(:,i)','tail','right');
% end
% 
% for i=1:10
%     [p_prob(i),h_prob(i),stats_prob{i}]=ranksum(prob_data(:,i)',mean_prob_sh(:,i)','tail','left');
% end


%Figure for each session
% figure
% for i=1:length(waves)
% subplot(4,4,i)
% plot(2:clus,prob_data(i,2:end),'linewidth',2);
% hold on
% errorbar(2:clus,mean_prob_sh(i,2:end),sem_prob_sh(i,2:end),'linewidth',2);
% set(gca,'fontsize',20)
% xlabel('Wave length (ensembles)');
% ylabel('Probability');
% box off
% % legend('Data',['Shuffle - Mean ',char(177),' SD'])
% % legend boxoff 
% xticks([2 6 10])
% axis square
% axis([2 10 0 1]);
% end




% 
% cc=parula(13);
% figure
% for i=1:length(waves)
% hold on
% plot(2:clus,prob_data(i,2:end),'linewidth',2);
% errorbar(2:clus,mean_prob_sh(i,2:end),std_prob_sh(i,2:end),'--','linewidth',2);
% set(gca,'fontsize',20)
% xlabel('Wave length (ensembles)');
% ylabel('Probability');
% box off
% % legend('Data',['Shuffle - Mean ',char(177),' SD'])
% % legend boxoff 
% axis square
% axis([2 10 0 1]);
% end



% figure
% for i=1:length(waves)
% hold on
% plot(2:clus,prob_data(i,2:end),'linewidth',2);
% errorbar(2:clus,mean_prob_sh(i,2:end),std_prob_sh(i,2:end),'--','linewidth',2);
% set(gca,'fontsize',20)
% xlabel('Wave length (ensembles)');
% ylabel('Probability');
% box off
% % legend('Data',['Shuffle - Mean ',char(177),' SD'])
% % legend boxoff 
% axis square
% set(gca, 'YScale', 'log')
% axis([2 8 -inf inf]);
% end

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

count=length(waves);

figure
for i=1:count
    subplot(4,4,i)
    imagesc(TM(:,:,i));
    colormap hot
%     colorbar
    caxis([0 0.1])
    set(gca,'fontsize',20)
    xticks([]);
    yticks([]);
    axis square
end
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Unfiltered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Unfiltered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_IDENTIFY6.svg');
close all

% figure
% for i=1:count
%     subplot(4,4,i)
%     imagesc(TM(:,:,i));
%     colormap hot
%     colorbar
%     caxis([0 0.1])
%     axis square
% end


figure
for i=1:count
    subplot(4,4,i)
    imagesc(TM_th(:,:,i));
    colormap hot
%     colorbar   
    caxis([0 0.1])
  set(gca,'fontsize',20)
    xticks([]);
    yticks([]);
    axis square
end
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Filtered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Filtered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_IDENTIFY6.svg');
close all
% c=colorbar;
% set(c,'XTick',[0 0.08])
% set(gca,'fontsize',20)
% ylabel('Ensemble #');
% xlabel('Ensemble #');

figure
for i=1:count
    subplot(4,4,i)
    imagesc(TM_sh_allsessions{1,i}(:,:,18));
    colormap hot
%     colorbar
    caxis([0 0.1])
    axis square
    xticks([]);
    yticks([]);
end
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Unfiltered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_Shuffle_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Unfiltered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_Shuffle_IDENTIFY6.svg');
close all

figure
for i=1:count
    subplot(4,4,i)
    imagesc(TM_sh_th_allsessions{1,i}(:,:,18));
    colormap hot
%     colorbar
    caxis([0 0.08])
    axis square
end
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Filtered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_Shuffle_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Filtered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_Shuffle_IDENTIFY6.svg');
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Graphs

figure
for i=1:count
    G{i}=digraph(TM_th(:,:,i));
    subplot(4,4,i)
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
end
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Filtered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Filtered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_IDENTIFY6.svg');
close all

figure
for i=1:count   
    G{i}=digraph(TM(:,:,i));
    subplot(4,4,i)
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
    
end
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Unfiltered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Unfiltered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_IDENTIFY6.svg');
close all

figure
for i=1:count   
    G{i}=digraph(TM_sh_allsessions{1,i}(:,:,18));
    subplot(4,4,i)
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
    
end
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Filtered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_Shuffle_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Filtered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_Shuffle_IDENTIFY6.svg');
close all

figure
for i=1:count   
    G{i}=digraph(TM_sh_th_allsessions{1,i}(:,:,18));
    subplot(4,4,i)
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
    
end
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Unfiltered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_Shuffle_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Unfiltered_Smoothing_NOAdjusted_dt_AllSessions_13Panels_Shuffle_IDENTIFY6.svg');
close all

%% Graph theory measures 

% Path from 1 to 10
%lengths=[5,6,7,8,9,10];

for i=1:count
    L_s=1./(TM_th(:,:,i));
    alfa=TM_th(:,:,i);
    alfa(alfa>0)=1;
    edges_num(i)=sum(sum(alfa))/(10*10-10);
    [dist_WS(:,:,i),B_S(:,:,i)]=distance_wei((L_s));
    
    co=0;
    for j=1:10
        for l=j+5:10
            co=co+1;
            paths(i,co)=B_S(j,l,i)./(l-j+1);
            length_path(i,co)=(l-j+1);            
        end
    end
    
    %one more control
    count_sh=0;
    for sh=1:250
        TM_sh_c=TM_th(randperm(10),randperm(10),i);
        L_s_sh_c=1./(TM_sh_c);
        [dist_WS_sh_c,B_S_sh_c]=distance_wei((L_s_sh_c));

        co=0;
        for j=1:10
            for l=j+5:10
                count_sh=count_sh+1;
                co=co+1;
                paths_sh_c(i,co,sh)=B_S_sh_c(j,l)./(l-j+1);
                length_path_sh_c(i,co,sh)=(l-j+1);
                path_sh_c_all(count_sh)=B_S_sh_c(j,l)./(l-j+1);
            end
        end
        clear dist_WS_sh_c B_S_sh_c
    end
    
    mean_paths_sh_c(i)=mean(path_sh_c_all);
    
    
    clear L_s alfa path_sh_c_all
end

for i=1:count
    for sh=1:size(TM_sh_th_allsessions{1,i},3)
        th_sh=TM_sh_th_allsessions{1,i}(:,:,sh)';
        L_s=1./th_sh;
        alfa=th_sh;
        alfa(alfa>0)=1;
        [dist_WS_sh(:,:,sh),B_S_sh(:,:,sh)]=distance_wei((L_s));
        
        co=0;
        for j=1:10
            for l=j+5:10
                co=co+1;
                paths_sh(i,co,sh)=B_S_sh(j,l,sh)./(l-j+1);
                length_path_sh(i,co,sh)=(l-j+1);
            end
        end
%         max_dist_sh(i,sh)=max(B_S_sh(1,10,sh),B_S_sh(10,1,sh));
    end
    clear dist_WS_sh dist_WS_sh B_S_sh
end

%Mean for each sessions
path_length=mean(paths,2);
path_length_sem=std(paths,[],2)/sqrt(size(paths,2));

%Mean pooling sessions
mean_path_length=mean(path_length);
sem_path_length=std(path_length)/sqrt(length(waves));

sig_99_spl=zeros(1,length(waves));
sig_999_spl=zeros(1,length(waves));
sig_95_spl=zeros(1,length(waves));
for i=1:length(waves)
    conc_pathlength=[];
    for sh=1:N_sh/2
        conc_pathlength=[conc_pathlength,paths_sh(i,:,sh)];
        path_length_sh_null_dist(i,sh)=mean(paths_sh(i,:,sh));
    end
        path_length_sh(i)=mean(conc_pathlength);
        path_length_sem_sh(i)=std(conc_pathlength)/sqrt(length(conc_pathlength));
        thre95_spl(i)=prctile(path_length_sh_null_dist(i,:),95);
        thr99_spl(i)=prctile(path_length_sh_null_dist(i,:),99);
        thr999_spl(i)=prctile(path_length_sh_null_dist(i,:),99.9);
        
        if (path_length(i)>thr99_spl(i)) sig_99_spl(i)=1; end
        if (path_length(i)>thr999_spl(i)) sig_999_spl(i)=1; end
        if (path_length(i)>thre95_spl(i)) sig_95_spl(i)=1; end
end

path_length_sem_sh_pooled=std(path_length_sh)/sqrt(length(waves));

figure
b=bar([1,2],[mean(path_length(:)),mean(path_length_sh(:))],'Linewidth',1.5);
b.FaceColor = 'flat';
b.CData(1,:) = [0, 0.4470, 0.7410]*1.2;
b.CData(2,:) = [0.8500, 0.3250, 0.0980]*1.2;
hold on
er=errorbar([1,2],[mean(path_length(:)),mean(path_length_sh(:))],[sem_path_length,path_length_sem_sh_pooled]);
er.LineStyle='none';
er.Color='k';
er.LineWidth=1.5;
xticklabels({'Data','Shuffle'})
ylabel({'Normalized';'shortest path length'});
set(gca,'fontsize',20)
axis([0.3 2.7 0 0.7])
box off

figure
b=bar([1,2,3],[mean(path_length(:)),mean(mean_paths_sh_c), mean(path_length_sh(:))],'Linewidth',1.5);
b.FaceColor = 'flat';
b.CData(1,:) = [0, 0.4470, 0.7410]*1.2;
b.CData(2,:) = [0.8500, 0.4470, 0.7410]*1.2;
b.CData(3,:) = [0.8500, 0.3250, 0.0980]*1.2;
hold on
er=errorbar([1,2,3],[mean(path_length(:)),mean(mean_paths_sh_c),mean(path_length_sh(:))],[sem_path_length,std(mean_paths_sh_c)/sqrt(15),path_length_sem_sh_pooled]);
er.LineStyle='none';
er.Color='k';
er.LineWidth=1.5;
xticklabels({'Data','Shuffle 1','Shuffle 2'})
ylabel({'Normalized';'shortest path length'});
set(gca,'fontsize',20)
axis([0.3 3.7 0 0.7])
box off

% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\ShortestPathLength_all_sessions.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\ShortestPathLength_all_sessions.svg');
% close all 
[p,h,stats] = ranksum(path_length(:),path_length_sh(:));
[p,h,stats] = ranksum(paths(:),paths_sh(:));

figure
histogram(path_length(:))
hold on
histogram(path_length_sh(:))

mean_MEC=path_length;

%% figure

mean_mec=[0.3156
    0.5000
    0.4208
    0.3716
    0.5000
    0.4570
    0.4727
    0.5000
    0.5141
    0.5461
    0.7010
    0.4347
    0.4132
    0.3152
    0.5000];

mean_v1=[   0.2249
    0.3260
    0.2736
    0.2952
    0.3137
    0.3067
    0.4351
    0.3090
    0.2836
    0.1888
    0.3601
    0.2479
    0.3014
    0.2514
    0.3186
    0.2178
    0.3248
    0.2425
    0.3053];
mean_pas=[0.2265
    0.2351
    0.3735
    0.2807
    0.2884
    0.2551
    0.2692
    0.2356
    0.2298
    0.2592
    0.1747
    0.2446
    0.3383
    0.2409
    0.2261
    0.2736
    0.2342
    0.2526
    0.1865
    0.2507
    0.3108
    0.2764
    0.2404
    0.2441
    0.2637
    0.2534
    0.2680
    0.2497
    0.2530];

figure
b2=bar([1,2,3],[mean(mean_MEC),mean(mean_v1),mean(mean_pas)],'Linewidth',1.5);
b.FaceColor = 'flat';
b.CData(1,:) = [0, 0.4470, 0.7410]*1.2;
b.CData(2,:) = [0.8500, 0.3250, 0.0980]*1.2;
hold on
er=errorbar([1,2,3],[mean(mean_MEC),mean(mean_v1),mean(mean_pas)],...
    [std(mean_MEC)/sqrt(length(mean_MEC)),std(mean_v1)/sqrt(length(mean_v1)),std(mean_pas)/sqrt(length(mean_pas))]);
er.LineStyle='none';
er.Color='k';
er.LineWidth=1.5;
xticks([1,2,3])
xticklabels({'MEC','V1','PaS'})
ylabel({'Normalized';'shortest path length'});
set(gca,'fontsize',20)
axis([0.3 3.7 0 0.7])
box off

[h12,p12,stats_mec_v1] = ranksum(mean_MEC(:),mean_v1(:));
[h13,p13,stats13_mec_pas] = ranksum(mean_MEC(:),mean_pas(:));
[h23,p23,stats23_v1_pas] = ranksum(mean_v1(:),mean_pas(:));



%% Figures for L9M4Day17

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
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_L9M4Day17.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_L9M4Day17.svg');
close all

% figure
% for i=1:count
%     subplot(4,4,i)
%     imagesc(TM(:,:,i));
%     colormap hot
%     colorbar
%     caxis([0 0.1])
%     axis square
% end


figure
imagesc(TM_th(:,:,i));
colormap hot
%     colorbar
caxis([0 0.1])
set(gca,'fontsize',20)
xticks([]);
yticks([]);
axis square
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Filtered_L9M4Day17.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Filtered_L9M4Day17.svg');
close all
% c=colorbar;
% set(c,'XTick',[0 0.08])
% set(gca,'fontsize',20)
% ylabel('Ensemble #');
% xlabel('Ensemble #');

figure
imagesc(TM_sh_allsessions{1,i}(:,:,18));
colormap hot
%     colorbar
caxis([0 0.1])
axis square
xticks([]);
yticks([]);
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Unfiltered_L9M4Day17_SH.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Unfiltered_L9M4Day17_SH.svg');
close all

figure
imagesc(TM_sh_th_allsessions{1,i}(:,:,18));
colormap hot
%     colorbar
caxis([0 0.08])
axis square
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Filtered_L9M4Day17_SH.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\TM_Filtered_L9M4Day17_SH.svg');
close all

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
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Filtered_L9M4Day17.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Filtered_L9M4Day17.svg');
close all

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

saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Unfiltered_L9M4Day17.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Unfiltered_L9M4Day17.svg');
close all

figure
G{i}=digraph(TM_sh_allsessions{1,i}(:,:,18));
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

saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Filtered_L9M4Day17_SH.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Filtered_L9M4Day17_SH.svg');
close all

figure
G{i}=digraph(TM_sh_th_allsessions{1,i}(:,:,18));
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


saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Unfiltered_L9M4Day17_SH.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 3\Graph_Unfiltered_L9M4Day17_SH.svg');
close all


%% Figures for L8M2Day19

i=2;

figure
imagesc(TM(:,:,i));
colormap hot
%     colorbar
caxis([0 0.1])
set(gca,'fontsize',20)
xticks([]);
yticks([]);
axis square
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\TM_L8M2Day19_MEC.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\TM_L8M2Day19_MEC.svg');
close all

% figure
% for i=1:count
%     subplot(4,4,i)
%     imagesc(TM(:,:,i));
%     colormap hot
%     colorbar
%     caxis([0 0.1])
%     axis square
% end


figure
imagesc(TM_th(:,:,i));
colormap hot
%     colorbar
caxis([0 0.1])
set(gca,'fontsize',20)
xticks([]);
yticks([]);
axis square
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\TM_Filtered_L8M2Day19_MEC.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\TM_Filtered_L8M2Day19_MEC.svg');
close all
% c=colorbar;
% set(c,'XTick',[0 0.08])
% set(gca,'fontsize',20)
% ylabel('Ensemble #');
% xlabel('Ensemble #');

figure
imagesc(TM_sh_allsessions{1,i}(:,:,18));
colormap hot
%     colorbar
caxis([0 0.1])
axis square
xticks([]);
yticks([]);
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\TM_Unfiltered_L8M2Day19_SH_MEC.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\TM_Unfiltered_L8M2Day19_SH_MEC.svg');
close all

figure
imagesc(TM_sh_th_allsessions{1,i}(:,:,18));
colormap hot
%     colorbar
caxis([0 0.08])
axis square
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\TM_Filtered_L8M2Day19_SH_MEC.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\TM_Filtered_L8M2Day19_SH_MEC.svg');
close all

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
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\Graph_Filtered_L8M2Day19_MEC.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\Graph_Filtered_L8M2Day19_MEC.svg');
close all

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

saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\Graph_Unfiltered_L8M2Day19_MEC.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\Graph_Unfiltered_L8M2Day19_MEC.svg');
close all

figure
G{i}=digraph(TM_sh_allsessions{1,i}(:,:,18));
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

saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\Graph_Unfiltered_L8M2Day19_SH_MEC.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\Graph_Unfiltered_L8M2Day19_SH_MEC.svg');
close all

figure
G{i}=digraph(TM_sh_th_allsessions{1,i}(:,:,18));
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


saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\Graph_Filtered_L8M2Day19_SH_MEC.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\Graph_Filtered_L8M2Day19_SH_MEC.svg');
close all
% 
% %Node degree
% for i=1:count
%     alfa=TM(:,:,i);
%     alfa(alfa>0)=1;
%     L_s=1./(alfa);
%     G{i}=digraph(alfa);
%     degree_out(i,:) = outdegree(G{i})';
%     degree_in(i,:) = indegree(G{i})';    
% end
% 
% for i=1:count
%     for sh=1:size(TM_sh{i},3)
%         
%         alfa=TM_sh{i}(:,:,sh);
%         alfa(alfa>0)=1;
%         L_s=1./(alfa);
%         G_sh=digraph(alfa);
%                
%         degree_out_sh(i,:,sh) = outdegree(G_sh);
%         degree_in_sh(i,:,sh) = indegree(G_sh);
%         
%         clear G_sh
%     end
%     degree_out_sh_full(i,:)=reshape(degree_out_sh(i,1:10,:),1,clus*size(TM_sh{i},3));
%     degree_in_sh_full(i,:)=reshape(degree_in_sh(i,1:10,:),1,clus*size(TM_sh{i},3));    
% end
% 
% 
% degree_sh_prctile_90=prctile(degree_out_sh_full,90,2);
% degree_sh_prctile_10=prctile(degree_out_sh_full,10,2);
% 
% 
% h=histogram(degree_out(:),0:10,'Normalization','Probability');
% p=h.Values;
% 
% h=histogram(degree_out_sh_full(:),0:10,'Normalization','Probability');
% p_sh=h.Values;
% 
% figure
% plot(1:10,p,'-*','linewidth',2);
% hold on
% plot(1:10,p_sh,'-*','linewidth',2);
% alpha 0.7
% ylabel('Probability');
% xlabel('Node degree (ensembles)');
% legend({'Data','Shuffled'});
% legend boxoff
% box off
% set(gca,'fontsize',20);
% yticks([0 0.25 0.4])
% xticks([1 5 10])
% axis([0 11 0 0.4])
% 
% mean_degree=mean(degree_out,2);
% std_degree=std(degree_out,[],2);
% cv_degree=mean_degree./std_degree;
% mean_degree_shuffle=mean(degree_out_sh_full,2);
% std_degree_shuffle=std(degree_out_sh_full,[],2);
% cv_degree_shuffle=mean_degree_shuffle./std_degree_shuffle;
% 
% spar=discretize(big_table(waves,6),3);
% figure
% hold on
% for i=1:count
% %     plot([1,2],[mean_degree(i),mean_degree_shuffle(i)],'*-','color',cc(spar(i),:),'linewidth',1.5);
%     plot([1,2],[mean_degree(i),mean_degree_shuffle(i)],'*-','color','k','linewidth',1.5);
% %     errorbar([1,2],[mean_degree(i),mean_degree_shuffle(i)],[0,std_degree_shuffle(i)],'*-','color','k','linewidth',1.5);
% 
% end
% axis([0.5 2.5 2 10])
% set(gca,'fontsize',20)
% colormap copper(3)
% % c=colorbar;
% caxis([min(big_table(waves,6)) max(big_table(waves,6))])
% % set(c,'ticks',[0.2 0.3 0.4])
% ylabel('Mean degree');
% xticks([1 2])
% xticklabels({'Data';'Shuffle'})