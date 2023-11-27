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
    load([save_data_path ['WS_Osc_14_',mice(m,:),'.mat']]);
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

%% PaS sessions

PaS_ses1=find(big_table(:,10)<=3); %Sessions inthe in MEC
PaS_ses2=find(big_table(:,10)>0); %Sessions inthe in MEC
PaS_ses=intersect(PaS_ses1,PaS_ses2);
adults=find(big_table(:,3)>15); %Sessions inthe of adults
sessions=intersect(PaS_ses,adults);
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
    N1=big_table(sessions(repeated_days(w)),12);
    N2=big_table(sessions(repeated_days(w)+1),12);
    
    if (N1>N2)
        count=count+1;
        index_to_delete(count)=repeated_days(w)+1;
    elseif  (N1<N2)
        count=count+1;
        index_to_delete(count)=repeated_days(w);
    else
        print('Other case');
    end
end

sessions(index_to_delete)=[];
big_table_2=big_table(sessions,:);
PaS_sessions=sessions;

%% Calculates the transition matrices
N_sh=500;
clus=10;
count=0;

for w=1:length(PaS_sessions)
    row_w=PaS_sessions(w);
    disp(w)
    
    mouse=['L',num2str(big_table(PaS_sessions(w),1)),'M',num2str(big_table(PaS_sessions(w),2))];
    day=big_table(PaS_sessions(w),3);
    s=big_table(PaS_sessions(w),4);
    munit=big_table(PaS_sessions(w),5);
    
    dt=floor(big_table(PaS_sessions(w),8));
    if isinteger(dt)
    else
        dt=floor(dt);
    end
    
    if isnan(dt)
        disp('Error')
        continue
    else
        count=count+1;
        num_clus_discr=10;
        downsample_factor=1;
        load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
        if isfield(dates,'actual_day')==1
            day=find(dates.actual_day==day);
        end
        file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        load(file_name_spk,'-mat');
        spikes=full(spikes_d_s);
        
        % There is no identification of waves in this sessions, because there are
        % no waves
        
        % Doesn't isolate wave epoch, instead it uses the entire spike matrix
        
        spikes_w=spikes;
        
        %Calculates the probability of wave length
        
%                 if dt>96
% % % %                     dt=96;
%                 end
        dt=66;
        
        [prob_data(w,:)] = comp_wave_prob(spikes_w,floor(dt),clus);
        for i=1:N_sh
            mat_sh=shuffle(spikes_w')';
            [prob_data_sh(i,:)]=comp_wave_prob(mat_sh,floor(dt),clus);
            clear mat_sh
        end
        prob_data_sh_allsessions{count}=prob_data_sh;
        
        %TM
        [TM(:,:,count),adj_sh] = comp_TM_quick(spikes_w,floor(dt),clus,N_sh);
        TM_sh_allsessions{count}=adj_sh;
        
    end
    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2
    
end

%Recalculating threshold
% for w=1:length(med)
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

ind=sum(prob_data,2);
ind2=find(ind==0);
prob_data_copy=prob_data;
prob_data(ind2,:)=[];

%% Threshond transition matrices

prctile_th=95;

for count=1:length(PaS_sessions)
    adj_sh=TM_sh_allsessions{count};

    subset_size=ceil(N_sh/2);
    subset2_sh=adj_sh(:,:,subset_size:N_sh);
    thr_mat=prctile(subset2_sh,prctile_th,3);
    
    aux_mat=zeros(clus,clus);
    aux_mat_2=TM(:,:,count);
    aux_mat(aux_mat_2>thr_mat)=aux_mat_2(aux_mat_2>thr_mat);
    TM_th(:,:,count)=aux_mat;
    
    for sh=1:subset_size
        aux_mat_sh=zeros(clus,clus);
        aux2=adj_sh(:,:,sh);
        aux_mat_sh(aux2>thr_mat)=aux2(aux2>thr_mat);
        adj_sh_th(:,:,sh)=aux_mat_sh;
        clear aux_mat_sh
    end
    
    TM_sh_th_allsessions{count}=adj_sh_th;
    
end
%% Plot prob wave length
ses_med=count;

% Pooling sessions of recorded data
for i=1:length(PaS_sessions)
    mean_data=mean(prob_data,1);
%     mean_sh=mean(mean_prob_sh,1);
    sem_data=std(prob_data,[],1)./sqrt(length(PaS_sessions));
    std_data=std(prob_data,[],1);
end

prob_sh=[];
for i=1:length(PaS_sessions)
    prob_sh=[prob_sh;prob_data_sh_allsessions{i}];
    mean_prob_sh(i,:)=mean(prob_data_sh_allsessions{i});
end

mean_sh=mean(prob_sh,1);
sem_sh=std(prob_sh,[],1)./sqrt(size(prob_sh,1));
std_sh=std(prob_sh,[],1);

%S.D. log scale
figure
errorbar(2:clus,mean_sh(2:end),std_sh(2:end),'linewidth',2,'Color',[0.8500, 0.3250, 0.0980]	);
hold on
errorbar(2:clus,mean_data(2:end),std_data(2:end),'linewidth',2,'Color',[0, 0.4470, 0.7410]);
alpha 0.8
axis([2 10 0 1]);
axis square
box off
xlabel('Number of ensembles');
ylabel({'Probability of';'sequential activation'});
set(gca,'fontsize',16,'Ycolor','k','Xcolor','k')
set(gca,'Yscale','log')
axis([2 10 0.00001 1])
yticks([0.00001 1])
yticklabels({'0','1'})

for i=1:10
    [p_prob(i),h_prob(i),stats_prob{i}]=ranksum(prob_data(:,i)',prob_sh(:,i)','tail','right');
end

for i=1:10
    [p_prob(i),h_prob(i),stats_prob{i}]=ranksum(prob_data(:,i)',prob_sh(:,i)');
end

% for c=1:clus
%     variance(c,1)=var(prob_data(:,c));
%     variance(c,2)=var(prob_sh(:,c));
%     [h_prob_wave(c),p_prob_wave(c),ci,stat_prob_wave(c)]=ttest2(prob_data(:,c),prob_sh(:,c));
%     [h_prob_wave_uv(c),p_prob_wave_uv(c),ci_uv,stat_prob_wave_uv(c)]=ttest2(prob_data(:,c),prob_sh(:,c),'Vartype','unequal');
% end
% 
% %Mean distribution per session
% sig_99=zeros(length(PaS_sessions),clus);
% sig_999=zeros(length(PaS_sessions),clus);
% sig_95=zeros(length(PaS_sessions),clus);
% 
% for i=1:length(PaS_sessions)
% mean_prob_sh(i,:)=mean(prob_data_sh_allsessions{i},1);
% sem_prob_sh(i,:)=std(prob_data_sh_allsessions{i},[],1)./sqrt(N_sh);
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



% Figure with S.D with n=25 both for data and shuffle

% meansh=mean(mean_prob_sh);
% stdsh=std(mean_prob_sh);
% 
% figure
% errorbar(2:clus,meansh(2:end),stdsh(2:end),'linewidth',2,'Color',[0.8500, 0.3250, 0.0980]	);
% hold on
% errorbar(2:clus,mean_data(2:end),std_data(2:end),'linewidth',2,'Color',[0, 0.4470, 0.7410]);
% axis([2 10 0 1]);
% axis square
% box off
% xlabel('Number of ensembles');
% ylabel({'Probability of';'sequential activation'});
% set(gca,'fontsize',20, 'YColor','k', 'XColor','k')
% legend('Data','Shuffle')
% legend boxoff
% set(gca,'fontsize',16,'Ycolor','k','Xcolor','k')
% set(gca,'Yscale','log')
% axis([2 10 0.00001 1])
% yticks([0.00001 1])
% yticklabels({'0','1'})
% 
% % for i=1:10
% %     [p_prob(i),h_prob(i),stats_prob{i}]=ranksum(prob_data(:,i),mean_prob_sh(:,i),'tail','right');
% % end
% 
% for i=1:10
%     [p_prob(i),h_prob(i),stats_prob{i}]=ranksum(prob_data(:,i),mean_prob_sh(:,i));
% end
% 
% %S.D. linear scale
% figure
% errorbar(2:clus,mean_data(2:end),std_data(2:end),'linewidth',2);
% hold on
% errorbar(2:clus,mean_sh(2:end),std_sh(2:end),'linewidth',2);
% axis([2 10 0 1]);
% axis square
% box off
% xlabel('Number of ensembles');
% ylabel({'Probability of';'sequential activation'});
% set(gca,'fontsize',20)
% % set(gca,'Yscale','log')
% axis([2 10 0.00001 1])
% 
% 
% 
% %S.E.M.
% figure
% for i=1:length(PaS_sessions)
% subplot(6,5,i)
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

% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Prob of wave length\Prob_wave_length_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_OneCurve.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Prob of wave length\Prob_wave_length_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_OneCurve.svg');
% close all

%% Plot TM and graphs

figure
for i=1:count
    subplot(6,5,i)
    imagesc(TM(:,:,i));
    colormap hot
%     colorbar
    caxis([0 0.1])
    set(gca,'fontsize',20)
    xticks([]);
    yticks([]);
    axis square
end
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\TM_Unfiltered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\TM_Unfiltered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels.svg');
% close all

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
    subplot(6,5,i)
    imagesc(TM_th(:,:,i));
    colormap hot
%     colorbar   
    caxis([0 0.1])
  set(gca,'fontsize',20)
    xticks([]);
    yticks([]);
    axis square
end
% c=colorbar;
% set(c,'XTick',[0 0.08])
% set(gca,'fontsize',20)
% ylabel('Ensemble #');
% xlabel('Ensemble #');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\TM_Filtered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\TM_Filtered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels.svg');
% close all

figure
for i=1:count
    subplot(6,5,i)
    imagesc(TM_sh_allsessions{1,i}(:,:,18));
    colormap hot
%     colorbar
    caxis([0 0.1])
    axis square
    xticks([]);
    yticks([]);
end
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\TM_Unfiltered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels_Shuffle.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\TM_Unfiltered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels_Shuffle.svg');
% close all

figure
for i=1:count
    subplot(6,5,i)
    imagesc(TM_sh_th_allsessions{1,i}(:,:,18));
    colormap hot
%     colorbar
    caxis([0 0.08])
    axis square
end
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\TM_Filtered_Smoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels_Shuffle.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\TM_Filtered_Smoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels_Shuffle.svg');
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Graphs

figure
for i=1:count
    G{i}=digraph(TM_th(:,:,i));
    subplot(6,5,i)
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
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\Graph_Filtered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\Graph_Filtered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels.svg');
% close all

figure
for i=1:count   
    G{i}=digraph(TM(:,:,i));
    subplot(6,5,i)
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
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\Graph_Unfiltered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\Graph_Unfiltered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels.svg');
% close all

figure
for i=1:count   
    G{i}=digraph(TM_sh_allsessions{1,i}(:,:,18));
    subplot(6,5,i)
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
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\Graph_Unfiltered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels_Shuffle.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\Graph_Unfiltered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels_Shuffle.svg');
% close all

figure
for i=1:count   
    G{i}=digraph(TM_sh_th_allsessions{1,i}(:,:,19));
    subplot(6,5,i)
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
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\Graph_Filtered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels_Shuffle.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\TM\Graph_Filtered_NOSmoothing_NOAdjusted_dt_MEDMEC-AllSessions_Panels_Shuffle.svg');
% close all

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
path_length=mean(paths,2); %mean across paths for each session
path_length_sem=std(paths,[],2)/sqrt(size(paths,2)); %sem across paths for each session

%Mean pooling sessions
mean_path_length=mean(path_length); %mean across sessions
sem_path_length=std(path_length)/sqrt(length(PaS_sessions)); %sem across sessions

sig_99_spl=zeros(1,length(PaS_sessions));
sig_999_spl=zeros(1,length(PaS_sessions));
sig_95_spl=zeros(1,length(PaS_sessions));
for i=1:length(PaS_sessions)
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

%Mean of SH - Poolong sessions
path_length_sem_sh_pooled=std(path_length_sh)/sqrt(length(PaS_sessions)); %mean scross sessions - Shuffle


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
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\ShortestPathLength_PaS_all_sessions_dt96.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July-August 2021\Result 5\ShortestPathLength_PaS_all_sessions_dt96.svg');
close all 
[p,h,stats] = ranksum(path_length(:),path_length_sh(:)');
[p,h,stats] = ranksum(paths(:),paths_sh(:));

figure
histogram(path_length(:))
hold on
histogram(path_length_sh(:))


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
% spar=discretize(big_table(med,6),3);
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
% caxis([min(big_table(med,6)) max(big_table(med,6))])
% % set(c,'ticks',[0.2 0.3 0.4])
% ylabel('Mean degree');
% xticks([1 2])
% xticklabels({'Data';'Shuffle'})


% Path from 1 to 10

%lengths=[5,6,7,8,9,10];
% 
% for i=1:count
%     L_s=1./(TM(:,:,i));
%     alfa=TM(:,:,i);
%     alfa(alfa>0)=1;
%     [dist_WS(:,:,i),B_S(:,:,i)]=distance_wei((L_s));
%     
%     co=0;
%     for j=1:10
%         for l=j+5:10
%             co=co+1;
%             paths(i,co)=B_S(j,l,i)./(l-j+1);
%             length_path(i,co)=(l-j+1);            
%         end
%     end
% end
% 
% 
% for i=1:count
%     for sh=1:size(TM_sh_allsessions{1,i},3)
%         th_sh=TM_sh_allsessions{1,i}(:,:,sh)';
%         L_s=1./th_sh;
%         alfa=th_sh;
%         alfa(alfa>0)=1;
%         [dist_WS_sh(:,:,sh),B_S_sh(:,:,sh)]=distance_wei((L_s));
%         
%         co=0;
%         for j=1:10
%             for l=j+5:10
%                 co=co+1;
%                 paths_sh(i,co,sh)=B_S_sh(j,l,sh)./(l-j+1);
%                 length_path_sh(i,co,sh)=(l-j+1);
%             end
%         end
% %         max_dist_sh(i,sh)=max(B_S_sh(1,10,sh),B_S_sh(10,1,sh));
%     end
%     clear dist_WS_sh dist_WS_sh B_S_sh
% end
% 
% path_length=mean(paths');
% path_length_std=std(paths');
% 
% path_length_sh=mean(mean(paths_sh,3)');
% path_length_sh_std=std(mean(paths_sh,3)');
% 
% mean_tot=mean(paths(:));
% std_tot=std(paths(:));
% mean_tot_shuffle=mean(paths_sh(:));
% std_tot_shuffle=std(paths_sh(:));
% 
% path_sh_re=reshape(paths_sh,count,15*N_sh);
% 
% figure
% b=bar([1,2],[mean(paths(:)),mean(path_sh_re(:))],'Linewidth',1.5);
% b.FaceColor = 'flat';
% b.CData(1,:) = [0, 0.4470, 0.7410]*1.2;
% b.CData(2,:) = [0.8500, 0.3250, 0.0980]*1.2;
% hold on
% er=errorbar([1,2],[mean(paths(:)),mean(path_sh_re(:))],[std(paths(:)),std(path_sh_re(:))]);
% er.LineStyle='none';
% er.Color='k';
% er.LineWidth=1.5;
% xticklabels({'Data','Shuffle'})
% ylabel({'Normalized';'shortest path length'});
% set(gca,'fontsize',20)
% axis([0.3 2.7 0 0.7])
% box off
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Shortest path length\ShortestPathLength_NOSmoothing_NOAdjusted_dt_Pooled data_MEDMEC-AllSessions.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Shortest path length\ShortestPathLength_NOSmoothing_NOAdjusted_dt_Pooled data_MEDMEC-AllSessions.svg');
% close all 
% 
% figure
% b=bar([1,2],[mean(path_length),mean(path_length_sh)],'Linewidth',1.5);
% b.FaceColor = 'flat';
% b.CData(1,:) = [0, 0.4470, 0.7410]*1.2;
% b.CData(2,:) = [0.8500, 0.3250, 0.0980]*1.2;
% hold on
% er=errorbar([1,2],[mean(path_length),mean(path_length_sh)],[std(path_length),std(path_length_sh)]);
% er.LineStyle='none';
% er.Color='k';
% er.LineWidth=1.5;
% xticklabels({'Data','Shuffle'})
% ylabel({'Normalized';'shortest path length'});
% set(gca,'fontsize',20)
% axis([0.3 2.7 0 0.7])
% box off
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Shortest path length\ShortestPathLength_NOSmoothing_NOAdjusted_dt_Avg of Avg_MEDMEC-AllSessions.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Shortest path length\ShortestPathLength_NOSmoothing_NOAdjusted_dt_Avg of Avg_MEDMEC-AllSessions.svg');
% close all 
%     
% 
% 
% %Now with the SEM
% 
% 
% samples_data=length(paths(:));
% samples_shuffle=length(path_sh_re(:));
% 
% figure
% b=bar([1,2],[mean(paths(:)),mean(path_sh_re(:))],'Linewidth',1.5);
% b.FaceColor = 'flat';
% b.CData(1,:) = [0, 0.4470, 0.7410]*1.2;
% b.CData(2,:) = [0.8500, 0.3250, 0.0980]*1.2;
% hold on
% er=errorbar([1,2],[mean(paths(:)),mean(path_sh_re(:))],[std(paths(:))/sqrt(samples_data),std(path_sh_re(:))/sqrt(samples_shuffle)]);
% er.LineStyle='none';
% er.Color='k';
% er.LineWidth=1.5;
% xticklabels({'Data','Shuffle'})
% ylabel({'Normalized';'shortest path length'});
% set(gca,'fontsize',20)
% axis([0.3 2.7 0 0.7])
% box off
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Shortest path length\ShortestPathLength_NOSmoothing_NOAdjusted_dt_Pooled data_MEDMEC-AllSessions_SEM.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Shortest path length\ShortestPathLength_NOSmoothing_NOAdjusted_dt_Pooled data_MEDMEC-AllSessions_SEM.svg');
% close all 
% 
% 
% samples_data=length(path_length);
% samples_shuffle=length(path_length_sh);
% 
% figure
% b=bar([1,2],[mean(path_length),mean(path_length_sh)],'Linewidth',1.5);
% b.FaceColor = 'flat';
% b.CData(1,:) = [0, 0.4470, 0.7410]*1.2;
% b.CData(2,:) = [0.8500, 0.3250, 0.0980]*1.2;
% hold on
% er=errorbar([1,2],[mean(path_length),mean(path_length_sh)],[std(path_length)/sqrt(samples_data),std(path_length_sh)/sqrt(samples_shuffle)]);
% er.LineStyle='none';
% er.Color='k';
% er.LineWidth=1.5;
% xticklabels({'Data','Shuffle'})
% ylabel({'Normalized';'shortest path length'});
% set(gca,'fontsize',20)
% axis([0.3 2.7 0 0.7])
% box off
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Shortest path length\ShortestPathLength_NOSmoothing_NOAdjusted_dt_Avg of Avg_MEDMEC-AllSessions_SEM.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Shortest path length\ShortestPathLength_NOSmoothing_NOAdjusted_dt_Avg of Avg_MEDMEC-AllSessions_SEM.svg');
% close all 
% 



%% Figures for L8M4Day17

i=14;

figure
imagesc(TM(:,:,i));
colormap hot
%     colorbar
caxis([0 0.1])
set(gca,'fontsize',20)
xticks([]);
yticks([]);
axis square
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\TM_Unfiltered_L8M4Day17_dt66.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\TM_Unfiltered_L8M4Day17_dt66.svg');
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
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\TM_Filtered_L8M4Day17_dt66.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\TM_Filtered_L8M4Day17_dt66.svg');
close all
% c=colorbar;
% set(c,'XTick',[0 0.08])
% set(gca,'fontsize',20)
% ylabel('Ensemble #');
% xlabel('Ensemble #');

figure
imagesc(TM_sh_allsessions{1,i}(:,:,19));
colormap hot
%     colorbar
caxis([0 0.1])
axis square
xticks([]);
yticks([]);
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\TM_Unfiltered_L8M4Day17_SH_dt66.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\TM_Unfiltered_L8M4Day17_SH_dt66.svg');
close all

figure
imagesc(TM_sh_th_allsessions{1,i}(:,:,19));
colormap hot
%     colorbar
caxis([0 0.08])
axis square
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\TM_Filtered_L8M4Day17_SH_dt66.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\TM_Filtered_L8M4Day17_SH_dt66.svg');
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
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\Graph_Filtered_L8M4Day17_dt66.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\Graph_Filtered_L8M4Day17_dt66.svg');
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

saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\Graph_Unfiltered_L8M4Day17_dt66.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\Graph_Unfiltered_L8M4Day17_dt66.svg');
close all

figure
G{i}=digraph(TM_sh_allsessions{1,i}(:,:,19));
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

saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\Graph_Unfiltered_L8M4Day17_SH_dt66.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\Graph_Unfiltered_L8M4Day17_SH_dt66.svg');
close all

figure
G{i}=digraph(TM_sh_th_allsessions{1,i}(:,:,19));
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


saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\Graph_Filtered_L8M4Day17_SH_dt66.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures November\Result 5\Graph_Filtered_L8M4Day17_SH_dt66.svg');
close all