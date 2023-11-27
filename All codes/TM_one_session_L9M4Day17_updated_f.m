%% Calculates transition probabilities.
clear all

rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
WS_path='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

mouse_name='L09M4';
day=17;
s=1; % number of session out of dates.sesnum(day) serssions

if mouse_name(2)=='0'
    mouse=[mouse_name(1),mouse_name(3:5)];
else
    mouse=mouse_name;
end

N_sh=1000;
clus=10;

load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
munit=dates.ses{day}(s);

file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load([WS_path ['WS_Osc_15_sf7p73II_',mouse_name,'.mat']]);

% dt
dt=WS_stat.dt(day,s);
if isinteger(dt)
else
    dt=floor(dt);
end

% Calculate sorting
load(file_name_spk,'-mat');
spikes=full(spikes_d_s);
[~,sorting,~]=get_sorting_smoothed(spikes,dt);

% Identification of waves
num_clus_discr=10;
make_fig=1;
[table_u,N,T]=identify_waves_latestversion_6(mouse,day,num_clus_discr,dt,make_fig,spikes);


%Prepares new spike matrix by keeping frames with waves only
spikes_w=[];
for i=1:size(table_u,1)
    spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
end

%Calculates the probability of wave length

if dt>96
    dt=96;
end

[prob_data,prob_data_amp] = comp_wave_prob(spikes_w,floor(dt),clus);

for i=1:N_sh
    mat_sh=shuffle(spikes_w')';
    [prob_data_sh(i,:),prob_data_amp_sh(i,:)]=comp_wave_prob(mat_sh,dt,clus);
    clear mat_sh
end        

%Calculates the transition matrix
[TM,TM_sh] = comp_TM_quick(spikes_w,dt,clus,N_sh);

%Apply threshold to TM
prctile_th=80;
subset_size=ceil(N_sh/2);
subset2_sh=TM_sh(:,:,subset_size:N_sh);
thr_mat=prctile(subset2_sh,prctile_th,3);

aux_mat=zeros(clus,clus);
aux_mat(TM>thr_mat)=TM((TM>thr_mat));
TM_th=aux_mat;

for sh=1:subset_size
    aux_mat_sh=zeros(clus,clus);
    aux2=TM_sh(:,:,sh);
    aux_mat_sh(aux2>thr_mat)=aux2(aux2>thr_mat);
    TM_sh_th(:,:,sh)=aux_mat_sh;
    clear aux_mat_sh
end

%% Figures: Prob of wave length

mean_prob_sh=mean(prob_data_sh,1);
std_prob_sh=std(prob_data_sh,[],1);

% mean_prob_amp_sh=mean(prob_data_amp_sh,1);
% std_prob_amp_sh=std(prob_data_amp_sh,[],1);

figure
plot(2:clus,prob_data(2:end),'linewidth',2);
hold on
errorbar(2:clus,mean_prob_sh(2:end),std_prob_sh(2:end),'linewidth',2);
set(gca,'fontsize',20)
xlabel('Wave length (ensembles)');
ylabel('Probability');
box off
legend('Data',['Shuffle - Mean ',char(177),' SD'])
legend boxoff 
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\Prob wave length\Prob_wave_length_L9M4Day17_Smoothing_Adjusted_dt_Linear_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\Prob wave length\Prob_wave_length_L9M4Day17_Smoothing_Adjusted_dt_Linear_IDENTIFY6.svg');
close all

%Log scale
figure
plot(2:clus,prob_data(2:end),'linewidth',2);
hold on
errorbar(2:clus,mean_prob_sh(2:end),std_prob_sh(2:end),'linewidth',2);
set(gca,'fontsize',20)
xlabel('Wave length (ensembles)');
ylabel('Probability');
box off
legend('Data',['Shuffle - Mean ',char(177),' SD'])
legend boxoff 
set(gca, 'YScale', 'log')
axis([2 8 0.00001 0.85])
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\Prob wave length\Prob_wave_length_L9M4Day17_Smoothing_Adjusted_dt_Log_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\Prob wave length\Prob_wave_length_L9M4Day17_Smoothing_Adjusted_dt_Log_IDENTIFY6.svg');
close all

% This looks at the amplitude, but it's not informative at all since in the
% shuffled case the signal oscillates back and forth between 1 and 10
% figure
% plot(1:clus,prob_data_amp(1:end),'linewidth',2);
% hold on
% errorbar(2:clus,mean_prob_amp_sh(2:end),std_prob_amp_sh(2:end),'linewidth',2);
% set(gca,'fontsize',20)
% xlabel('Wave length (ensembles)');
% ylabel('Probability');
% box off
% legend('Data',['Shuffle - Mean ',char(177),' SD'])
% legend boxoff 

%% Figures: Graphs and TM


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sorted unfiltered
figure
imagesc(TM);
colormap hot
colorbar
caxis([0 0.1])
axis square
c=colorbar;
set(c,'XTick',[0 0.1])
set(gca,'fontsize',20)
ylabel('Ensemble #');
xlabel('Ensemble #');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\TM_Unfiltered_Smoothing_Adjusted_dt_L9M4Day17_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\TM_Unfiltered_Smoothing_Adjusted_dt_L9M4Day17_IDENTIFY6.svg');
close all

G=digraph(TM);
figure
LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
p=plot(G,'LineWidth',LWidths,'Layout','circle');
p.EdgeColor= 'k';
p.MarkerSize = 12;
cc=parula(10);
for i=1:10
    highlight(p,i,'NodeColor',cc(i,:))
end
axis off
axis square
box off
labelnode(p,1:10,'');
% title('Sorted Unfiltered')
colormap(parula(10));
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\Graph_Unfiltered_Smoothing_Adjusted_dt_L9M4Day17_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\Graph_Unfiltered_Smoothing_Adjusted_dt_L9M4Day17_IDENTIFY6.svg');
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sorted filtered
figure
imagesc(TM_th);
colormap hot
colorbar
caxis([0 0.1])
axis square
c=colorbar;
set(c,'XTick',[0 0.1])
set(gca,'fontsize',20)
ylabel('Ensemble #');
xlabel('Ensemble #');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\TM_Filtered_Smoothing_Adjusted_dt_L9M4Day17_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\TM_Filtered_Smoothing_Adjusted_dt_L9M4Day17_IDENTIFY6.svg');
close all

Gs=digraph(TM_th);
figure
LWidths = 5*Gs.Edges.Weight/max(G.Edges.Weight);
p=plot(Gs,'LineWidth',LWidths,'Layout','circle');
p.EdgeColor= 'k';
p.MarkerSize = 12;
cc=parula(10);
for i=1:10
    highlight(p,i,'NodeColor',cc(i,:))
end
axis off
axis square
box off
labelnode(p,1:10,'');
% title('Sorted Filtered')
colormap(parula(10));
% c=colorbar;
% c.Box='off';
% c.Ticks=[0 1];
% c.TickLabels={'1','10'};
% set(gca,'fontsize',16)
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\Graph_Filtered_Smoothing_Adjusted_dt_L9M4Day17_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\Graph_Filtered_Smoothing_Adjusted_dt_L9M4Day17_IDENTIFY6.svg');
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Random unfiltered
figure
imagesc(TM_sh(:,:,9));
colormap hot
colorbar
caxis([0 0.1])
axis square
c=colorbar;
set(c,'XTick',[0 0.1])
set(gca,'fontsize',20)
ylabel('Ensemble #');
xlabel('Ensemble #');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\TM_Unfiltered_Smoothing_Adjusted_dt_L9M4Day17_Shuffle_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\TM_Unfiltered_Smoothing_Adjusted_dt_L9M4Day17_Shuffle_IDENTIFY6.svg');
close all

Gsh=digraph(TM_sh(:,:,9));
figure
LWidths = 5*Gsh.Edges.Weight/max(G.Edges.Weight);
p=plot(Gsh,'LineWidth',LWidths,'Layout','circle');
p.EdgeColor= 'k';
p.MarkerSize = 12;
cc=parula(10);
for i=1:10
    highlight(p,i,'NodeColor',cc(i,:))
end
axis off
axis square
box off
labelnode(p,1:10,'');
% title('Random Unfiltered')
colormap(parula(10));
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\Graph_Unfiltered_Smoothing_Adjusted_dt_L9M4Day17_Shuffle_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\Graph_Unfiltered_Smoothing_Adjusted_dt_L9M4Day17_Shuffle_IDENTIFY6.svg');
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Random filtered
figure
imagesc(TM_sh_th(:,:,9));
colormap hot
colorbar
caxis([0 0.1])
axis square
c=colorbar;
set(c,'XTick',[0 0.1])
set(gca,'fontsize',20)
ylabel('Ensemble #');
xlabel('Ensemble #');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\TM_Filtered_Smoothing_Adjusted_dt_L9M4Day17_Shuffle_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\TM_Filtered_Smoothing_Adjusted_dt_L9M4Day17_Shuffle_IDENTIFY6.svg');
close all

Gs_sh=digraph(TM_sh_th(:,:,9));
figure
LWidths = 5*Gs_sh.Edges.Weight/max(G.Edges.Weight);
p=plot(Gs_sh,'LineWidth',LWidths,'Layout','circle');
p.EdgeColor= 'k';
p.MarkerSize = 12;
cc=parula(10);
for i=1:10
    highlight(p,i,'NodeColor',cc(i,:))
end
axis off
axis square
box off
labelnode(p,1:10,'');
% title('Random Filtered')
colormap(parula(10));
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\Graph_Filtered_Smoothing_Adjusted_dt_L9M4Day17_Shuffle_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\TM\Graph_Filtered_Smoothing_Adjusted_dt_L9M4Day17_Shuffle_IDENTIFY6.svg');
close all


%% Graph theory measures
% 
% % 1) Node degree
% alfa=TM;
% alfa(alfa>0)=1;
% L_s=1./(alfa);
% % [dist_S,B_S]=distance_wei(L_s);
% G=digraph(alfa);
% degree_out = outdegree(G)';
% degree_in = indegree(G)';
% clear G
% 
% for sh=1:size(TM_sh,3)    
%     alfa=TM_sh(:,:,sh);
%     alfa(alfa>0)=1;
%     G_sh=digraph(alfa);
%         
%     degree_out_sh(sh,:) = outdegree(G_sh);
%     degree_in_sh(sh,:) = indegree(G_sh);
%     
%     clear G_sh
% end
% degree_out_sh_full=reshape(degree_out_sh,1,clus*size(TM_sh,3));
% degree_in_sh_full=reshape(degree_in_sh,1,clus*size(TM_sh,3));
%     
% 
% h=histogram(degree_out(:),0:10,'Normalization','Probability');
% % h=histogram(degree_out(:),0:10);
% p=h.Values;
% 
% h=histogram(degree_out_sh_full(:),0:10,'Normalization','Probability');
% % h=histogram(degree_out_sh_full(:),0:10);
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
% yticks([0 0.25 0.5])
% xticks([1 5 10])
% axis([0 11 0 0.5])
% 
% 
% deg_mean=mean(degree_out_sh,1);
% deg_std=std(degree_out_sh,[],1);
% 
% figure
% plot(1:10,degree_out,'-o','Markersize',5,'MarkerfaceColor',[0, 0.4470, 0.7410],'Linewidth',1.5);
% hold on
% errorbar(1:10,deg_mean,deg_std,'Linewidth',1.5);
% axis([0.5 10.5 0 10])
% box off
% ylabel('Node degree');
% xlabel('Node #');
% set(gca,'fontsize',20)
% legend({'Data';'Shuffle'})
% legend boxoff

% 2) Path length
lengths=10;
L_s=1./(TM);
alfa=TM;
alfa(alfa>0)=1;
% [dist_S]=distance_bin(triu(alfa));
[dist_WS,B_S]=distance_wei((L_s));

co=0;
for j=1:10
    for l=j+5:10  %Paths of length >=5
        co=co+1;
        paths(co)=B_S(j,l)./(l-j+1);
        length_path(co)=(l-j+1);        
    end
end

% max_dist=max(B_S(1,10),B_S(10,1));
% clear th

for sh=1:size(TM_sh,3)
    
    th_sh=TM_sh(:,:,sh);
    L_s=1./th_sh;
    alfa=th_sh;
    alfa(alfa>0)=1;
%     [dist_S_sh(:,:,sh)]=distance_bin(triu(alfa));
    [dist_WS_sh(:,:,sh),B_S_sh(:,:,sh)]=distance_wei((L_s));
    
    co=0;
    for j=1:10
        for l=j+5:10
            co=co+1;
            paths_sh(co,sh)=B_S_sh(j,l,sh)./(l-j+1);
            length_path_sh(co,sh)=(l-j+1);
        end
    end
    
%     max_dist_sh(sh)=max(B_S_sh(1,10,sh),B_S_sh(10,1,sh));
    
end
clear dist_WS_sh dist_WS_sh B_S_sh


path_sh_re=reshape(paths_sh,1,15*N_sh);
% figure
% hold on
% errorbar([1,2],[mean(paths),mean(path_sh_re)],[std(paths),std(path_sh_re)]);


figure
b=bar([1,2],[mean(paths),mean(paths_sh(:))],'Linewidth',1.5);
b.FaceColor = 'flat';
b.CData(1,:) = [0, 0.4470, 0.7410]*1.2;
b.CData(2,:) = [0.8500, 0.3250, 0.0980]*1.2;
hold on
er=errorbar([1,2],[mean(paths),mean(paths_sh(:))],[std(paths),std(paths_sh(:))]);
er.LineStyle='none';
er.Color='k';
er.LineWidth=1.5;
xticklabels({'Data','Shuffle'})
ylabel({'Normalized';'shortest path length'});
set(gca,'fontsize',20)
axis([0.3 2.7 0 0.7])
box off
% xlabel({'First line';'Second line'})
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\Shortest path length\ShortestPhatLength_Smoothing_Adjusted_dt_L9M4Day17_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\Shortest path length\ShortestPhatLength_Smoothing_Adjusted_dt_L9M4Day17_IDENTIFY6.svg');

%Now with the SEM


samples_data=length(paths(:));
samples_shuffle=length(path_sh_re(:));

figure
b=bar([1,2],[mean(paths(:)),mean(path_sh_re(:))],'Linewidth',1.5);
b.FaceColor = 'flat';
b.CData(1,:) = [0, 0.4470, 0.7410]*1.2;
b.CData(2,:) = [0.8500, 0.3250, 0.0980]*1.2;
hold on
er=errorbar([1,2],[mean(paths(:)),mean(path_sh_re(:))],[std(paths(:))/sqrt(samples_data),std(path_sh_re(:))/sqrt(samples_shuffle)]);
er.LineStyle='none';
er.Color='k';
er.LineWidth=1.5;
xticklabels({'Data','Shuffle'})
ylabel({'Normalized';'shortest path length'});
set(gca,'fontsize',20)
axis([0.3 2.7 0 0.7])
box off
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\Shortest path length\ShortestPhatLength_Smoothing_Adjusted_dt_L9M4Day17_SEM_IDENTIFY6.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 3\To use\Shortest path length\ShortestPhatLength_Smoothing_Adjusted_dt_L9M4Day17_SEM_IDENTIFY6.svg');
close all 

