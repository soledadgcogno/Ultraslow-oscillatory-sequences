%% Sorting correlation and characterization of the sorting
clear all
close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';


mice_number=6;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];
% mice=['L8M1';'L8M2';'L8M3';'L8M4';'L9M1';'L9M4'];


m=6;

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

day=17;
s=1;
disp(s)
munit=dates.ses{day}(s);
file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];


disp(day)
load(file_name,'-mat');
spikes_d=full(spikes_d_s);
[N,T]=size(spikes_d);



%Xcorr
maxlag=10;
downsampling_factor=4;
sorting_corr=sorting_xcorr_f(spikes_d,maxlag,downsampling_factor,117); %check the smoothing kernel
% aux=find(sorting_corr==ini);
% sorting_corr_ini=circshift(sorting_corr,-(aux-1));
% aux_w=find(sorting_corr==ini_w);
% sorting_corr_ini_w=circshift(sorting_corr,-(aux_w-1));
figure
spy(spikes_d(sorting_corr,:),'k')
pbaspect([23,2,1])
% figure
% spy(spikes_d(sorting_corr_ini,:),'k')
% pbaspect([23,2,1])



%% Distribution of time lags



spikes=spikes_d;

% Lagged correlation 
maxlag=480;%480;%250 - 520;
% for i=1:N
%     FR(i,:)=full(fire_rate(spikes(i,:),117,'g')); %Smoothing 
% end
downsampling_factor=4;
sf=7.73;
dt=117;
% FRp = spikes_downsample(FR,N,downsampling_factor);
% FRp = spikes_downsample(spikes,N,downsampling_factor);

% for i=1:N
%     FR(i,:)=full(fire_rate(FRp(i,:),dt,'g')); %Smoothing 
% end



FRp = spikes_downsample(spikes,N,downsampling_factor);
for i=1:N
    FR(i,:)=full(fire_rate(FRp(i,:),dt,'g')); 
end

% FR=FRp;

tic
ini=481;
fin=721;%755
count=0;
for i=1:N
    for j=1:N
        [val,time]=(xcorr(FR(i,:),FR(j,:),maxlag)); %Check whether I need to zscore
%         ini=find(time==0);
%         fin=find(time==320);
%        val=zscore(val);
        [v,in]=max(val(ini:fin));       
%         [v,in]=max(val);        

%         if time(in)>=0
%             Corr_mat(i,j)=v;
%             Corr_mat(j,i)=v;
            lag_mat(i,j)=time(in+ini);
%             lag_mat(j,i)=-time(in+ini);            
%         else
%             Corr_mat(i,j)=v;
%             Corr_mat(j,i)=v;
%            lag_mat(i,j)=time(in+ini);
%             lag_mat(j,i)=-time(in+ini);
%         end        
        clear val time
    end
end
toc

clear distance_sorting time_lag
count=0;
for i=1
    for j=i+1:N
        count=count+1;
        
        first=sorting_corr(i);
        second=sorting_corr(j);
        
        distance_sorting(count)=j-i;
        time_lag(count)=(lag_mat(first,second));
        
    end
end

bins_size=(1/sf)*4;
p=[0.4940 0.1840 0.5560];
gr=[0.4660 0.6740 0.1880];
or=[[0.8500 0.3250 0.0980]];
% 
% figure
% scatter(distance_sorting,time_lag*bins_size,'filled');
% alpha 0.5
% yticks([0,62,124]);
% yticklabels([0,floor(120*bins_size),floor(240*bins_size)]);
% xlabel({'Distance between cells';'in the sorting (number of cells)'});
% ylabel({'Time lag that maximizes the correlation' ;'betwen the activity of cell pairs (s)'});
% hold on 
% scatter(distance_sorting(1),time_lag(1)*bins_size,50,or,'filled');
% scatter(distance_sorting(199),time_lag(199)*bins_size,50,p,'filled');
% scatter(distance_sorting(401),time_lag(401)*bins_size,50,gr,'filled');
% set(gca,'fontsize',14);
% [p,S] = polyfit(distance_sorting(100:400),time_lag(100:400)*bins_size,1);
% [f,delta] = polyval(p,distance_sorting,S);
% hold on
% plot(distance_sorting,f,'--m','linewidth',2);
% axis([0 480 0 124])
% hold on
% % plot(distance_sorting,f+2*delta,'m--',distance_sorting,f-2*delta,'m--')
% 

figure
scatter(distance_sorting,time_lag*bins_size,'filled');
alpha 0.5
yticks([0,62,124]);
yticklabels([0,floor(120*bins_size),floor(240*bins_size)]);
xlabel({'Distance between cells';'in the sorting (number of cells)'});
ylabel({'Time lag that maximizes the correlation' ;'betwen the activity of cell pairs (s)'});
hold on 
% scatter(distance_sorting(1),time_lag(1)*bins_size,50,or,'filled');
% scatter(distance_sorting(199),time_lag(199)*bins_size,50,p,'filled');
% scatter(distance_sorting(401),time_lag(401)*bins_size,50,gr,'filled');
set(gca,'fontsize',14,'Ycolor','k','Xcolor','k');
[p,S] = polyfit(distance_sorting(100:400),time_lag(100:400)*bins_size,1);
[f,delta] = polyval(p,distance_sorting,S);
hold on
plot(distance_sorting,f,'--m','linewidth',2);
axis([0 480 0 124])
hold on

mdl = fitlm(distance_sorting(100:400),time_lag(100:400)*bins_size);
[rho,pval]=corr(distance_sorting',time_lag','type','Spearman');


% 
% figure
% scatter(distance_sorting,time_lag*bins_size,'filled');
% alpha 0.5
% yticks([0,62,124]);
% yticklabels([0,floor(120*bins_size),floor(240*bins_size)]);
% xlabel({'Distance between cells';'in the sorting (number of cells)'});
% ylabel({'Time lag that maximizes the correlation' ;'betwen the activity of cell pairs (s)'});
% hold on 
% % scatter(distance_sorting(1),time_lag(1)*bins_size,50,or,'filled');
% % scatter(distance_sorting(199),time_lag(199)*bins_size,50,p,'filled');
% % scatter(distance_sorting(401),time_lag(401)*bins_size,50,gr,'filled');
% set(gca,'fontsize',14);
% hold on
% plot(distance_sorting,23.661  + 0.15899*distance_sorting,'--m','linewidth',2);
% axis([0 480 0 124])
% hold on

% 3 example cross-correlations

p=[0.4940 0.1840 0.5560];
gr=[0.4660 0.6740 0.1880];
or=[[0.8500 0.3250 0.0980]];


bins_size=(1/sf)*4;
figure
%Example of xcorr 1
maxlag=480;%480;%250 - 520;
first=sorting_corr(1);
second=sorting_corr(31);
[val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
aux=zscore(val);
[v,in]=max(aux(ini:fin));
timelag_cell1=time(in+ini)*bins_size;
plot(time*bins_size,zscore(val),'Color',or,'linewidth',2);
hcell1=line([timelag_cell1,timelag_cell1],[0,v]);
hcell1.LineStyle=':';
hcell1.Color=or;
hcell1.LineWidth=2.5;
hcell1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
first=sorting_corr(1);
second=sorting_corr(200);
clear val
[val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
plot(time*bins_size,zscore(val),'Color',p,'linewidth',2);
aux=zscore(val);
[v,in]=max(aux(ini:fin));
timelag_cell2=time(in+ini)*bins_size;
hcell2=line([timelag_cell2,timelag_cell2],[0,v]);
hcell2.LineStyle=':';
hcell2.Color=p;
hcell2.LineWidth=2.5;
hcell2.Annotation.LegendInformation.IconDisplayStyle = 'off';
first=sorting_corr(1);
second=sorting_corr(402);
[val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
plot(time*bins_size,zscore(val),'Color',gr,'linewidth',2);
aux=zscore(val);
[v,in]=max(aux(ini:fin));
timelag_cell3=time(in+ini)*bins_size;
hcell3=line([timelag_cell3,timelag_cell3],[0,v]);
hcell3.LineStyle=':';
hcell3.Color=gr;
hcell3.LineWidth=2.5;
box off
set(gca,'fontsize',16,'Ycolor','k','Xcolor','k')
legend boxoff
ylabel('zscored cross-correlation');
xlabel('Time lag (s)');
h1=refline(0,0);
h1.LineStyle='--';
h1.Color='k';
h2=line([0,0],[-2,3]);
h2.LineStyle='--';
h2.Color='k';
axis([-240 240 -2 3]);
legend('distance=30','distance=199','distance=401');


% Another plot of 3 example cross-correlations

p=[0.4940 0.1840 0.5560];
gr=[0.4660 0.6740 0.1880];
or=[[0.8500 0.3250 0.0980]];

for C=[6,10,12,14,17,20]
bins_size=(1/sf)*4;
figure
%Example of xcorr 1
maxlag=480;%480;%250 - 520;
first=sorting_corr(1);
second=sorting_corr(C);
[val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
aux=zscore(val);
[v,in]=max(aux(ini:fin));
timelag_cell1=time(in+ini)*bins_size;
plot(time*bins_size,zscore(val),'Color',or,'linewidth',2);
hcell1=line([timelag_cell1,timelag_cell1],[0,v]);
hcell1.LineStyle=':';
hcell1.Color=or;
hcell1.LineWidth=2.5;
hcell1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
first=sorting_corr(1);
second=sorting_corr(200);
clear val
[val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
plot(time*bins_size,zscore(val),'Color',p,'linewidth',2);
aux=zscore(val);
[v,in]=max(aux(ini:fin));
timelag_cell2=time(in+ini)*bins_size;
hcell2=line([timelag_cell2,timelag_cell2],[0,v]);
hcell2.LineStyle=':';
hcell2.Color=p;
hcell2.LineWidth=2.5;
hcell2.Annotation.LegendInformation.IconDisplayStyle = 'off';
first=sorting_corr(1);
second=sorting_corr(402);
[val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
plot(time*bins_size,zscore(val),'Color',gr,'linewidth',2);
aux=zscore(val);
[v,in]=max(aux(ini:fin));
timelag_cell3=time(in+ini)*bins_size;
hcell3=line([timelag_cell3,timelag_cell3],[0,v]);
hcell3.LineStyle=':';
hcell3.Color=gr;
hcell3.LineWidth=2.5;
box off
set(gca,'fontsize',16,'Ycolor','k','Xcolor','k')
legend boxoff
ylabel('zscored cross-correlation');
xlabel('Time lag (s)');
h1=refline(0,0);
h1.LineStyle='--';
h1.Color='k';
h2=line([0,0],[-2,3]);
h2.LineStyle='--';
h2.Color='k';
axis([-240 240 -2 3]);
legend('distance=30','distance=199','distance=401');
title(C)
end

% Another plot of 3 example cross-correlations

p=[0.4940 0.1840 0.5560];
gr=[0.4660 0.6740 0.1880];
or=[[0.8500 0.3250 0.0980]];

for C=[6]
bins_size=(1/sf)*4;
figure
%Example of xcorr 1
maxlag=480;%480;%250 - 520;
first=sorting_corr(1);
second=sorting_corr(C);
[val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
aux=zscore(val);
[v,in]=max(aux(ini:fin));
timelag_cell1=time(in+ini)*bins_size;
plot(time*bins_size,zscore(val),'Color',or,'linewidth',2);
hcell1=line([timelag_cell1,timelag_cell1],[0,v]);
hcell1.LineStyle=':';
hcell1.Color=or;
hcell1.LineWidth=2.5;
hcell1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
first=sorting_corr(1);
second=sorting_corr(200);
clear val
[val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
plot(time*bins_size,zscore(val),'Color',p,'linewidth',2);
aux=zscore(val);
[v,in]=max(aux(ini:fin));
timelag_cell2=time(in+ini)*bins_size;
hcell2=line([timelag_cell2,timelag_cell2],[0,v]);
hcell2.LineStyle=':';
hcell2.Color=p;
hcell2.LineWidth=2.5;
hcell2.Annotation.LegendInformation.IconDisplayStyle = 'off';
first=sorting_corr(1);
second=sorting_corr(402);
[val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
plot(time*bins_size,zscore(val),'Color',gr,'linewidth',2);
aux=zscore(val);
[v,in]=max(aux(ini:fin));
timelag_cell3=time(in+ini)*bins_size;
hcell3=line([timelag_cell3,timelag_cell3],[0,v]);
hcell3.LineStyle=':';
hcell3.Color=gr;
hcell3.LineWidth=2.5;
box off
set(gca,'fontsize',16,'Ycolor','k','Xcolor','k')
legend boxoff
ylabel('zscored cross-correlation');
xlabel('Time lag (s)');
h1=refline(0,0);
h1.LineStyle='--';
h1.Color='k';
h2=line([0,0],[-2,3]);
h2.LineStyle='--';
h2.Color='k';
axis([-240 240 -2 3]);
legend('distance=5','distance=199','distance=401');
% title(C)
end

%Example of xcorr 1
maxlag=480;%480;%250 - 520;
first=sorting_corr(1);
second=sorting_corr(2);
 [val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
figure
plot(time/2,val,'k','linewidth',2);

%Example of xcorr 2
first=sorting_corr(1);
second=sorting_corr(200);
 [val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
figure
plot(time/2,val,'k','linewidth',2);

%Example of xcorr 3
first=sorting_corr(1);
second=sorting_corr(402);
 [val,time]=(xcorr(FR(first,:),FR(second,:),maxlag)); %Check whether I need to zscore
figure
plot(time/2,val,'k','linewidth',2);

%% Plot raster plot

%1 wave  17792       19171
% another wave   14375       15689


[coeff,score]=pca(spikes_d');
[sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
% [~,sorting_w,~]=get_sorting_smoothed(spikes_d,117);
sorting_pca=sorting_descend;
figure
spy(spikes_d(sorting_pca,:),'k')
pbaspect([23,2,1])
ini=sorting_pca(1);
aux=find(sorting_corr==ini);
sorting_corr_ini=circshift(sorting_corr,-(aux-1));

% figure
% spy(spikes_d(sorting_corr_ini,:),'k')
% pbaspect([23,2,1])

%corr
sorting_corr_ini=flip(sorting_corr_ini);
sorting_corr=flip(sorting_corr);

fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_corr(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
set(gcf, 'Renderer', 'opengl');
hold on
for i=1:size(spikes_d,1)
    scatter([14375:15689]./8,i*spikes_d((sorting_corr(i)),[14375:15689]),5,[0.6350 0.0780 0.1840],'filled')
    alpha 0.3
end
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - Corr sorting.fig');
% saveas(gcf,'C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\Result 5\Rasterplots with Dim Red and Correlation\Rasterplot L8M4Day17 - Corr sorting.svg');
% close all

% 24352       25287
% 25376       26424
% 26425       27391



fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spikes_d,2))./8,i*spikes_d((sorting_corr(i)),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([200 400])
set(gcf, 'Renderer', 'opengl');
hold on
for i=1:size(spikes_d,1)
    scatter([24352:25287]./8,i*spikes_d((sorting_corr(i)),[24352:25287]),5,[0.6350 0.0780 0.1840],'filled')
    alpha 0.3
end