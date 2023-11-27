%L9M4

clear all
% close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Behavior L8\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';


% 9 seconds for L8M2 Day 19 ; 20 seconds for L9M4 Day 17 
mouse='L9M4';
day=17;
munit=0;

load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
file_name=[dpath_spikes ['spikes_120ms_Do_THR1_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name,'-mat');
spikes=full(spikes_d_s);
N=size(spikes,1);
sf=7.73;
sf_d=7.73/4;

[N,T]=size(spikes);

ensembles=N;
% FRp = spikes_downsample(spikes,N,8);
FRp = spikes_downsample(spikes,N,4);
FRp2 = spikes_downsample(spikes,N,1);

%LEM
for i=1:ensembles
    FR(i,:)=(full(fire_rate(FRp(i,:),5*29,'g')));  % 9 seconds for L8M2 Day 19 ; 15 seconds for L9M4 Day 17 
end
clear score
[coeff,score]=pca(FR');

addpath("C:\Users\xscogno\drtoolbox.tar");
[Y,X] = laplacian_eigen([score(:,1:5)],2,15,2);

% [Y,X] = laplacian_eigen([score(:,1),score(:,2),score(:,3),score(:,4),score(:,5)],2,15,2);
% [Y,X] = laplacian_eigen([score(:,1),score(:,2),score(:,3)],2,15,2);
figure
% plot3((1:size(Y,1))./8,-Y(1:end,1),-Y(1:end,2),'Linewidth',3)
plot3((1:size(Y,1)),-Y(1:end,1),-Y(1:end,2),'Linewidth',3)
% axis([0 2000 -inf inf -inf inf])
colormap viridis
cla
patch([1:size(Y,1) nan],[-Y(:,1)' nan],[-Y(:,2)' nan],[1:size(Y,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('Time (s)')
ylabel('LEM - Dimension 1');
zlabel('LEM - Dimension 2');
yticks([]);
zticks([]);
set(gca,'fontsize',16)
xticks([0 3865 7730])
xticklabels([0 3865/sf_d 7730/sf_d])
axis([0 8000 -inf inf -inf inf])

%LEM2
for i=1:ensembles
    FR(i,:)=(full(fire_rate(FRp2(i,:),4*floor(117),'g')));  % 9 seconds for L8M2 Day 19 ; 15 seconds for L9M4 Day 17 
end
clear score
[coeff,score]=pca(FR);

addpath("C:\Users\xscogno\drtoolbox.tar");
[Y,X] = laplacian_eigen([score(:,1:5)],2,15,2);

% [Y,X] = laplacian_eigen([score(:,1),score(:,2),score(:,3),score(:,4),score(:,5)],2,15,2);
% [Y,X] = laplacian_eigen([score(:,1),score(:,2),score(:,3)],2,15,2);
figure
% plot3((1:size(Y,1))./8,-Y(1:end,1),-Y(1:end,2),'Linewidth',3)
plot3((1:size(Y,1)),-Y(1:end,1),-Y(1:end,2),'Linewidth',3)
% axis([0 2000 -inf inf -inf inf])
colormap viridis
cla
patch([1:size(Y,1) nan],[-Y(:,1)' nan],[-Y(:,2)' nan],[1:size(Y,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('Time (s)')
ylabel('LEM - Dimension 1');
zlabel('LEM - Dimension 2');
yticks([]);
zticks([]);
set(gca,'fontsize',16)
xticks([0 3865 7730])
xticklabels([0 3865/sf_d 7730/sf_d])
axis([0 8000 -inf inf -inf inf])

rmpath("C:\Users\xscogno\drtoolbox.tar");

% clear score
% [coeff,score]=pca(FR');
figure
plot3((1:size(score,1)),score(1:end,1),score(1:end,2),'Linewidth',3)
axis([0 8000 -inf inf -inf inf])
% axis([0 2000 -inf inf -inf inf])
colormap viridis
cla
patch([1:size(score,1) nan],[score(:,1)' nan],[score(:,2)' nan],[1:size(score,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('Time (s)')
ylabel('PC1');
zlabel('PC2');
yticks([]);
zticks([]);
set(gca,'fontsize',16)
xticks([0 3865 7730])
xticklabels([0 3865/sf_d 7730/sf_d])
axis([0 8000 -inf inf -inf inf])

% c=colorbar( 'XTick', [400 1400]);
c=colorbar( 'XTick', [400 3000]);
% cticks([200 1400])
set(gca,'fontsize',16)



[Y,X] = laplacian_eigen([score(:,1:5)],2,15,2);
figure
plot(-Y(1:end,1),-Y(1:end,2),'Linewidth',3)
colormap viridis
cla
patch([-Y(:,1)' nan],[-Y(:,2)' nan],[1:size(Y,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('LEM - Dimension 1');
ylabel('LEM - Dimension 2');
yticks([]);
xticks([]);
set(gca,'fontsize',16)
box off
% xticks([0 3865 7730])
% xticklabels([0 3865/sf_d 7730/sf_d])
% axis([0 8000 -inf inf -inf inf])

%% LEM for one V1 session and one PaS session

clear all
close all

rec_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
spath = 'C:\Users\xscogno\Dropbox\SfN2019\Development\Rasterplots\L9M1\';

%PaS
mouse='L8M4';
load([rec_path,strcat('recording_dates_',mouse,'.mat')]);
day=17;
s=2;
munit=dates.ses{day}(s);
file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name,'-mat');
spikes=full(spikes_d_s);
[N,T]=size(spikes);
sf=7.73;
sf_d=7.73/4;
T_PaS=T;

% FRp = spikes_downsample(spikes,N,4);
FRp = spikes_downsample(spikes,N,1);

%LEM
for i=1:N
%     FR(i,:)=(full(fire_rate(FRp(i,:),29,'g')));  
    FR(i,:)=(full(fire_rate(FRp(i,:),4*floor(66/1),'g')));

end
clear score
[coeff,score]=pca(FR');


addpath("C:\Users\xscogno\MATLAB\drtoolbox.tar");
[Y,X] = laplacian_eigen([score(:,1:5)],2,15,2);

figure
plot3((1:size(Y,1)),-Y(1:end,1),-Y(1:end,2),'Linewidth',3)
colormap viridis
cla
patch([1:size(Y,1) nan],[-Y(:,1)' nan],[-Y(:,2)' nan],[1:size(Y,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('Time (s)')
ylabel('LEM - Dimension 1');
zlabel('LEM - Dimension 2');
yticks([]);
zticks([]);
set(gca,'fontsize',16)
xticks([0 1546 3092])
xticklabels([0 1546/sf_d 3092/sf_d])
axis([0 3500 -inf inf -inf inf])

%LEM projection PaS

figure
plot(-Y(1:end,1),-Y(1:end,2),'Linewidth',3)
colormap viridis
cla
patch([-Y(:,1)' nan],[-Y(:,2)' nan],[1:size(Y,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('LEM - Dimension 1');
ylabel('LEM - Dimension 2');
yticks([]);
xticks([]);
set(gca,'fontsize',16)
box off

rmpath("C:\Users\xscogno\MATLAB\drtoolbox.tar");

figure
plot3((1:size(score,1)),score(1:end,1),score(1:end,2),'Linewidth',3)
axis([0 8000 -inf inf -inf inf])
% axis([0 2000 -inf inf -inf inf])
colormap viridis
cla
patch([1:size(score,1) nan],[score(:,1)' nan],[score(:,2)' nan],[1:size(score,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('Time (s)')
ylabel('PC1');
zlabel('PC2');
yticks([]);
zticks([]);
set(gca,'fontsize',16)
xticks([0 1546 3092])
xticklabels([0 1546/sf_d 3092/sf_d])
axis([0 3500 -inf inf -inf inf])
caxis([0 3500]);
% cticks([200 1400])
set(gca,'fontsize',16)
cbh = colorbar ; %Create Colorbar
cbh.Ticks = [0 3092] ; %Create 8 ticks from zero to 1
cbh.TickLabels = [0/sf_d 3092/sf_d] ;
set(gca,'fontsize',16)

%PCA projection PaS
figure
plot(score(1:end,1),score(1:end,2),'Linewidth',3)
    alfa=max([abs(score(:,2));abs(score(:,1))]);
    axis([-alfa-0.05 alfa+0.05 -alfa-0.05 alfa+0.05]);
% axis([-40 40 -50 50])
% axis([0 2000 -inf inf -inf inf])
colormap viridis
cla
patch([score(:,1)' nan],[score(:,2)' nan],[1:size(score,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('PC1');
ylabel('PC2');
xticks([]);
yticks([]);
set(gca,'fontsize',16)
box off
% c=colorbar( 'XTick', [400 1400]);
cbh = colorbar ; %Create Colorbar
cbh.Ticks = [1 size(score(:,1),1)] ; %Create 8 ticks from zero to 1
% cbh.TickLabels = [0/sf_d floor(size(score(:,1),1)/sf_d)] ;
cbh.TickLabels = [0/sf_d floor(size(score(:,1),1)/7.73)-1] ;
set(gca,'fontsize',16)


clear X Y spikes spikes_d_s score mouse i FRp FR coeff cells_d T


%V1
mouse='92229';
load([rec_path,strcat('recording_dates_',mouse,'.mat')]);
day=7;
s=1;
munit=dates.ses{day}(s);
file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name,'-mat');
spikes=full(spikes_d_s);
[N,T]=size(spikes);
sf=7.73;
sf_d=7.73/4;

spikes_full=spikes;
clear spikes
spikes=spikes_full(:,1:T_PaS);



% FRp = spikes_downsample(spikes,N,4);
FRp = spikes_downsample(spikes,N,1);
%LEM
for i=1:N
%     FR(i,:)=(full(fire_rate(FRp(i,:),29,'g')));  
    FR(i,:)=(full(fire_rate(FRp(i,:),4*floor(66/1),'g')));

end
clear score
[coeff,score]=pca(FR');

addpath("C:\Users\xscogno\drtoolbox.tar");
[Y,X] = laplacian_eigen([score(:,1:5)],2,15,2);

figure
plot3((1:size(Y,1)),-Y(1:end,1),-Y(1:end,2),'Linewidth',3)
colormap viridis
cla
patch([1:size(Y,1) nan],[-Y(:,1)' nan],[-Y(:,2)' nan],[1:size(Y,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('Time (s)')
ylabel('LEM - Dimension 1');
zlabel('LEM - Dimension 2');
yticks([]);
zticks([]);
set(gca,'fontsize',16)
xticks([0 1546 3092])
xticklabels([0 1546/sf_d 3092/sf_d])
axis([0 3500 -inf inf -inf inf])
caxis([0 3500]);

%LEM projection V1

figure
plot(-Y(1:end,1),-Y(1:end,2),'Linewidth',3)
colormap viridis
cla
patch([-Y(:,1)' nan],[-Y(:,2)' nan],[1:size(Y,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('LEM - Dimension 1');
ylabel('LEM - Dimension 2');
yticks([]);
xticks([]);
set(gca,'fontsize',16)
box off

rmpath("C:\Users\xscogno\drtoolbox.tar");

figure
plot3((1:size(score,1)),score(1:end,1),score(1:end,2),'Linewidth',3)
axis([0 8000 -inf inf -inf inf])
% axis([0 2000 -inf inf -inf inf])
colormap viridis
cla
patch([1:size(score,1) nan],[score(:,1)' nan],[score(:,2)' nan],[1:size(score,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('Time (s)')
ylabel('PC1');
zlabel('PC2');
yticks([]);
zticks([]);
set(gca,'fontsize',16)
xticks([0 1546 3092])
xticklabels([0 1546/sf_d 3092/sf_d])
axis([0 3500 -inf inf -inf inf])
caxis([0 3500]);
cbh = colorbar ; %Create Colorbar
cbh.Ticks = [0 3092] ; %Create 8 ticks from zero to 1
cbh.TickLabels = [0/sf_d 3092/sf_d] ;
set(gca,'fontsize',16)

%PCA projection V1
figure
plot(score(1:end,1),score(1:end,2),'Linewidth',3)
    alfa=max([abs(score(:,2));abs(score(:,1))]);
    axis([-alfa-0.05 alfa+0.05 -alfa-0.05 alfa+0.05]);
% axis([-40 120 -40 120])
% axis([0 2000 -inf inf -inf inf])
colormap viridis
cla
patch([score(:,1)' nan],[score(:,2)' nan],[1:size(score,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('PC1');
ylabel('PC2');
xticks([]);
yticks([]);
set(gca,'fontsize',16)
box off
% c=colorbar( 'XTick', [400 1400]);
cbh = colorbar ; %Create Colorbar
cbh.Ticks = [1 size(score(:,1),1)] ; %Create 8 ticks from zero to 1
% cbh.TickLabels = [0/sf_d floor(size(score(:,1),1)/sf_d)] ;
cbh.TickLabels = [0/sf_d floor(size(score(:,1),1)/7.73)-1] ;
set(gca,'fontsize',16)
%% 2D Projection L8M2

figure
plot(atan(score(:,1)./score(:,2)));

win=[1151:1194];
x_mean=sum(score(win,1))./length(score(win,1));
y_mean=sum(score(win,2))./length(score(win,2));
figure
plot(score(win,1)-x_mean,score(win,2)-y_mean,'k','linewidth',3);
box off
axis square
axis([-50 50 -50 50])
ylabel('PC2');
xlabel('PC1');
set(gca,'fontsize',16)

win=[1202:1245];
figure
plot(score(win,1),score(win,2),'k','linewidth',3);
box off
axis square
axis([-40 40 -40 40])

win=[1245:1295];
figure
plot(score(win,1),score(win,2),'k','linewidth',3);
box off
axis square
axis([-40 40 -40 40])

%% 2D Projection L9M4

clear all
% close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Behavior L8\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';


% 9 seconds for L8M2 Day 19 ; 20 seconds for L9M4 Day 17 
mouse='L9M4';
day=17;
munit=0;

load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
file_name=[dpath_spikes ['spikes_120ms_Do_THR1_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name,'-mat');
spikes=full(spikes_d_s);
N=size(spikes,1);
sf=7.73;
sf_d=7.73/4;

[N,T]=size(spikes);

ensembles=N;
% FRp = spikes_downsample(spikes,N,8);
FRp = spikes_downsample(spikes,N,4);
for i=1:ensembles
    FR(i,:)=(full(fire_rate(FRp(i,:),4*29,'g')));  % 9 seconds for L8M2 Day 19 ; 15 seconds for L9M4 Day 17 
end
[coeff,score]=pca(FR');

%PCA projection
figure
plot(score(1:end,1),score(1:end,2),'Linewidth',3)
axis([-100 100 -80 80])
% axis([0 2000 -inf inf -inf inf])
colormap viridis
cla
patch([score(:,1)' nan],[score(:,2)' nan],[1:size(score,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('PC1');
ylabel('PC2');
xticks([]);
yticks([]);
set(gca,'fontsize',16)
box off
% c=colorbar( 'XTick', [400 1400]);
cbh = colorbar ; %Create Colorbar
cbh.Ticks = [1 size(score(:,1),1)] ; %Create 8 ticks from zero to 1
cbh.TickLabels = [0/sf_d floor(size(score(:,1),1)/sf_d)] ;
set(gca,'fontsize',16)

%LEM
addpath("C:\Users\xscogno\MATLAB\drtoolbox.tar");
[Y,X] = laplacian_eigen([score(:,1:5)],2,15,2);
figure
plot(-Y(1:end,1),-Y(1:end,2),'Linewidth',3)
colormap viridis
cla
patch([-Y(:,1)' nan],[-Y(:,2)' nan],[1:size(Y,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('LEM - Dimension 1');
ylabel('LEM - Dimension 2');
yticks([]);
xticks([]);
set(gca,'fontsize',16)
box off
% xticks([0 3865 7730])
% xticklabels([0 3865/sf_d 7730/sf_d])
% axis([0 8000 -inf inf -inf inf])

rmpath("C:\Users\xscogno\drtoolbox.tar");

%Only one cycle

% figure
% plot(atan(score(:,1)./score(:,2)));

win=[1:6958];
x_mean=sum(score(win,1))./length(score(win,1));
y_mean=sum(score(win,2))./length(score(win,2));
figure
plot(score(win,1)-x_mean,score(win,2)-y_mean,'k','linewidth',3);
box off
axis square
axis([-100 100 -100 100])
ylabel('PC2');
xlabel('PC1');
set(gca,'fontsize',18)
xticks([])
yticks([])
grid on
hold on
plot([-100:10:100],ones(1,21)*0,'--','color',[0.45, 0.45, 0.45],'linewidth',2);
plot(ones(1,21)*0,[-100:10:100],'--','color',[0.45, 0.45, 0.45],'linewidth',2);

figure
% win=[1:6900];
sf=7.73;
bins_size=(1/sf)*4;
win=[1:6958];
duration=length(win)*bins_size;
section(:,1)=score(win,1);
section(:,2)=score(win,2);
section_i(:,1) = interp1(1:length(win),section(:,1),1:0.1:length(win)) ;
section_i(:,2) = interp1(1:length(win),section(:,2),1:0.1:length(win)) ;
cc=viridis(size(section_i,1));
x_mean=sum(score(win,1))./length(score(win,1));
y_mean=sum(score(win,2))./length(score(win,2));
figure
hold on
for i=1:size(section_i,1)
scatter(section_i(i,1)-x_mean,section_i(i,2) -y_mean,10,cc(i,:),'filled');
end
box off
axis square
axis([-100 100 -100 100])
ylabel('PC2');
xlabel('PC1');
set(gca,'fontsize',18)
xticks([])
yticks([])
grid on
hold on
plot([-100:10:100],ones(1,21)*0,'--','color',[0.45, 0.45, 0.45],'linewidth',2);
plot(ones(1,21)*0,[-100:10:100],'--','color',[0.45, 0.45, 0.45],'linewidth',2);
colormap viridis
c=colorbar( 'XTick', [0 1],'TickLabels',[0 duration]);

sf=7.73;
bins_size=(1/sf)*4;
win=[1:6958];
duration=length(win)*bins_size;
section(:,1)=score(win,1);
section(:,2)=score(win,2);
cc=viridis(size(section,1));
x_mean=sum(score(win,1))./length(score(win,1));
y_mean=sum(score(win,2))./length(score(win,2));
figure
hold on
for i=1:size(section,1)
scatter(section(i,1)-x_mean,section(i,2) -y_mean,10,cc(i,:),'filled');
end
box off
axis square
axis([-150 150 -150 150])
ylabel('PC2');
xlabel('PC1');
set(gca,'fontsize',18)
xticks([])
yticks([])
grid on
hold on
plot([-150:10:150],ones(1,31)*0,'--','color',[0.45, 0.45, 0.45],'linewidth',2);
plot(ones(1,31)*0,[-150:10:150],'--','color',[0.45, 0.45, 0.45],'linewidth',2);
colormap viridis
c=colorbar( 'XTick', [0 1],'TickLabels',[0 floor(duration)]);
axis square

%% 2D Projection L9M4 NEW

clear all
% close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Behavior L8\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';


% 9 seconds for L8M2 Day 19 ; 20 seconds for L9M4 Day 17 
mouse='L9M4';
day=17;
munit=0;

load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
file_name=[dpath_spikes ['spikes_120ms_Do_THR1_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name,'-mat');
spikes=full(spikes_d_s);
N=size(spikes,1);
sf=7.73;
sf_d=7.73/4;

[N,T]=size(spikes);

ensembles=N;
% FRp = spikes_downsample(spikes,N,8);
FRp = spikes_downsample(spikes,N,4);
for i=1:ensembles
    FR(i,:)=(full(fire_rate(spikes(i,:),1*117,'g')));  % 9 seconds for L8M2 Day 19 ; 15 seconds for L9M4 Day 17 
end
[coeff,score]=pca(FR');

%PCA projection
figure
plot(score(1:end,1),score(1:end,2),'Linewidth',3)
axis([-100 100 -80 80])
% axis([0 2000 -inf inf -inf inf])
colormap viridis
cla
patch([score(:,1)' nan],[score(:,2)' nan],[1:size(score,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('PC1');
ylabel('PC2');
xticks([]);
yticks([]);
set(gca,'fontsize',16)
box off
% c=colorbar( 'XTick', [400 1400]);
cbh = colorbar ; %Create Colorbar
cbh.Ticks = [1 size(score(:,1),1)] ; %Create 8 ticks from zero to 1
cbh.TickLabels = [0/sf_d floor(size(score(:,1),1)/sf_d)] ;
set(gca,'fontsize',16)

%LEM
addpath("C:\Users\xscogno\MATLAB\drtoolbox.tar");
[Y,X] = laplacian_eigen([score(:,1:5)],2,15,2);
figure
plot(-Y(1:end,1),-Y(1:end,2),'Linewidth',3)
colormap viridis
cla
patch([-Y(:,1)' nan],[-Y(:,2)' nan],[1:size(Y,1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
axis on
xlabel('LEM - Dimension 1');
ylabel('LEM - Dimension 2');
yticks([]);
xticks([]);
set(gca,'fontsize',16)
box off
% xticks([0 3865 7730])
% xticklabels([0 3865/sf_d 7730/sf_d])
% axis([0 8000 -inf inf -inf inf])

rmpath("C:\Users\xscogno\drtoolbox.tar");

%Only one cycle

% figure
% plot(atan(score(:,1)./score(:,2)));

win=[1:6958];
x_mean=sum(score(win,1))./length(score(win,1));
y_mean=sum(score(win,2))./length(score(win,2));
figure
plot(score(win,1)-x_mean,score(win,2)-y_mean,'k','linewidth',3);
box off
axis square
axis([-100 100 -100 100])
ylabel('PC2');
xlabel('PC1');
set(gca,'fontsize',18)
xticks([])
yticks([])
grid on
hold on
plot([-100:10:100],ones(1,21)*0,'--','color',[0.45, 0.45, 0.45],'linewidth',2);
plot(ones(1,21)*0,[-100:10:100],'--','color',[0.45, 0.45, 0.45],'linewidth',2);

figure
% win=[1:6900];
sf=7.73;
bins_size=(1/sf)*4;
win=[1:6958];
duration=length(win)*bins_size;
section(:,1)=score(win,1);
section(:,2)=score(win,2);
section_i(:,1) = interp1(1:length(win),section(:,1),1:0.1:length(win)) ;
section_i(:,2) = interp1(1:length(win),section(:,2),1:0.1:length(win)) ;
cc=viridis(size(section_i,1));
x_mean=sum(score(win,1))./length(score(win,1));
y_mean=sum(score(win,2))./length(score(win,2));
figure
hold on
for i=1:size(section_i,1)
scatter(section_i(i,1)-x_mean,section_i(i,2) -y_mean,10,cc(i,:),'filled');
end
box off
axis square
axis([-100 100 -100 100])
ylabel('PC2');
xlabel('PC1');
set(gca,'fontsize',18)
xticks([])
yticks([])
grid on
hold on
plot([-100:10:100],ones(1,21)*0,'--','color',[0.45, 0.45, 0.45],'linewidth',2);
plot(ones(1,21)*0,[-100:10:100],'--','color',[0.45, 0.45, 0.45],'linewidth',2);
colormap viridis
c=colorbar( 'XTick', [0 1],'TickLabels',[0 duration]);

sf=7.73;
bins_size=(1/sf)*4;
win=[1:6958];
duration=length(win)*bins_size;
section(:,1)=score(win,1);
section(:,2)=score(win,2);
cc=viridis(size(section,1));
x_mean=sum(score(win,1))./length(score(win,1));
y_mean=sum(score(win,2))./length(score(win,2));
figure
hold on
for i=1:size(section,1)
scatter(section(i,1)-x_mean,section(i,2) -y_mean,10,cc(i,:),'filled');
end
box off
axis square
axis([-150 150 -150 150])
ylabel('PC2');
xlabel('PC1');
set(gca,'fontsize',18)
xticks([])
yticks([])
grid on
hold on
plot([-150:10:150],ones(1,31)*0,'--','color',[0.45, 0.45, 0.45],'linewidth',2);
plot(ones(1,31)*0,[-150:10:150],'--','color',[0.45, 0.45, 0.45],'linewidth',2);
colormap viridis
c=colorbar( 'XTick', [0 1],'TickLabels',[0 floor(duration)]);
axis square

%% Ring for all wave sessions dutin uninterrupted oscillations


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

%Wave sessions


figpath='C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Raster Plots all MEC sessions\';
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

% Phase and ring
figure
countf=0;
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

    load(file_name_spk,'-mat'); %Spike times
    spikes_d=full(spikes_d_s);
    [N,T]=size(spikes_d);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Condition on having waves
    sf=7.73;
    sf_d=7.73/4;
    num_clus_discr=10;
    make_fig=0;
    dt=floor(big_table(row_w,8));
    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes_d);

    for i=2:size(table_u,1)
        ICI_ses(i-1)=(table_u(i,1)-table_u(i-1,2));
    end

    epoch=[];
    count=0;
    for i=1:size(ICI_ses,2)
        disp(i)
        if ICI_ses(i)==1
            count=count+1;
            epoch(count,1)=i;
            for j=i+1:size(ICI_ses,2)
                if ICI_ses(j)==1
                else
                    epoch(count,2)=j+1;
                    break
                    %                     i=j+1;
                end
            end

        end
    end

    aux=find(epoch(:,2)==0);
    epoch(aux,2)=size(table_u,1);
    [a,r]=max(abs(epoch(:,2)-epoch(:,1)))

    FRp = spikes_downsample(spikes_d,N,4);

    for i=1:N
        FR(i,:)=(full(fire_rate(FRp(i,:),2*floor(dt/4),'g')));
    end
    clear score
    [coeff,score]=pca(FR');
    window=(floor(table_u(epoch(r,1),1)/4):floor(table_u(epoch(r,2),2)/4));
%     subplot(4,4,w);
%     plot(atan2(score(window,2),score(window,1)))
    window_1cycle=(floor(table_u(epoch(r,1)+1,1)/4):floor(table_u(epoch(r,1)+1,2)/4));
    countf=countf+1;
    %PCA projection PaS
    subplot(5,6,countf);
    plot(smooth(score(window_1cycle,1),4*floor(dt/4)),smooth(score(window_1cycle,2),4*floor(dt/4)),'Linewidth',3)
    alfa=max([score(window_1cycle,2);score(window_1cycle,1)]);
    axis([-alfa-5 alfa+5 -alfa-5 alfa+5]);
    % axis([0 2000 -inf inf -inf inf])
    colormap viridis
    cla
    patch([score(window_1cycle,1)' nan],[score(window_1cycle,2)' nan],[1:size(score(window_1cycle,:),1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
    axis on
    xlabel('PC1');
    ylabel('PC2');
    xticks([]);
    yticks([]);
    set(gca,'fontsize',16)
    box off
    % c=colorbar( 'XTick', [400 1400]);
    cbh = colorbar ; %Create Colorbar
    cbh.Ticks = [1 size(score(window,1),1)] ; %Create 8 ticks from zero to 1
    cbh.TickLabels = [0/sf_d floor(size(score(window,1),1)/sf_d)] ;
    set(gca,'fontsize',16)
    title(num2str(epoch(r,2)-epoch(r,1)));
    countf=countf+1;
    subplot(5,6,countf);
    plot(smooth(score(window,1),4*floor(dt/4)),smooth(score(window,2),4*floor(dt/4)),'Linewidth',3)
    alfa=max([score(window,2);score(window,1)]);
    axis([-alfa-5 alfa+5 -alfa-5 alfa+5]);
    % axis([0 2000 -inf inf -inf inf])
    colormap viridis
    cla
    patch([score(window,1)' nan],[score(window,2)' nan],[1:size(score(window,:),1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
    axis on
    xlabel('PC1');
    ylabel('PC2');
    xticks([]);
    yticks([]);
    set(gca,'fontsize',16)
    box off
    % c=colorbar( 'XTick', [400 1400]);
    cbh = colorbar ; %Create Colorbar
    cbh.Ticks = [1 size(score(window,1),1)] ; %Create 8 ticks from zero to 1
    cbh.TickLabels = [0/sf_d floor(size(score(window,1),1)/sf_d)] ;
    set(gca,'fontsize',16)
    title(num2str(epoch(r,2)-epoch(r,1)));

    clear X Y spikes spikes_d_s score mouse i FRp FR coeff cells_d T window score FR FRp epoch ICI_ses table_u spikes spikes_d ICI_ses_t

end

%% %% Ring for all wave sessions for all oscillations


clear all
% close all

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

%Wave sessions

figpath='C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Raster Plots all MEC sessions\';
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

% Phase and ring
figure
countf=0;
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

    load(file_name_spk,'-mat'); %Spike times
    spikes_d=full(spikes_d_s);
%     [~,~]=size(spikes_d);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Condition on having waves
    sf=7.73;
%     sf_d=7.73/4;
    num_clus_discr=10;
    make_fig=0;
    dt=floor(big_table(row_w,8));
    [table_u,N,~]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes_d);


    spikes_r=[];
    for wa=1:size(table_u,1)
        spikes_r=[spikes_r,spikes_d(:,table_u(wa,1):table_u(wa,2))];
        %phase_r=[phase_r;phase_f(table_u(wa,1):table_u(wa,2))];
    end

    ds=1;
    sf_d=7.73/ds;
    time_points=[];
    for i=1:size(table_u,1)
        if (floor(table_u(i,1)/ds)==0)
            time_points=[time_points,1:floor(table_u(i,2)/ds)];

        else
            time_points=[time_points,floor(table_u(i,1)/ds):floor(table_u(i,2)/ds)];
        end
    end

    FRp = spikes_downsample(spikes_d,N,ds);

    if w==7
        for i=1:N
            FR(i,:)=(full(fire_rate(FRp(i,:),2*floor(dt/1),'g')));
        end
    else
        for i=1:N
            FR(i,:)=(full(fire_rate(FRp(i,:),4*floor(dt/1),'g')));
        end
    end

    clear score
    [~,score]=pca(FR');

%     score(:,1)=smooth(score(:,1),4*floor(dt/1));
%     score(:,2)=smooth(score(:,2),4*floor(dt/1));

    countf=countf+1;


%     seq2show=size(table_u,1)-3;

    %PCA projection PaS

    subplot(4,4,countf)
%     figure
    plot(score(time_points,1),score(time_points,2),'Linewidth',3);
    alfa=max([abs(score(:,2));abs(score(:,1))]);
    axis([-alfa-0.05 alfa+0.05 -alfa-0.05 alfa+0.05]);
    colormap viridis
    cla
    patch([score(time_points,1)' nan],[score(time_points,2)' nan],[1:size(score(time_points,:),1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
    axis on
    xlabel('PC1');
    ylabel('PC2');
    xticks([]);
    yticks([]);
    set(gca,'fontsize',16,'YColor','k','XColor','k')
    box off
    cbh = colorbar ; %Create Colorbar
    cbh.Ticks = [1 length(score(time_points,1))] ; %Create 8 ticks from zero to 1
    cbh.TickLabels = [0/sf_d floor(length(score(time_points,1))/sf_d/60)] ;
%     title(num2str(epoch(r,2)-epoch(r,1)));
    axis square
%     saveas(gcf,['C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures March\Rings\Ring_',num2str(mouse),'day',num2str(day),'.fig']);
%     saveas(gcf,['C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures March\Rings\Ring_',num2str(mouse),'day',num2str(day),'.svg']);

%     close all


    clear X Y spikes spikes_d_s score mouse i FRp FR coeff cells_d T window score FR FRp epoch ICI_ses table_u spikes spikes_d ICI_ses_t
 clear time_points  spikes_r
end


%% PC1 PC2 ring shaped manifold for example session


clear all
% close all

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

%Wave sessions

figpath='C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Raster Plots all MEC sessions\';
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

% Phase and ring
figure
countf=0;
for w=8%1:length(waves)
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

    load(file_name_spk,'-mat'); %Spike times
    spikes_d=full(spikes_d_s);
%     [~,~]=size(spikes_d);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Condition on having waves
    sf=7.73;
%     sf_d=7.73/4;
    num_clus_discr=10;
    make_fig=0;
    dt=floor(big_table(row_w,8));
    [table_u,N,~]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes_d);


    spikes_r=[];
    for wa=1:size(table_u,1)
        spikes_r=[spikes_r,spikes_d(:,table_u(wa,1):table_u(wa,2))];
        %phase_r=[phase_r;phase_f(table_u(wa,1):table_u(wa,2))];
    end

    ds=1;
    sf_d=7.73/ds;
    time_points=[];
    for i=1:size(table_u,1)
        if (floor(table_u(i,1)/ds)==0)
            time_points=[time_points,1:floor(table_u(i,2)/ds)];

        else
            time_points=[time_points,floor(table_u(i,1)/ds):floor(table_u(i,2)/ds)];
        end
    end

    FRp = spikes_downsample(spikes_d,N,ds);

    if w==7
        for i=1:N
            FR(i,:)=(full(fire_rate(FRp(i,:),2*floor(dt/1),'g')));
        end
    else
        for i=1:N
            FR(i,:)=(full(fire_rate(FRp(i,:),4*floor(dt/1),'g')));
        end
    end

    clear score
    [~,score]=pca(FR');

%     score(:,1)=smooth(score(:,1),4*floor(dt/1));
%     score(:,2)=smooth(score(:,2),4*floor(dt/1));

    countf=countf+1;


%     seq2show=size(table_u,1)-3;

    %PCA projection PaS

    figure
    plot(score(:,1),score(:,2),'Linewidth',3);
    alfa=max([abs(score(:,2));abs(score(:,1))]);
    axis([-alfa-0.05 alfa+0.05 -alfa-0.05 alfa+0.05]);
    colormap viridis
    cla
    patch([score(:,1)' nan],[score(:,2)' nan],[1:size(score(:,:),1) nan],'Linewidth',5,'EdgeColor','interp','FaceColor','none')
    axis on
    xlabel('PC1');
    ylabel('PC2');
    xticks([]);
    yticks([]);
    set(gca,'fontsize',16,'YColor','k','XColor','k')
    box off
    cbh = colorbar ; %Create Colorbar
    cbh.Ticks = [1 length(score(:,1))] ; %Create 8 ticks from zero to 1
    cbh.TickLabels = [0/sf_d floor(length(score(:,1))/sf_d)] ;
%     title(num2str(epoch(r,2)-epoch(r,1)));
    axis square
%     saveas(gcf,['C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures March\Rings\Ring_',num2str(mouse),'day',num2str(day),'.fig']);
%     saveas(gcf,['C:\Users\xscogno\Dropbox\Waves\Manuscript iteration\Figures March\Rings\Ring_',num2str(mouse),'day',num2str(day),'.svg']);

%     close all


    clear X Y spikes spikes_d_s score mouse i FRp FR coeff cells_d T window score FR FRp epoch ICI_ses table_u spikes spikes_d ICI_ses_t
 clear time_points  spikes_r
end
