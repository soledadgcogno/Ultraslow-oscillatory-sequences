%Fig 3a,c,e
%Extended data Fig 7g
%Extended data Fig 10i
%Extended dara Fig 10 j,k
%Extended data Fig 7d
%Fig 3d
%Extended data Fig10c,d

%% Loads data
clear all
close all

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

%Load files
load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
munit=dates.ses{day}(s);
file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']]; %SNR
file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']]; %Spikes
file_name_dff=[dpath ['DFF_120ms_Do_SNRH','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']]; %DFF
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']]; %Anat
file_dt=[dpath ['WS_Osc_14_',mouse_name,'.mat']]; %dt


SNR=load(file_name_snr,'-mat'); %SNR
snr=SNR.SNR_dff; %SNR
load(file_name_anat,'-mat'); %Anatomical information
load(file_name_dff,'-mat'); %DFF
load(file_name_spk,'-mat'); %Spike times
load(file_dt,'-mat'); %dt
spikes_d=full(spikes_d_s);
[N,T]=size(spikes_d);
dff=signal_dff;


%% Calculate and condition spikes and phase of wave on wave epochs

dt=WS_stat.dt(day,s);
if isinteger(dt)
else
    dt=floor(dt);
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikes=spikes_d;
phase=phase_r;
spikes_d=[];
spikes_d=spikes_r;

% figure
% subplot(2,1,1)
% plot(phase_f)
% hold on
% for i=1:size(table_u)
%     plot(table_u(i,1):table_u(i,2),phase_f(table_u(i,1):table_u(i,2)),'r','linewidth',1.5);
% end
% subplot(2,1,2)
% plot(phase_f)
% hold on
% for i=1:size(table_u)
%     plot(table_u_corrected(i,1):table_u_corrected(i,2),phase_f(table_u_corrected(i,1):table_u_corrected(i,2)),'r','linewidth',1.5);
% end

%% Compute tuning curves and locking using deconvolved spikes

N_sh=10; %in the paper we used 1000               
for i=1:N
    p=phase(find(spikes_d(i,:)));
    var_p(i)=circ_var(p);
    mean_p(i)=circ_mean(p);
    std_p(i)=circ_std(p);
    MVL(i) = circ_r(p);
    prob_phase_firing(i,:)=histcounts(p,-pi:2*pi/40:pi,'Normalization','Probability');
    clear p H
end
i=1;
p=phase(find(spikes_d(i,:)));
H=histogram(p,-pi:2*pi/40:pi,'Normalization','Probability');
phase_TC=H.BinEdges;
close all

for i=1:N
    disp(i)
    for sh=1:N_sh
        p_sh=phase(randperm(length(phase_r),length(find(spikes_d(i,:)))));
        var_p_sh(i,sh)=circ_var(p_sh);
        std_p_sh(i,sh)=circ_std(p_sh);
        mean_p_sh(i,sh)=circ_mean(p_sh);
        MVL_sh(i,sh) = circ_r(p_sh);
        if sh==1
            prob_phase_firing_sh(i,:)=histcounts(p_sh,-pi:2*pi/40:pi,'Normalization','Probability');
            clear H
        end
%         H=histogram(p_sh,-pi:2*pi/40:pi,'Normalization','Probability');
%         prob_phase_firing_sh_tot(i,:,sh)=H.Values;
%         clear H
                clear p_sh
    end
end
                
locking=MVL;
locking_sh=MVL_sh;
locking_sh_mean=mean(MVL_sh,2);
locking_sh_99=prctile(MVL_sh,99,2);
locking_sh_1=prctile(MVL_sh,1,2);

locked=find(MVL>locking_sh_99');
fraction_locked=length(locked)/N;
not_locked=find(MVL<=locking_sh_99');
fraction_notlocked=length(not_locked)./N;

%Keep only locked cells
N_locked=length(locked);
N_not_locked=length(not_locked);

mean_p_locked=mean_p(locked);
std_p_locked=std_p(locked);

mean_p_sh_locked=mean_p_sh(locked,:);
std_p_sh_locked=std_p_sh(locked,:);

prob_phase_firing_locked=prob_phase_firing(locked,:);
prob_phase_firing_sh_locked=prob_phase_firing_sh(locked,:);

tuning_curve=prob_phase_firing_locked;
tuning_curve_sh=prob_phase_firing_sh_locked;

%% Figures: Locking and tuning curves (Figure 3a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pie chart of locked and not locked

figure
labels = {'Locked','Not locked'};
pie([fraction_locked,fraction_notlocked],labels)
colormap jet(4)
set(gca,'fontsize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MVL

[a,b]=sort(MVL,'ascend');
figure
hold on
plot(locking_sh_99(b),'.','MarkerSize',15,'Color','r');
plot(MVL(b),'.','MarkerSize',15,'Color','k');
alpha 0.6
ylabel('Locking to phase');
xlabel('Neuron #')
set(gca,'fontsize',18)
box off
legend('Shuffle - 99th percentile', 'Data')
legend boxoff
axis([0 530 0 1])
xticks([100 400])
yticks([0 0.5 1])
axis([-inf inf 0 1])
axis square


for k=1:length(not_locked)
    index_notlocked(k)=find(not_locked(k)==b);
end

[a,b]=sort(MVL,'ascend');
figure
hold on
plot(locking_sh_99(b),'.','MarkerSize',15,'Color',[160 160 160]/255);
plot(MVL(b),'.','MarkerSize',15,'Color','k');
plot(index_notlocked,MVL(not_locked),'.','MarkerSize',15,'Color','r');
alpha 0.6
ylabel('Locking to phase');
xlabel('Neuron #')
set(gca,'fontsize',18)
box off
legend('Shuffle - 99th percentile', 'Data')
legend boxoff
axis([0 530 0 1])
xticks([100 400])
yticks([0 0.5 1])
axis([-inf inf 0 1])
axis square
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preferred phase (mean phase of locked cells)

[~,sorting_mean_angle]=sort(mean_p_locked,'ascend');
N_n=size(locked,2);


figure
x=1:size(mean_p_locked,2);
x2=[x,fliplr(x)];
inBetween = [mean_p_locked(sorting_mean_angle)+std_p_locked(sorting_mean_angle), fliplr(mean_p_locked(sorting_mean_angle)-std_p_locked(sorting_mean_angle))];
fill(x2, inBetween, [252 215 215]./255);
hold on
plot(1:N_n,mean_p_locked(sorting_mean_angle)+std_p_locked(sorting_mean_angle),'w','linewidth',1.5);
hold on
plot(mean_p_locked(sorting_mean_angle)-std_p_locked(sorting_mean_angle),'w','linewidth',1.5);
plot(mean_p_locked(sorting_mean_angle),'k','linewidth',3);
axis([2 N_n-1 -5 5]);
yticks([-3.14,3.14])
yticklabels({'-\pi','\pi'});
xticks([100 400])
ylabel({'Preferred phase (rad)'})
xlabel('(Locked) Neuron #')
set(gca,'fontsize',18)
box off



% figure
% plot(mean_p_locked, std_p_locked,'k.','Markersize',15);
% ylabel({'Variance of phase';'(Width of tuning curve)'});
% xlabel({'Mean phase';'(Peak of tuning curve)'});
% set(gca,'fontsize',18)
% box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tuning curves of single cell to the phase using as sorting criterion the mean phase at which cells fired


[~,sorting_mean_angle]=sort(mean_p_locked,'ascend');
figure
imagesc(prob_phase_firing_locked(sorting_mean_angle,:));
xticks([1,20,40])
xticklabels({'-3.14','0','3.14'})
xlabel('Phase (rad)');
ylabel('Cell');
c=colorbar;
c.Label.String='Probability';
% caxis([0 0.2]);
% colorbar('Ticks',[0,0.1,0.2]);
set(gca,'fontsize',18)
colormap inferno
% co=colorbar;
caxis([0 0.15]);
co=colorbar('Ticks',[0,0.15,0.2]);
set(gca,'fontsize',18)
co.Label.String = 'Fraction of event counts';


[~,sorting_mean_angle_sh]=sort(mean_p_sh_locked(:,1),'ascend');
figure
imagesc(prob_phase_firing_sh(sorting_mean_angle_sh,:));
xticks([1,20,40])
xticklabels({'-3.14','0','3.14'})
xlabel('Phase (rad)');
ylabel('Cell');
colormap inferno
caxis([0 0.15]);
co=colorbar('Ticks',[0,0.15]);
%                 c.Label.String='Probability';
set(gca,'fontsize',18)
co.Label.String = 'Fraction of event counts';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mean phase and maximum of tuning curve
figure
plot(b(sorting_max_tc),mean_p_locked(sorting_max_tc),'*');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distribution of phase that maximizes the tuning curve in anatomical cortex

% tuning_curve=prob_phase_firing;
% not_locked=find(locking<=locking_sh_99');
[~,b]=max(prob_phase_firing');
% [~,b2]=sort(b,'ascend');
cc=inferno(40);
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_anat,'-mat');
figure
hold on
pixel_size_new=1.18185;
for i=1:N
    x=Anat.pos{1,i}(:,1);
    y=Anat.pos{1,i}(:,2);
    ov=Anat.overlap{1,i};
    
    if ismember(i,not_locked)==1
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
    else
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(b(i),:),'filled');
    end
    clear x y ov
end
box on
xlabel('X [\mum]');
ylabel('Y [\mum]');
axis([0 480 20 500])
yticks([120 420]);
yticklabels([100 400]);
xticks([100 400]);
set(gca,'fontsize',16)
axis square
colormap inferno(40)
co=colorbar('XTick',0:1,'XTickLabel',{'-\pi','\pi'});
co.Label.String = 'Preferred phase';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distribution of mean phase in anatomical cortex

% tuning_curve=prob_phase_firing;
% not_locked=find(locking<=locking_sh_99');
% [~,b2]=sort(mean_p,'ascend');
[b,E]=discretize(mean_p,40);
cc= cividis(40);
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_anat,'-mat');
pixel_size_new=1.18185;

figure
hold on
for i=1:N
    x=Anat.pos{1,i}(:,1);
    y=Anat.pos{1,i}(:,2);
    ov=Anat.overlap{1,i};
    
    if ismember(i,not_locked)==1
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
    else
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(b(i),:),'filled');
    end
        clear x y ov
        
end
box on
xlabel('X (\mum)');
ylabel('Y (\mum)');
axis([0 580 20 600])
yticks([120 420]);
yticklabels([100 400]);
xticks([100 400]);
set(gca,'fontsize',16)
axis square
colormap  cividis(40)%gmt_nighttime(40)
co=colorbar('XTick',0:1,'XTickLabel',{'-\pi','\pi'});
co.Label.String = 'Preferred phase';

clear b E

% colormap puor
% colormap gmt_ocean;
% colormap cividis
% colormap gmt_nighttime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% firing rate

% tuning_curve=prob_phase_firing;
% not_locked=find(locking<=locking_sh_99');
% [~,b2]=sort(mean_p,'ascend');
mean_spike_number=sum(spikes_d,2)/(T/8);
[b,E]=discretize(mean_spike_number,10);
cc= gmt_nighttime(10);
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_anat,'-mat');
pixel_size_new=1.18185;

figure
hold on
for i=1:N
    x=Anat.pos{1,i}(:,1);
    y=Anat.pos{1,i}(:,2);
    ov=Anat.overlap{1,i};
    
    if ismember(i,not_locked)==1
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
    else
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(b(i),:),'filled');
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
colormap  gmt_nighttime(10)
co=colorbar('XTick',0:1,'XTickLabel',{'0.06','0.4'});
co.Label.String = 'Mean Deconvolved spikes (Hz)';

clear b E


%% Examples of individual cells (Extended data Fig 7g)

wave_epochs=-3*zeros(1,T);
for w=1:size(table_u,1)
    wave_epochs(table_u(w,1):table_u(w,2))=1;
end

col_dots=[0.9283    0.4730    0.3261];%[0.6107    0.0902    0.6200];
% col_area=[214 214 245]./255;%[0.0504,0.0298,0.5280]*0.5;
col_area=[255 190 51]./255;%[0.0504,0.0298,0.5280]*0.5;
% col_area=[255 92 51]./255;%[0.0504,0.0298,0.5280]*0.5;

cell_high_locking=18; %MVL=0.94
cell_high_locking_2=30; %MVL=0.91
cell_high_locking_3=474; %MVL=0.90
cell_med_locking=97; %MVL=0.62
cell_low_locking_not_locked=43; %MVL=0.1
cell_low_locking=217; %MVL=0.4
cell_low_locking_not_locked2=449; %MVL=0.15
cell_med_locking_2=86; %MVL=0.59



%%%%%%%%%%%%%%%%%Figures with more than one example - 1
[sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
ind_cell3=find(sorting_ascend==cell_high_locking_3);
ind_cell2=find(sorting_ascend==cell_high_locking_2);
ind_cell1=find(sorting_ascend==cell_high_locking);

figure
subplot(4,1,1)
hold on
for i=1:size(spikes,1)
    %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_descend(i),:),5,'k','filled')
    scatter((1:size(spikes,2))./8,i*spikes(sorting_ascend(i),:),5,'k','filled')
    alpha 0.1
end
scatter((1:size(spikes,2))./8,ind_cell3.*spikes(cell_high_locking_3,:),40,[183,212,219]/255,'filled')
scatter((1:size(spikes,2))./8,ind_cell2.*spikes(cell_high_locking_2,:),40,[149,125,173]/255,'filled')
scatter((1:size(spikes,2))./8,ind_cell1.*spikes(cell_high_locking,:),40,[244,209,181]/255,'filled')
axis([(T/2)/8 inf 2 inf])
ylabel('Neurons #');
yticks([100 400])
set(gca,'XColor', 'none')
set(gca,'fontsize',18)
subplot(4,1,2)
axis([(T/2)/8  inf -2 14])
aux=(abs(phase_f-mean_p(cell_high_locking)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[244,209,181]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking),:)),'Color','k','linewidth',1.5);
set(gca,'XColor', 'none','YColor','none')
subplot(4,1,3)
axis([(T/2)/8 inf -2 14])
aux=(abs(phase_f-mean_p(cell_high_locking_3)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[183,212,219]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking_3),:)),'Color','k','linewidth',1.5);
yticks([])
set(gca,'XColor', 'none','YColor','none')
clear aux vals timepoints
subplot(4,1,4)
axis([floor(T/2) inf -2 14])
aux=(abs(phase_f-mean_p(cell_high_locking_2)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[149,125,173]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking_2),:)),'Color','k','linewidth',1.5);
xticks([1,floor(T/2),T]);
xticklabels([0,floor(T/2)/8,T/8]);
set(gca,'YColor','none')
set(gca,'fontsize',18)
xlabel('Time (s)')

%Figures with more than one example - 2
[sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
ind_cell3=find(sorting_ascend==cell_high_locking_3);
ind_cell2=find(sorting_ascend==cell_high_locking_2);
ind_cell1=find(sorting_ascend==cell_high_locking);

figure
subplot(4,1,1)
hold on
for i=1:size(spikes,1)
    %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_descend(i),:),5,'k','filled')
    scatter((1:size(spikes,2))./8,i*spikes(sorting_ascend(i),:),5,'k','filled')
    alpha 0.1
end
scatter((1:size(spikes,2))./8,ind_cell3.*spikes(cell_high_locking_3,:),40,[183,212,219]/255,'filled')
scatter((1:size(spikes,2))./8,ind_cell2.*spikes(cell_high_locking_2,:),40,[149,125,173]/255,'filled')
scatter((1:size(spikes,2))./8,ind_cell1.*spikes(cell_high_locking,:),40,[244,209,181]/255,'filled')
axis([0 inf 2 inf])
ylabel('Neurons #');
yticks([100 400])
set(gca,'XColor', 'none')
set(gca,'fontsize',16);
subplot(4,1,2)
axis([-inf inf -2 18])
ylabel('zscore( DF/F )')
aux=(abs(phase_f-mean_p(cell_high_locking)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[244,209,181]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking),:)),'Color','k','linewidth',1.5);
% set(gca,'XColor', 'none','YColor','none')
set(gca,'XColor', 'none')
yticks([0 8 16])
set(gca,'fontsize', 16);
subplot(4,1,3)
axis([-inf inf -2 18])
aux=(abs(phase_f-mean_p(cell_high_locking_3)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[183,212,219]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking_3),:)),'Color','k','linewidth',1.5);
ylabel('zscore( DF/F )')
yticks([0 8 16])
set(gca,'XColor', 'none');
set(gca,'fontsize', 16);
clear aux vals timepoints
subplot(4,1,4)
axis([-inf inf -2 18])
aux=(abs(phase_f-mean_p(cell_high_locking_2)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[149,125,173]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking_2),:)),'Color','k','linewidth',1.5);
xticks([1,floor(T/2),T]);
xticklabels([0,floor(T/2)/8,T/8]);
ylabel('zscore( DF/F )');
yticks([0 8 16]);
% set(gca,'fontsize', 16);
% set(gca,'XColor', 'none','YColor','none')
% set(gca,'YColor','none')
set(gca,'fontsize',18)
xlabel('Time (s)');

%%%%%%%%%%%%%%%%%%%% Figures with more than one example - 3
[sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
ind_cell3=find(spikes(cell_high_locking_3,:)==1);
ind_cell2=find(spikes(cell_high_locking_2,:)==1);
ind_cell1=find(spikes(cell_high_locking,:)==1);


figure
subplot(4,1,1)
plot((1:length(phase_f))/8,phase_f,'k');
hold on
scatter(ind_cell3./8,mean_p(cell_high_locking_3)*ones(1,length(ind_cell3)),35,[183,212,219]/255,'filled')
scatter(ind_cell2./8,mean_p(cell_high_locking_2)*ones(1,length(ind_cell2)),35,[149,125,173]/255,'filled')
scatter(ind_cell1./8,mean_p(cell_high_locking)*ones(1,length(ind_cell1)),35,[244,209,181]/255,'filled')
set(gca,'XColor', 'none');
ylabel('Phase (rad)');
yticks([-3.14,0,3.14]);
yticklabels([{'-\pi'},{'0'},{'\pi'}]);
set(gca,'fontsize',18)
axis([-inf inf -3.15 3.15])
subplot(4,1,2)
axis([-inf inf -2 18])
aux=(abs(phase_f-mean_p(cell_high_locking)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[244,209,181]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking),:)),'Color','k','linewidth',1.5);
set(gca,'XColor', 'none','YColor','none')
subplot(4,1,3)
axis([-inf inf -2 18])
aux=(abs(phase_f-mean_p(cell_high_locking_3)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
set(gca,'XColor', 'none','YColor','none')
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[183,212,219]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking_3),:)),'Color','k','linewidth',1.5);
yticks([])
set(gca,'XColor', 'none','YColor','none')
clear aux vals timepoints
subplot(4,1,4)
axis([-inf inf -2 18])
aux=(abs(phase_f-mean_p(cell_high_locking_2)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[149,125,173]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking_2),:)),'Color','k','linewidth',1.5);
xticks([1,floor(T/2),T]);
xticklabels([0,floor(T/2)/8,T/8]);
set(gca,'YColor','none')
set(gca,'fontsize',18)
xlabel('Time (s)')



%% Distribution of cells in ensembles (Extended data Fig 10i)

clear Ens
% [~,sorting,~]=get_sorting(spikes_d);
[~,sorting_w,~]=get_sorting_smoothed(spikes_d,dt);

cells_per_ens=floor(N/10);
for e=1:10
    ens(e,:)=sorting_w((e-1)*cells_per_ens+1 : e*cells_per_ens);   
    ense_n(e,:)=e*ones(1,cells_per_ens);
end
Ens(:,1)=reshape(ens,cells_per_ens*10,1);
Ens(:,2)=reshape(ense_n,cells_per_ens*10,1);

if (N-10*cells_per_ens)>0
    Ens(10*cells_per_ens+1:N,1)=sorting_w(10*cells_per_ens+1:N);
    Ens(10*cells_per_ens+1:N,2)=ones(length(sorting_w(10*cells_per_ens+1:N)),1)*10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distribution of ensembles in anatomical cortex

[b,E]=discretize(Ens(:,2),10);
cc=parula(10);
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_anat,'-mat');
pixel_size_new=1.18185;
figure
hold on
for i=1:N
    x=Anat.pos{1,Ens(i,1)}(:,1);
    y=Anat.pos{1,Ens(i,1)}(:,2);
    ov=Anat.overlap{1,Ens(i,1)};    
    if ismember(Ens(i,1),not_locked)==1
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
    else
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(b(i),:),'filled');
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
colormap parula(10)
co=colorbar('XTick',0:1,'XTickLabel',{'1','10'});
co.Label.String = 'Ensemble';

[~,index_ascend]=sort(Ens(:,1),'ascend');
Ens_a=Ens(index_ascend,:);
[b,E]=discretize(Ens_a(:,2),10);
cc=parula(10);

figure
hold on
for i=1:N
    x=Anat.pos{1,Ens_a(i,1)}(:,1);
    y=Anat.pos{1,Ens_a(i,1)}(:,2);
    ov=Anat.overlap{1,Ens_a(i,1)};
    
    if ismember(Ens_a(i,1),not_locked)==1
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
    else
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(Ens_a(i,2),:),'filled');
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
colormap parula(10)
co=colorbar('XTick',0:1,'XTickLabel',{'1','10'});
co.Label.String = 'Ensemble';



%% Distance within and across ensemble, need to come after distribution of cells in ensembles
%Extended dara Fig 10 j,k

pixel_size_new=1.18185;
%Delta position on the tissue
d_within=[];
d_across=[];
c_w=0;%counter within
c_a=0;%counter across
count=0;
for i=1:N
    for j=i+1:N
        count=count+1;
        r_i=[Anat.med{1,i}(1)*pixel_size_new,Anat.med{1,i}(2)*pixel_size_new];
        r_j=[Anat.med{1,j}(1)*pixel_size_new,Anat.med{1,j}(2)*pixel_size_new];
        delta_tissue(i,j)=norm(r_j-r_i);
        delta_tissue(j,i)=norm(r_j-r_i);
        delta_tissue_vec(count)=norm(r_j-r_i);
    end
end

% Delta phases based on the mean phase of each neuron
count=0;
for i=1:N
    for j=i+1:N
        count=count+1;
        delta_phase_mean(i,j)=angdiff(mean_p(i),mean_p(j));
        delta_phase_mean(j,i)=angdiff(mean_p(j),mean_p(i));
        delta_phase_mean_vec(count)=angdiff(mean_p(i),mean_p(j));        
    end
end

% Calculation of distances

for i=1:10 %Loop on ensembles
    within_aux=(find(Ens(:,2)==i));
    within=Ens(within_aux,1);

    across=1:N;
    across(within)=[];
    a_aux_a=0;
    a_auc_w=0;
   
   for l=1:length(within)
       for j=l+1:length(within)
           c_w=c_w+1;
           a_auc_w=a_auc_w+1;
           d_within(c_w)=delta_tissue(within(l),within(j));
           aux_within(a_auc_w)=delta_tissue(within(l),within(j));
       end
   end
   
   for l=1:length(within)
       for j=1:length(across)
           c_a=c_a+1;
           a_aux_a=a_aux_a+1;
           d_across(c_a)=delta_tissue(within(l),across(j));
           aux_across(a_aux_a)=delta_tissue(within(l),across(j));
       end
   end
   
   d_within_ens(i,1)=mean(aux_within);
   d_within_ens(i,2)=std(aux_within);
   d_within_ens(i,3)=std(aux_within)./sqrt(length(aux_within));

   d_across_ens(i,1)=mean(aux_across);
   d_across_ens(i,2)=std(aux_across);
   d_across_ens(i,3)=std(aux_across)./sqrt(length(aux_across));
   
   within_ens{i}=aux_within;
   across_ens{i}=aux_across;

    clear within across aux_within aux_across within_aux
end


%box plot
%I build the matrix for the box plot
mat_within_across=nan(22464,20); %22464 is the maximum number of pairs
count=0;
for i=1:2:20
    count=count+1;
    mat_within_across(1:length(within_ens{count}),i)=within_ens{count}';
    mat_within_across(1:length(across_ens{count}),i+1)=across_ens{count}';

end


figure
boxplot(mat_within_across);
ylabel('Distance (\mum)')
xlabel('Ensemble #')
set(gca,'fontsize',16,'YColor','k','XColor','k');
xticks([1.5:2:19.5])
xticklabels({'1','2','3','4','5','6','7','8','9','10'});
ylim([0 750])
% yticks([])
box off


for i=1:10
    var_w(i,1)=var(within_ens{i});
    var_w(i,2)=var(across_ens{i});
    [h_ens(i),p_ens(i),stat_d{i}]=ranksum(within_ens{i},across_ens{i});
end


%Pooling ensembles
figure
b=bar([1,2],[mean(d_within_ens(:,1)),mean(d_across_ens(:,1))],0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [0.2422    0.1504    0.6603]*1.1;
b.CData(2,:) = [0.0704    0.7457    0.7258];
hold on
er=errorbar([1,2],[mean(d_within_ens(:,1)),mean(d_across_ens(:,1))],[std(d_within_ens(:,1))/sqrt(length(d_within_ens(:,1))),std(d_across_ens(:,1))/sqrt(length(d_across_ens(:,1)))]);
er.LineStyle='none';
er.Color='k';
er.LineWidth=1.5;
axis([0.5,2.5,0,350])
ylabel('Distance on tissue (\mum)')
xticks([1,2])
xticklabels({'Within ensembles', 'Across ensembles'});
box off
set(gca,'fontsize',16)
varw=var(d_within);
vara=var(d_across);

figure
boxplot([d_within_ens(:,1),d_across_ens(:,1)]);
ylabel('Distance (\mum)')
% xlabel('Ensemble #')
set(gca,'fontsize',16,'YColor','k','XColor','k');
xticks([1,2])
xticklabels({'Within ensembles','Across ensembles'});
ylim([0 400])
% yticks([])
box off
[h_allens,p_allens,stat_allens]=ranksum(d_within_ens(:,1),d_across_ens(:,1));

%figure of distance on the tissue agains distance of the phases

% figure
% plot(delta_tissue_vec,delta_phase_TC_vec,'k.');
% ylabel('Delta Preferred phase (rad)');
% xlabel('Delta position [um]');
% set(gca,'fontsize',18);

figure
plot(delta_tissue_vec,delta_phase_mean_vec,'k.');
ylabel('Delta Mean phase (rad)');
xlabel('Delta position [um]');
set(gca,'fontsize',18);


%% Raster plot dropping the cells with highest locking
% Extended dada Fig 7d

%Drop of 80%
p20=prctile(MVL,20);
cells_locked=find(MVL>p20);
cells_nlocked=find(MVL<p20);

[~,sorting_copy,~] = (get_sorting(spikes));
sorting_copy=flip(sorting_copy);
% sorting_copy=flip(sorting_w);
for j=1:length(cells_locked) 
    for i=1:N
        if sorting_copy(i)==cells_locked(j)
           sorting_copy(i)=0;
        end
    end
end
sorting_copy(sorting_copy ==0)=[];
       
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.3]);
hold on
for i=1:length(cells_nlocked) 
    scatter((1:size(spikes,2))./8,i*spikes(sorting_copy(i),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf]);
%title([mouse,' Day',num2str(day)]);
xlabel('Time (s)');
ylabel('Neurons #');
set(gca,'fontsize',18);

clear cells_locked cells_nlocked

%Drop of 60%
p40=prctile(MVL,40);
cells_locked=find(MVL>p40);
cells_nlocked=find(MVL<p40);

[~,sorting_copy,~] = (get_sorting(spikes));
sorting_copy=flip(sorting_copy);
% sorting_copy=flip(sorting_w);
for j=1:length(cells_locked) 
    for i=1:N
        if sorting_copy(i)==cells_locked(j)
           sorting_copy(i)=0;
        end
    end
end
sorting_copy(sorting_copy ==0)=[];
       
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.3]);
hold on
for i=1:length(cells_nlocked) 
    scatter((1:size(spikes,2))./8,i*spikes(sorting_copy(i),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf]);
%title([mouse,' Day',num2str(day)]);
xlabel('Time (s)');
ylabel('Neurons #');
set(gca,'fontsize',18);

clear cells_locked cells_nlocked

%Drop of 40%
p60=prctile(MVL,60);
cells_locked=find(MVL>p60);
cells_nlocked=find(MVL<p60);

[~,sorting_copy,~] = (get_sorting(spikes));
sorting_copy=flip(sorting_copy);
% sorting_copy=flip(sorting_w);
for j=1:length(cells_locked) 
    for i=1:N
        if sorting_copy(i)==cells_locked(j)
           sorting_copy(i)=0;
        end
    end
end
sorting_copy(sorting_copy ==0)=[];
       
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.3]);
hold on
for i=1:length(cells_nlocked) 
    scatter((1:size(spikes,2))./8,i*spikes(sorting_copy(i),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf]);
%title([mouse,' Day',num2str(day)]);
xlabel('Time (s)');
ylabel('Neurons #');
set(gca,'fontsize',18);

clear cells_locked cells_nlocked

%Drop of 20%

p80=prctile(MVL,80);
cells_locked=find(MVL>p80);
cells_nlocked=find(MVL<p80);

% sorting_copy=flip(sorting_w);
[~,sorting_copy,~] = (get_sorting(spikes));
sorting_copy=flip(sorting_copy);
for j=1:length(cells_locked) 
    for i=1:N
        if sorting_copy(i)==cells_locked(j)
           sorting_copy(i)=0;
        end
    end
end

sorting_copy(sorting_copy ==0)=[];
       
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.3]);
hold on
for i=1:length(cells_nlocked) 
    scatter((1:size(spikes,2))./8,i*spikes(sorting_copy(i),:),5,'k','filled')
    alpha 0.3
end
axis([-inf inf 1 inf]);
%title([mouse,' Day',num2str(day)]);
xlabel('Time (s)');
ylabel('Neurons #');
set(gca,'fontsize',18);
                
% figure
% spk_fil=spikes_d(sorting_copy,:);
% spy(spk_fil,'k');
% pbaspect([28 3 8])


%% Figures that I am not using

% table_cells(:,1)=MI_TB;
% table_cells(:,2)=locking;
% table_cells(:,3)=sum(spikes_d,2);
spikes_count_dis=discretize(sum(spikes_d,2),2);

figure
colormap viridis(2)
scatter(MI,locking,55,spikes_count_dis,'filled')
alpha 0.5
ylabel('Locking to phase');
xlabel('MI corrected for bias (bits)');
set(gca,'fontsize',18)
axis([-0.02 0.2 0 1])
box off
colorbar


figure
hold on
scatter(radius_f,MI,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255)
alpha 0.5
ylabel('MI corrected for bias (bits)');
xlabel('Radius PC1-PC2')
set(gca,'fontsize',20)
axis([0 0.2 0 0.2])

figure
hold on
scatter(radius_f,MI,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255)
alpha 0.5
ylabel('MI corrected for bias (bits)');
xlabel('Radius PC1-PC2')
set(gca,'fontsize',20)
axis([0 0.2 0 0.2])
hold on
scatter(radius_f(not_locked),MI(not_locked),40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor','r')
alpha 0.7
ylabel('MI corrected for bias (bits)');
xlabel('Radius PC1-PC2')
set(gca,'fontsize',18)
axis([0 0.2 0 0.2])
box off

%% Participation index (Fig 3d)
[~,sorting_w,~]=get_sorting_smoothed(spikes_d,dt);

%stats

[p,h,stat]=ranksum(nonzeros(PR(:,5)),nonzeros(PR(:,6)));

for i=1:length(locked)
    ind_locked_pr(i)=find(PR_mat(:,1)==locked(i));
end
pr_locked_method2=PR_mat(ind_locked_pr,2);
figure
H=histogram(pr_locked_method2,0:0.05:1);
x=H.BinEdges(2:end);
y=H.Values;
bar(x,y,'FaceColor',[137, 137, 255]./255,'Linewidth',1)
axis([0 1.0 -inf 120]);
ylabel('Counts');
xlabel('PI');
set(gca,'fontsize',18)
box off


figure
H=histogram(PR_mat(:,2),0:0.05:1);
x=H.BinEdges(2:end);
y=H.Values;
bar(x,y,'FaceColor','k','Linewidth',1)
axis([0 1 0 122]);
yticks([0 30 60 90 120])
ylabel('Counts');
xlabel('Participation index');
set(gca,'fontsize',16,'ycolor','k','xcolor','k')
box off



%% Coupling to Population and Sensitivity - Using Tuning curve  % CONDITIONED ON WAVES

N_sh=5; %In the paper we used 200 shuffled iterations

%Downsample and binarize
new_bin = 4;
num_bins=floor(size(spikes_d,2)/new_bin);
for i=1:num_bins %Length of phase conditioned on waves and previously downsampled to 4 bins (0.5 bin size) when computing MI
    sp_do(:,i)=sum(spikes_d(:,(i-1)*new_bin+1:i*new_bin),2);
end

for n=1:N
    sp_do(n,find(sp_do(n,:)>0))=1;
end

%Calculate coactivity
prob_k=rcoact_dist(sp_do);
% figure
% plot(prob_k)

for n=1:N %Loop on neurons
    p_1=find(sp_do(n,:)>0); %Frames in which the neuron fired
    p_0=find(sp_do(n,:)==0); %Frames in which the neuron *did not* fired   
    mean_fr=(length(p_1)/size(sp_do,2)); %Mean firing rate of cell "n"

    disp(n)
    for ense=1:10 %Loop on ensembles
        cells_ens_aux=find(Ens(:,2)==ense); 
        cells_ens=Ens(cells_ens_aux,1); %Cells that belong to ensemble ense
        cells=setdiff(cells_ens,n); %Removes from cells_ens cell *n* for which we are computing the tuning curve
        a=sum(sp_do(cells,:)); %Number of cells that fired at each time point in ensemble "ense"
        prob_ki=rcoact_dist(sp_do(cells,:)); %Probability of coactivity for the cells in ensemble "ense"
        
        %Computing the sensitivity and the tuning curve
        for k=0:N
            aux_k=find(a==k);
            if length(aux_k)>0
                TC2(n,k+1,ense)=length(intersect(p_1,aux_k))./length(aux_k);
                TC3(n,k+1,ense)=(length(intersect(p_1,aux_k))./length(aux_k))/mean_fr;
%                 Arg_sens2(n,k+1,ense)= TC2(n,k+1,ense)* TC2(n,k+1,ense)/(length(p_1)/size(sp_do,2))*prob_ki(k+1);
                Arg_sens2(n,k+1,ense)= TC2(n,k+1,ense)* TC2(n,k+1,ense)*prob_ki(k+1);
            else
                TC2(n,k+1,ense)=0;
                TC3(n,k+1,ense)=0;
                Arg_sens2(n,k+1,ense)=0;                
            end   
            clear aux_k
        end        
%         Sens(n,ense) = sqrt(nansum(Arg_sens2(n,:,ense))- (length(p_1)/size(sp_do,2))*(length(p_1)/size(sp_do,2)));
        Sens(n,ense) = sqrt(nansum(Arg_sens2(n,:,ense))- (mean_fr*mean_fr));
        AUC(n,ense)=sum(TC2(n,:,ense));        
        clear aux2 aux
    end
        
    [c,d]=max(Sens(n,:)); %c is the maximum sensitivity. d is the ensemble that maximizes the sensitivity
    [a,b]=max(AUC(n,:));  %a is the maximum AUC. b is the ensemble that maximizes the AUC
    inferred_ensemble(n)=b; %Inferred ensemble using the AUC of the tuning curve
    inferred_ensemble_sens(n)=d; %Inferred ensemble using the sensitivity
    aux3=find(Ens(:,1)==n);
    real_ensemble(n)=Ens(aux3,2);    
    clear p_1 p_0 cells_ens_aux cells_ens cells a prob_ki
end

%Shuffle spike times of one cell at the time
for n=1:N    
    for sh=1:N_sh
        p_1=randperm(size(sp_do,2),length(find(sp_do(n,:)>0))); %Frames in which the neuron fired
        mean_fr=(length(p_1)/size(sp_do,2)); %Mean firing rate of cell "n"
        disp(n)
        for ense=1:10
            cells_ens_aux=find(Ens(:,2)==ense);
            cells_ens=Ens(cells_ens_aux,1);
            cells=setdiff(cells_ens,n);
            a=sum(sp_do(cells,:));  %Number of cells that fired at each time point in ensemble "ense"
            prob_ki=rcoact_dist(sp_do(cells,:));
            
            %Computing the sensitivity and the tuning curve
            for k=0:N
                aux_k=find(a==k);
                if length(aux_k)>0
                    TC2_sh(n,k+1,ense)=length(intersect(p_1,aux_k))./length(aux_k);
                    TC3_sh(n,k+1,ense)=(length(intersect(p_1,aux_k))./length(aux_k))/mean_fr;
%                     Arg_sens2_sh(n,k+1,ense)= TC2_sh(n,k+1,ense)* TC2_sh(n,k+1,ense)/(length(p_1)/size(sp_do,2))*prob_ki(k+1);
                    Arg_sens2_sh(n,k+1,ense)= TC2_sh(n,k+1,ense)* TC2_sh(n,k+1,ense)*prob_ki(k+1);
                else
                    TC2_sh(n,k+1,ense)=0;
                    TC3_sh(n,k+1,ense)=0;
                    Arg_sens2_sh(n,k+1,ense)=0;
                end
                clear aux_k
            end
            Sens_sh(n,ense,sh) = sqrt(nansum(Arg_sens2_sh(n,:,ense))- (length(p_1)/size(sp_do,2))*(length(p_1)/size(sp_do,2)));
            AUC_sh(n,ense,sh)=sum(TC2_sh(n,:,ense));
            clear aux2 aux
        end
        
        [c,d]=max(Sens_sh(n,:,sh));
        [a,b]=max(AUC_sh(n,:,sh));
        inferred_ensemble_sh(n,sh)=b;
        inferred_ensemble_sens_sh(n,sh)=d;
        
    clear p_1 p_0 cells_ens_aux cells_ens cells a prob_ki
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures 

figure
scatter(inferred_ensemble_sens,real_ensemble,'o','filled');
axis([0 11 0 11])
alpha 0.1
xlabel('Inferred Ensemble - Sens')
ylabel('Real Ensemble')
set(gca,'fontsize',16)
hold on 
axis square
xticks([1 5 10])
yticks([1 5 10])

H=hist3([inferred_ensemble_sens',real_ensemble'],'Edges',{1:1:10 1:1:10});
figure
imagesc(H);
colormap magma
xticks([1 5 10])
yticks([1 5 10])
ylabel('');
xlabel('Real ensemble');
ylabel('Inferred ensemble')
axis square
set(gca,'fontsize',18)
colorbar
caxis([0 40])
co.Label.String = 'Counts';



H=hist3([inferred_ensemble_sens_sh(:,3),real_ensemble']);
figure
imagesc(H);
colormap magma
xticks([1 5 10])
yticks([1 5 10])
ylabel('');
xlabel('Real ensemble');
ylabel('Inferred ensemble')
axis square
set(gca,'fontsize',18)
co=colorbar
caxis([0 40])
co.Label.String = 'Counts';


%% Population coupling - Using Pearson correlation - Extended data Fig 10c


n_sh_pearson=100;
%Downsample spike matrix
new_bin = 4;
num_bins=floor(size(spikes_d,2)/new_bin);
for i=1:num_bins %Length of phase conditioned on waves and previously downsampled to 4 bins (0.5 bin size) when computing MI
    sp_do(:,i)=sum(spikes_d(:,(i-1)*new_bin+1:i*new_bin),2);
end

for n=1:N
    sp_do(n,find(sp_do(n,:)>0))=1;
end

for n=1:N
    disp(n)
    for ense=1:10
        cells_ens_aux=find(Ens(:,2)==ense);
        cells_ens=Ens(cells_ens_aux,1); %Cells that belong to ensemble ense
        cells=setdiff(cells_ens,n); %Removes from cells_ens cell *n* for which we are computing the tuning curve        
        pop_sum=sum(sp_do(cells,:));
        %         res(n,ense)=sum(zdff(n,:).*pop_sum);
        [rho,pval]=corr(sp_do(n,:)',pop_sum');
        res(n,ense)=rho;
        [rho2,pval2]=corrcoef(sp_do(n,:)',pop_sum');
        res2(n,ense)=rho2(1,2);        
        for sh=1:n_sh_pearson
            temp=randperm(length(pop_sum));
            dff_sh=sp_do(n,temp); %shuffle spike train of cell "n"
            res_sh(n,ense,sh)=corr(dff_sh',pop_sum');
        end
    end    
    delta_coupling(n)=max(res(n,:)) - min(res(n,:));
    for sh=1:n_sh_pearson
        delta_coupling_sh(n,sh)=max(res_sh(n,:,sh)) - min(res_sh(n,:,sh));
    end    
    clear  pop_sum cells_ens_aux cells_ens cells
end
%Real ensemble
for n=1:N
    aux3=find(Ens(:,1)==n);
    real_ensemble(n)=Ens(aux3,2);
end
%Inferred ensemble using pearson
for n=1:N
    [val,inf_ensemble(n)]=max(res2(n,:));
%     real_ensemble(n)=Ens(aux3,2);
end
for n=1:N
    for sh=1:n_sh_pearson
        [val,inf_ensemble_sh(n,sh)]=max(res_sh(n,:,sh));
        dist_sh(n,sh)=abs(inf_ensemble_sh(n,sh)-real_ensemble(n));
        if  dist_sh(n,sh)==9  dist_sh(n,sh)=1;
        elseif  dist_sh(n,sh)==8  dist_sh(n,sh)=2;
        elseif  dist_sh(n,sh)==7  dist_sh(n,sh)=3;
        elseif  dist_sh(n,sh)==6  dist_sh(n,sh)=4;
        end
        %     real_ensemble(n)=Ens(aux3,2);
    end
end
figure
scatter(real_ensemble,inf_ensemble)

% Test for significance
thr=prctile(delta_coupling_sh',99);
thr_locked=thr(locked);
delta_coupling_locked=delta_coupling(locked);
prop_locked_to_ensembles=length(find(delta_coupling_locked>thr_locked))/length(locked);
dist=abs(inf_ensemble-real_ensemble);
for n=1:N
   if dist(n)==9 dist(n)=1;
   elseif dist(n)==8 dist(n)=2;
   elseif dist(n)==7 dist(n)=3;
   elseif dist(n)==6 dist(n)=4;
   end       
end

figure
h=histcounts(dist,0:1:6,'Normalization','cdf');
bar(h)
hold on
h_sh=histcounts(dist_sh(:),0:1:6,'Normalization','cdf');
bar(h_sh)
alpha 0.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tuning curves for all cells
figure
imagesc(res(sorting_w,:));
ylabel('Neurons #');
xlabel('Ensemble')
colormap plasma
xticks([1,5,10])
caxis([-0.2 0.5])
set(gca,'fontsize',18)
colorbar
hold on;
line([0,12], [48,48], 'Color', 'k','LineStyle','--','LineWidth',3);
line([0,12], [48*2,48*2], 'Color', 'k','LineStyle','--','LineWidth',3);
line([0,12], [48*3,48*3], 'Color', 'k','LineStyle','--','LineWidth',3);
line([0,12], [48*4,48*4], 'Color', 'k','LineStyle','--','LineWidth',3);
line([0,12], [48*5,48*5], 'Color', 'k','LineStyle','--','LineWidth',3);
line([0,12], [48*6,48*6], 'Color', 'k','LineStyle','--','LineWidth',3);
line([0,12], [48*7,48*7], 'Color', 'k','LineStyle','--','LineWidth',3);
line([0,12], [48*8,48*8], 'Color', 'k','LineStyle','--','LineWidth',3);
line([0,12], [48*9,48*9], 'Color', 'k','LineStyle','--','LineWidth',3);
axis square




