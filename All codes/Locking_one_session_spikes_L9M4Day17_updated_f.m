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
% 
% for w=1:size(table_u,1)
%     %     spikes_r=[spikes_r,spikes_d(:,table_u(w,1):table_u(w,2))];
%     %     phase_r=[phase_r;phase_f(table_u(w,1):table_u(w,2))];
%     
%     if(phase_f(table_u(w,1))<-2 && phase_f(table_u(w,2))>2)
%         cont_min=0;
%         if w==1
%             while table_u(w,1)-cont_min-1>0
%                 if phase_f(table_u(w,1)-cont_min)>-3.12
%                     cont_min=cont_min+1;
%                 else
%                     break;
%                 end
%             end
%         else
%             while table_u(w,1)-cont_min-1>table_u(w-1,2)
%                 if phase_f(table_u(w,1)-cont_min)>-3.12
%                     cont_min=cont_min+1;
%                 else
%                     break;
%                 end
%             end
%         end
%         
%         cont_max=0;
%         if w==size(table_u,1)
%             while table_u(w,2)+cont_max+1<T
%                 if phase_f(table_u(w,2)+cont_max)<3.12
%                     cont_max=cont_max+1;
%                 else
%                     break;
%                 end
%             end
%         else
%             while table_u(w,2)+cont_max+1<table_u(w+1,1)
%                 if phase_f(table_u(w,2)+cont_max)<3.12
%                     cont_max=cont_max+1;
%                 else
%                     break;
%                 end
%             end
%         end
%         
%     elseif phase_f(table_u(w,1))<-2 && phase_f(table_u(w,2))<2
%         
%         cont_min=0;
%         if w==1
%             while table_u(w,1)-cont_min-1>0
%                 if phase_f(table_u(w,1)-cont_min)>-3.12
%                     cont_min=cont_min+1;
%                 else
%                     break;
%                 end
%             end
%         else
%             while table_u(w,1)-cont_min-1>table_u(w-1,2)
%                 if phase_f(table_u(w,1)-cont_min)>-3.12
%                     cont_min=cont_min+1;
%                 else
%                     break;
%                 end
%             end
%         end
%         
%         cont_max=0;
%         
%     elseif phase_f(table_u(w,1))>-2 && phase_f(table_u(w,2))>2
%         
%         cont_min=0;
%         
%         
%         cont_max=0;
%         if w==size(table_u,1)
%             while table_u(w,2)+cont_max+1<T
%                 if phase_f(table_u(w,2)+cont_max)<3.12
%                     cont_max=cont_max+1;
%                 else
%                     break;
%                 end
%             end
%         else
%             while table_u(w,2)+cont_max+1<table_u(w+1,1)
%                 if phase_f(table_u(w,2)+cont_max)<3.12
%                     cont_max=cont_max+1;
%                 else
%                     break;
%                 end
%             end
%         end
%         
%     elseif phase_f(table_u(w,1))>-2 && phase_f(table_u(w,2))<2
%         
%         cont_min=0;
%         cont_max=0;
%     end
%     
%     spikes_r=[spikes_r,spikes_d(:,table_u(w,1)-cont_min:table_u(w,2)+cont_max)];
%     phase_r=[phase_r;phase_f(table_u(w,1)-cont_min:table_u(w,2)+cont_max)];
%     
%     table_u_corrected(w,1)=table_u(w,1)-cont_min;
%     table_u_corrected(w,2)=table_u(w,2)+cont_max;
%     
% end
% 
% table_u_old=table_u;
% table_u=[];
% table_u=table_u_corrected;

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

N_sh=1000;
               
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

%% Computed radius and shuffled values for radius

N_sh=1000;

dt=WS_stat.dt(day,s);
if isinteger(dt)
else
    dt=floor(dt);
end

for i=1:N
    FRp_cor(i,:)=full(fire_rate(spikes_d(i,:),1*dt,'g')); %smooth using as kernel the dt chosen for each session
end

% [coeff,score,latent,tsquared,explained,mu]
% [coefft_cor,scoret_cor,~,~,explained,~] = pca(FRp_cor');
[coefft_cor,scoret_cor,~,~,explained,~] = pca(spikes_d');

% phase_f=(atan2(smooth(scoret(:,2),floor(1*dt)),smooth(scoret(:,1),floor(1*dt))));
radius_f_cor=sqrt(coefft_cor(:,2).*coefft_cor(:,2)+coefft_cor(:,1).*coefft_cor(:,1));
radius_f_cor2=sqrt((explained(2)*coefft_cor(:,2)).*(explained(2)*coefft_cor(:,2))+(explained(1)*coefft_cor(:,1)).*(explained(1)*coefft_cor(:,1)));

for sh=1:N_sh
    
    if( mod(sh,100)==0)
        disp(sh)
    end
    
    spikes_d_sh=shuffle(spikes_d);
    
    %     for i=1:N
    %         FRp_cor_sh(i,:)=full(fire_rate(spikes_d_sh(i,:),1*dt,'g')); %smooth using as kernel the dt chosen for each session
    %     end
    
    %     [coefft_cor_sh,scoret_cor_sh,~] = pca(FRp_cor_sh');
    [coefft_cor_sh,scoret_cor_sh,~] = pca(spikes_d_sh');
    
    %     [coefft_cor_sh,scoret_cor_sh,~,~,explained_cor_sh,~] = pca(FRp_cor_sh');
    % phase_f=(atan2(smooth(scoret(:,2),floor(1*dt)),smooth(scoret(:,1),floor(1*dt))));
    radius_f_cor_sh(:,sh)=sqrt(coefft_cor_sh(:,2).*coefft_cor_sh(:,2)+coefft_cor_sh(:,1).*coefft_cor_sh(:,1));
    
    clear spikes_d_sh FRp_cor_sh coefft_cor_sh scoret_cor_sh explained_cor_sh
    
end

radius_sh_99=prctile(radius_f_cor_sh,99,2);
locked_r=find(radius_f_cor>radius_sh_99);
fraction_locked_r=length(locked_r)/N;
not_locked_r=find(radius_f_cor<=radius_sh_99);
fraction_notlocked_r=length(not_locked_r)./N;

median_radius_data=median(radius_f_cor);
median_radius_data=median(radius_f_cor);



for k=1:length(not_locked_r)
    index_notlocked_r(k)=find(not_locked_r(k)==b);
end

[a,b]=sort(radius_f_cor,'ascend');
figure
hold on
plot(radius_sh_99(b),'.','MarkerSize',15,'Color',[160 160 160]/255);
plot(radius_f_cor(b),'.','MarkerSize',15,'Color','k');
plot(index_notlocked_r,radius_f_cor(not_locked_r),'.','MarkerSize',15,'Color','r');
alpha 0.6
ylabel('Radius');
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

figure
histogram(radius_f_cor_sh(:,1));

figure
histogram(radius_f_cor)

figure
scatter(coefft_cor_sh(:,1),coefft_cor_sh(:,2));

figure
scatter(coefft_cor(:,1),coefft_cor(:,2));


% Alternative shuffling, where I shuffle the positions in the PC1-PC2 plane

for sh=1:N_sh
    
    sh_cells=randperm(N);
    radius_f_cor_sh_cells(:,sh)=sqrt(coefft_cor(sh_cells,2).*coefft_cor(sh_cells,2)+coefft_cor(sh_cells,1).*coefft_cor(sh_cells,1));
    
end


radius_sh_99_cells=prctile(radius_f_cor_sh_cells,99,2);
locked_r_cells=find(radius_f_cor>radius_sh_99_cells);
fraction_locked_r_cells=length(locked_r_cells)/N;
not_locked_r_cells=find(radius_f_cor<=radius_sh_99_cells);
fraction_notlocked_r_cells=length(not_locked_r_cells)./N;

for k=1:length(not_locked_r_cells)
    index_notlocked_r_cells(k)=find(not_locked_r_cells(k)==b);
end

[a,b]=sort(radius_f_cor,'ascend');
figure
hold on
plot(radius_sh_99_cells(b),'.','MarkerSize',15,'Color',[160 160 160]/255);
plot(radius_f_cor(b),'.','MarkerSize',15,'Color','k');
plot(index_notlocked_r_cells,radius_f_cor(not_locked_r_cells),'.','MarkerSize',15,'Color','r');
alpha 0.6
ylabel('Radius');
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

%% Figures 1: Locking and tuning curves

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
                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Radius
figure
scatter(radius_f,MVL,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255)
alpha 0.7
ylabel('Locking to phase');
xlabel('Radius loading PC1- loading PC2');
set(gca,'fontsize',18)
axis([-0.02 0.2 0 1])
box off
[R,P] = corr(radius_f,MVL','Type','Pearson');

figure
% subplot(1,3,1)
scatter(radius_f,MVL,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255)
alpha 0.6
ylabel('Locking to phase');
xlabel('MI_b (bits)');
set(gca,'fontsize',18)
axis([-0.02 0.2 0 1])
box off
hold on
scatter(radius_f(not_locked),MVL(not_locked),40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor','r')
alpha 0.6
ylabel('Locking to phase');
xlabel('Radius loading PCA1 - loading PC2');
set(gca,'fontsize',18)
axis([-0.02 0.2 0 1])
box off


%
% table_cells(:,1)=radius_f;
% table_cells(:,2)=MVL;
% table_cells(:,3)=sum(spikes_d,2);
% spikes_count_dis=discretize(table_cells(:,3),2);
% figure
% colormap viridis(2)
% scatter(radius_f,MVL,55,spikes_count_dis,'filled')
% alpha 0.7
% ylabel('Locking to phase');
% xlabel('Radius PC1-PC2');
% set(gca,'fontsize',18)
% axis([-0.02 0.2 0 1])
% box off
% colorbar

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


[a,b]=max(tuning_curve');
[~,sorting_max_tc]=sort(b,'ascend');
figure
x=1:size(mean_p_locked,2);
x2=[x,fliplr(x)];
inBetween = [mean_p_locked(sorting_max_tc)+std_p_locked(sorting_max_tc), fliplr(mean_p_locked(sorting_max_tc)-std_p_locked(sorting_max_tc))];
fill(x2, inBetween, [252 215 215]./255);
hold on
plot(1:N_n,mean_p_locked(sorting_max_tc)+std_p_locked(sorting_max_tc),'w','linewidth',1.5);
hold on
plot(mean_p_locked(sorting_max_tc)-std_p_locked(sorting_max_tc),'w','linewidth',1.5);
plot(mean_p_locked(sorting_max_tc),'k','linewidth',3);
axis([2 N_n-1 -5 5]);
yticks([-3.14,3.14])
yticklabels({'-\pi','\pi'});
xticks([100 400])
ylabel({'Preferred phase (rad)'})
xlabel('Neuron #')
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
co.Label.String = 'Firing probability';


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
co.Label.String = 'Firing probability';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tuning curves of single cell to the phase using as sorting criterion the peak of the tuning curves



% not_locked=find(locking<=locking_sh_99');
% tuning_curve(not_locked,:)=[];

[a,b]=max(tuning_curve');
[~,sorting_max_tc]=sort(b,'ascend');
figure
imagesc(tuning_curve(sorting_max_tc,:));
xticks([1,20,40])
yticks([100,400])
xticklabels({'-\pi','0','\pi'})
xlabel('Phase (rad)');
ylabel('Neuron #');
colormap inferno
% co=colorbar;
caxis([0 0.15]);
co=colorbar('Ticks',[0,0.15,0.2]);
set(gca,'fontsize',18)
co.Label.String = 'Firing probability';


phase_dis=-pi:2*pi/40:pi;
figure
plot(sort(phase_dis(b),'ascend'),'.k','Markersize',15);
ylabel('Peak of the tuning curve');
xlabel('Neuron #')
set(gca,'fontsize',16);


% tuning_curve_sh(not_locked,:)=[];
[a,bsh]=max(tuning_curve_sh');
[~,sorting_max_tc_sh]=sort(bsh,'ascend');
figure
imagesc(tuning_curve_sh(sorting_max_tc_sh,:));
xticks([1,20,40])
yticks([100,400])
xlabel('Phase (rad)');
ylabel('Neuron #');
colormap inferno
c=colorbar;
caxis([0 0.15]);
co=colorbar('Ticks',[0,0.15]);
set(gca,'fontsize',18)
xticklabels({'\pi','0','\pi'})
co.Label.String = 'Firing probability';

% figure
% plot(sort(phase_dis(bsh),'ascend'),'.k','Markersize',15);
% ylabel('Peak of the tuning curve');
% xlabel('Neuron #')
% set(gca,'fontsize',16);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Examples of individual cells

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

% cell_high_locking=18; %MVL=0.94 mean_p=2.54
figure
aux=(abs(phase_f-mean_p(cell_high_locking)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((cell_high_locking),:)));
plot(zscore(dff((cell_high_locking),:)) - 1.1*min(zscore(dff((cell_high_locking),:))),'Color','k');
box off
yticks([])
axis off
hold on
axis([0 T 0 inf])
title(['MVL = ',num2str(MVL(cell_high_locking))]);
scatter(timepoints,2*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.8
h=area(1:T,(max_val_dff- 1.1*min(zscore(dff((cell_high_locking),:)))).*wave_epochs,'LineStyle','-');
h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
h(1).FaceColor = col_area;
h(1).FaceAlpha = 0.1;
clear aux vals rimepoints lox_local_minima max_val_dff

% cell_high_locking_2=30; %MVL=0.91 mean_p=-1.8139
figure
aux=(abs(phase_f-mean_p(cell_high_locking_2)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((cell_high_locking_2),:)));
plot(zscore(dff((cell_high_locking_2),:)) - 1.1*min(zscore(dff((cell_high_locking_2),:))),'Color','k');
box off
yticks([])
axis off
hold on
axis([0 T 0 inf])
title(['MVL = ',num2str(MVL(cell_high_locking_2))]);
scatter(timepoints,2*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.9
h=area(1:T,(max_val_dff- 1.1*min(zscore(dff((cell_high_locking_2),:)))).*wave_epochs,'LineStyle','-');
h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
h(1).FaceColor = col_area;
h(1).FaceAlpha = 0.1;
hold on
clear aux vals rimepoints lox_local_minima max_val_dff

% cell_high_locking_3=474; %MVL=0.90 mean_p=-0.003
figure
aux=(abs(phase_f-mean_p(cell_high_locking_3)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((cell_high_locking_3),:)));
plot(zscore(dff((cell_high_locking_3),:)) - 1.1*min(zscore(dff((cell_high_locking_3),:))),'Color','k');
box off
yticks([])
axis off
hold on
axis([0 T 0 inf])
title(['MVL = ',num2str(MVL(cell_high_locking_3))]);
scatter(timepoints,2.5*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.9
h=area(1:T,(max_val_dff- 1.1*min(zscore(dff((cell_high_locking_3),:)))).*wave_epochs,'LineStyle','-');
h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
h(1).FaceColor = col_area;
h(1).FaceAlpha = 0.1;
hold on
clear aux vals rimepoints lox_local_minima max_val_dff

% cell_med_locking=97; %MVL=0.62
figure
plot(zscore(dff((cell_med_locking),:))- 1.1*min(zscore(dff((cell_med_locking),:))),'Color','k');
box off
axis([-inf inf -inf inf])
yticks([])
axis off
hold on
aux=(abs(phase_f-mean_p(cell_med_locking)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((cell_med_locking),:)));
scatter(timepoints,2*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.9
h=area(1:T,(max_val_dff- 1.1*min(zscore(dff((cell_high_locking_2),:)))).*wave_epochs,'LineStyle','-');
h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
h(1).FaceColor = col_area;
h(1).FaceAlpha = 0.1;
axis([0 T -inf inf])
clear aux vals rimepoints lox_local_minima max_val_dff
title(['MVL = ',num2str(MVL(cell_med_locking))]);

% cell_low_locking_not_locked=43; %MVL=0.1
figure
plot(zscore(dff((cell_low_locking_not_locked),:))- 1.1*min(zscore(dff((cell_low_locking_not_locked),:))),'Color','k');
box 
axis([-inf inf -inf inf])
yticks([])
axis off
hold on
aux=(abs(phase_f-mean_p(cell_low_locking_not_locked)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((cell_low_locking_not_locked),:)));
scatter(timepoints,2*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.9
h=area(1:T,(max_val_dff- 1.1*min(zscore(dff((cell_high_locking_2),:)))).*wave_epochs,'LineStyle','-');
h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
h(1).FaceColor = col_area;
h(1).FaceAlpha = 0.1;
axis([0 T -inf inf])
clear aux vals rimepoints lox_local_minima max_val_dff
axis([0 T -inf inf])
clear aux vals rimepoints lox_local_minima max_val_dff
title(['MVL = ',num2str(MVL(cell_low_locking_not_locked))]);

% cell_low_locking=217; %MVL=0.4
figure
plot(zscore(dff((cell_low_locking),:))- 1.1*min(zscore(dff((cell_low_locking),:))),'Color','k');
box 
axis([-inf inf -inf inf])
yticks([])
axis off
hold on
aux=(abs(phase_f-mean_p(cell_low_locking)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((cell_low_locking),:)));
scatter(timepoints,2*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.9
h=area(1:T,(max_val_dff- 1.1*min(zscore(dff((cell_high_locking_2),:)))).*wave_epochs,'LineStyle','-');
h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
h(1).FaceColor = col_area;
h(1).FaceAlpha = 0.1;
axis([0 T -inf inf]);
clear aux vals rimepoints lox_local_minima max_val_dff
title(['MVL = ',num2str(MVL(cell_low_locking))]);

% cell_low_locking_not_locked2=449; %MVL=0.15
figure
plot(zscore(dff((cell_low_locking_not_locked2),:))- 1.1*min(zscore(dff((cell_low_locking_not_locked2),:))),'Color','k');
box 
axis([-inf inf -inf inf])
yticks([])
axis off
hold on
aux=(abs(phase_f-mean_p(cell_low_locking_not_locked2)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((cell_low_locking_not_locked2),:)));
scatter(timepoints,2*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.9
h=area(1:T,(max_val_dff- 1.1*min(zscore(dff((cell_high_locking_2),:)))).*wave_epochs,'LineStyle','-');
h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
h(1).FaceColor = col_area;
h(1).FaceAlpha = 0.1;
axis([0 T -inf inf]);
clear aux vals rimepoints lox_local_minima max_val_dff
title(['MVL = ',num2str(MVL(cell_low_locking_not_locked2))]);


%cell_med_locking_2=86; %MVL=0.59
figure
plot(zscore(dff((cell_med_locking_2),:))- 1.1*min(zscore(dff((cell_med_locking_2),:))),'Color','k');
box 
axis([-inf inf -inf inf])
yticks([])
axis off
hold on
aux=(abs(phase_f-mean_p(cell_med_locking_2)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((cell_med_locking_2),:)));
scatter(timepoints,3*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.9
h=area(1:T,(max_val_dff- 1.6*min(zscore(dff((cell_high_locking_2),:)))).*wave_epochs,'LineStyle','-');
h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
h(1).FaceColor = col_area;
h(1).FaceAlpha = 0.1;
axis([0 T -inf inf]);
clear aux vals rimepoints lox_local_minima max_val_dff
title(['MVL = ',num2str(MVL(cell_med_locking_2))]);

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

%%%%%%%%%%%%%%%%%%%% Figures with more than one example - 4
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
axis([floor(T/2)/8 T/8 -3.15 3.15])
yticks([-3.14,0,3.14]);
yticklabels([{'-\pi'},{'0'},{'\pi'}]);
set(gca,'fontsize',18)
subplot(4,1,2)
axis([floor(T/2) inf -2 14])
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
axis([floor(T/2) inf -2 14])
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


%%%%%%%%%%%%%%%%%%%% Figures with more than one example - 5
figure
subplot(2,1,1)
plot((1:length(phase_f))/8,phase_f,'k');
hold on
scatter(ind_cell3./8,mean_p(cell_high_locking_3)*ones(1,length(ind_cell3)),35,[183,212,219]/255,'filled')
scatter(ind_cell2./8,mean_p(cell_high_locking_2)*ones(1,length(ind_cell2)),35,[149,125,173]/255,'filled')
scatter(ind_cell1./8,mean_p(cell_high_locking)*ones(1,length(ind_cell1)),35,[244,209,181]/255,'filled')
set(gca,'XColor', 'none');
ylabel('Phase (rad)');
axis([floor(T/2)/8 T/8 -3.15 3.15])
yticks([-3.14,0,3.14]);
yticklabels([{'-\pi'},{'0'},{'\pi'}]);
set(gca,'fontsize',18)
subplot(2,1,2)
aux=(abs(phase_f-mean_p(cell_high_locking)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,14],'--','color',[244,209,181]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking),:)),'Color',[255 128 0]/255,'linewidth',1.5);
clear aux vals timepoints low_local_minima
aux=(abs(phase_f-mean_p(cell_high_locking_3)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,14],'--','color',[183,212,219]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking_3),:)),'Color',[0 0 204]/255,'linewidth',1.5);
clear aux vals timepoints
axis([floor(T/2) inf -2 14])
aux=(abs(phase_f-mean_p(cell_high_locking_2)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,14],'--','color',[149,125,173]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking_2),:)),'Color',[153 51 255]/255,'linewidth',1.5);
xticks([1,floor(T/2),T]);
xticklabels([0,floor(T/2)/8,T/8]);
set(gca,'YColor','none')
set(gca,'fontsize',18)
xlabel('Time (s)')





%%%%%%%%%%%%%%%%%%%% Figures with more than one example - 6
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
subplot(2,1,2)
aux=(abs(phase_f-mean_p(cell_high_locking)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,14],'--','color',[244,209,181]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking),:)),'Color',[255 128 0]/255,'linewidth',1.5);
clear aux vals timepoints low_local_minima
aux=(abs(phase_f-mean_p(cell_high_locking_3)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,14],'--','color',[183,212,219]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking_3),:)),'Color',[0 0 204]/255,'linewidth',1.5);
clear aux vals timepoints
axis([floor(T/2) inf -2 14])
aux=(abs(phase_f-mean_p(cell_high_locking_2)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,14],'--','color',[149,125,173]/255,'linewidth',2);
end
plot(zscore(dff((cell_high_locking_2),:)),'Color',[153 51 255]/255,'linewidth',1.5);
xticks([1,floor(T/2),T]);
xticklabels([0,floor(T/2)/8,T/8]);
set(gca,'YColor','none')
set(gca,'fontsize',18)
xlabel('Time (s)')


%% Distribution of cells in ensembles

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

% Delta phases based on the phase that maximizes the tuning curve
% [a,b]=max(prob_phase_firing');
% pref_phase=phase_TC(b);
% 
% count=0;
% for i=1:N
%     for j=i+1:N
%         count=count+1;
%         delta_phase_TC(i,j)=angdiff(pref_phase(i),pref_phase(j));
%         delta_phase_TC(j,i)=angdiff(pref_phase(j),pref_phase(i));
%         delta_phase_TC_vec(count)=angdiff(pref_phase(i),pref_phase(j));
%         
%     end
% end


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


figure
b=bar([d_within_ens(:,1),d_across_ens(:,1)],0.8,'FaceColor','flat');
b(1).CData = [0.2422    0.1504    0.6603]*1.1;
b(2).CData = [0.0704    0.7457    0.7258];
% b(2).CData = [ 0.9686    0.7216    0.1373];
hold on
y=[d_within_ens(:,1),d_across_ens(:,1)];
err=[d_within_ens(:,3),d_across_ens(:,3)];
ngroups = 10;
nbars = 2;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), err(:,i), 'k.','linewidth',1.5);
end
hold off
axis([0.5,10.5,0,400])
ylabel('Distance on tissue (\mum)')
xlabel('Ensemble #');
box off
set(gca,'fontsize',16)
legend({'Within ensemble', 'Across ensemble'})
legend boxoff

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


% anova1(within_ens{i})

%Pooling distances
% figure
% b=bar([1,2],[mean(d_within),mean(d_across)],0.5);
% b.FaceColor = 'flat';
% b.CData(1,:) = [0.2422    0.1504    0.6603]*1.1;
% b.CData(2,:) = [0.0704    0.7457    0.7258];
% hold on
% er=errorbar([1,2],[mean(d_within),mean(d_across)],[std(d_within)/sqrt(length(d_within)),std(d_across)/sqrt(length(d_across))]);
% er.LineStyle='none';
% er.Color='k';
% er.LineWidth=1.5;
% axis([0.5,2.5,0,350])
% ylabel('Distance on tissue (\mum)')
% xticks([1,2])
% xticklabels({'Within ensembles', 'Across ensembles'});
% box off
% set(gca,'fontsize',16)
% varw=var(d_within);
% vara=var(d_across);
% [h_all,p_all]=ttest2(d_within,d_across);

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
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
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
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
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
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
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
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
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

    

%% Distribution of MI for each ensemble - FIGURES

%Now I take the 30% of most informative cells and I see how they are
%distributed among ensambles
% MI_t=prctile(MI_TB,70);
% cells_hmi=find(MI>=MI_t);
% radius_t=prctile(radius,70);
% cells_hr=find(radius>=radius_t);
locking_t=prctile(MVL,70);
cells_klocking=find(MVL>=locking_t);
% 
% ens_hmi=[];
% count=0;
% for i=1:length(cells_hmi)
%     count=count+1;
%     aux=find(Ens(:,1)==cells_hmi(i));
%     ens_hmi(count)=Ens(aux,2);
% end
% y=histogram(ens_hmi,1:11);
% figure
% bar(y.BinEdges(1:end-1),y.Values,'facecolor',[0.4667    0.6745    0.1882]);
% ylabel({'Counts';'30% with largest MI'});
% xlabel('Ensemble #');
% set(gca,'fontsize',16)
% box off
% 
% figure
% ens_hr=[];
% count=0;
% for i=1:length(cells_hr)
%     count=count+1;
%     aux=find(Ens(:,1)==cells_hr(i));
%     ens_hr(count)=Ens(aux,2);
% end
% y=histogram(ens_hr,1:11);
% figure
% bar(y.BinEdges(1:end-1),y.Values,'facecolor',[0.4667    0.6745    0.1882]);
% ylabel({'Counts';'30% with largest Radius'});
% xlabel('Ensemble #');
% set(gca,'fontsize',16)
% box off

figure
ens_hr=[];
count=0;
for i=1:length(cells_klocking)
    count=count+1;
    aux=find(Ens(:,1)==cells_klocking(i));
    ens_hlocking(count)=Ens(aux,2);
end
y=histogram(ens_hlocking,1:11);
figure
bar(y.BinEdges(1:end-1),y.Values,'facecolor',[0.4667    0.6745    0.1882]);
ylabel({'Counts';'30% of cells with largest locking'});
xlabel('Ensemble #');
set(gca,'fontsize',16)
box off

 

%% Compute MI using deconvolved spikes

% Number of bins for discretization

disc_phase=10; %Number of bins for the phase
% disc_spk=4; %Number of bins for the spikes

% Downsampling of spikes and phase
new_bin = 4;

phase_down=downsample(phase,new_bin);
% phase_di=discretize(phase_down,-pi:2*pi/disc_phase:pi);

for i=1:length(phase_down)
    spk_do(:,i)=sum(spikes_d(:,(i-1)*new_bin+1:i*new_bin),2);
end

%Calculate information
cells=1:N;
N_sh=1000;
for ind=1:N
    disp(ind)
    i=cells(ind);
    
    %Preprocess signal
      
    calc=spk_do(ind,:);
%     calc_di=discretize(calc,disc_spk);
    
    table(:,1)=calc';
    table(:,2)=phase_down;
                
    edges_spk=0:1:max(calc)+1;
    edges_ph=-pi:2*pi/disc_phase:pi;%0:max(full_speed_wave)/bins_speed:max(full_speed_wave);
        
    MI_b=compute_MI_SGC(table,edges_spk,edges_ph);
       
    for sh=1:N_sh
        table_sh(:,1)=table(:,1);
        table_sh(:,2)=table(randperm(size(table,1)),2);
        MI_sh(sh)=compute_MI_SGC(table_sh,edges_spk,edges_ph);
        clear table_sh
    end
    MI(ind)= MI_b - mean(MI_sh);
    MI_withbias(ind)=MI_b;
    bias(ind)=mean(MI_sh);
    
    clear table calc MI_sh sig
end

cells_inf=find(MI>0);
cells_no_inf=find(MI<=0);
% fraction_cells_inf=length(cells_mi_calc)./N;
% fraction_cells_no_inf=length(cells_mi_calc)./length(cells);

% non_inf=find(p>0.001);
% infor=find(p<=0.001);
% figure
% labels = {'Inf','Non Inf'};
% pie([length(infor)/N,length(non_inf)/N],labels)
% colormap jet(4)
% set(gca,'fontsize',16)


for k=1:length(cells_no_inf)
    index_no_inf(k)=find(cells_no_inf(k)==b);
end

[a,b]=sort(MI_withbias,'ascend');
figure
hold on
plot(bias(b),'.','MarkerSize',15,'Color',[160 160 160]/255);
plot(MI_withbias(b),'.','MarkerSize',15,'Color','k');
% plot(index_no_inf,MI_withbias(cells_no_inf),'.','MarkerSize',15,'Color','r');
alpha 0.6
ylabel('MI (bits)');
xlabel('Neuron #')
set(gca,'fontsize',18)
box off
legend('<MI_S_h_u_f_f_l_e>', 'Data')
legend boxoff
% axis([0 490 0 0.2])
xticks([100 400])
yticks([0 0.1 0.2])
axis([-inf 490 0 0.2])
axis square


figure
% subplot(1,3,1)
scatter(MI,locking,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255)
alpha 0.6
ylabel('Locking to phase');
xlabel('MI corrected for bias (bits)');
axis([-0.02 0.2 0 1])
box off
hold on
scatter(MI(not_locked),locking(not_locked),40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor','r')
alpha 0.7
ylabel('Locking to phase');
xlabel('MI corrected for bias (bits)');
set(gca,'fontsize',18)
axis([-0.02 0.2 0 1])
box off
[R,P] = corr(MI',locking','Type','Pearson');


figure
% subplot(1,3,1)
scatter(MI,locking,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255)
alpha 0.6
ylabel('Locking to phase');
xlabel('MI_b (bits)');
set(gca,'fontsize',18)
axis([-0.02 0.2 0 1])
box off
hold on
scatter(MI(cells_no_inf),locking(cells_no_inf),40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor','g')
alpha 0.6
ylabel('Locking to phase');
xlabel('MI corrected for bias (bits)');
set(gca,'fontsize',18)
axis([-0.02 0.2 0 1])
box off


figure
hold on
scatter(radius_f,MI,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255)
alpha 0.5
ylabel('MI corrected for bias (bits)');
xlabel('Radius PC1-PC2')
set(gca,'fontsize',20)
axis([0 0.2 0 0.2])
[R,P] = corr(MI',radius_f,'Type','Pearson');

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

%% Participation index
[~,sorting_w,~]=get_sorting_smoothed(spikes_d,dt);

% [table_u,N,T]=identify_waves_latestversion_5(mouse,day,num_clus_discr,dt,make_fig,spikes_d);      
% [~,sorting_angle,~]=get_sorting(spikes_d);

% Recurrent algorithms for participation ratio
spikes_sorted=spikes(sorting_w,:);
new_mat=spikes_sorted;

num_steps=6;
step=0;
mean_PR(1)=0;
std_PR(1)=0.5;
mean_PR_I(1)=0;
std_PR_I(1)=0.5;
n_waves=size(table_u,1);
% for i=1:N
%     new_mat(i,:)=spikes_d(sorting_angle(i),:);
% end
    
cell_wave_part=zeros(N,size(table_u,1));

for j=1:num_steps

    step=step+1;
    for i=1:size(new_mat,1) %Loop on elements of mat
        calc=(new_mat(i,:)); %spike train of row i

        for w=1:size(table_u,1) %Spikes of each cell (or row of new_mat) per wave
            spikes_per_wave(i,w)=sum(calc(table_u(w,1):table_u(w,2)));
        end
        
        aux=spikes_per_wave(i,:);
        aux2=find((cumsum(sort(aux,'descend'))./sum(aux))>0.90,1);%Number of waves needed to account for 95% of the spikes
        
        if j==1 %This indicates the waves in which a neuron fired, based on the criterion of the waves for which 95% of spiking is captured
            [ind1,ind2]=sort(aux,'descend');
            cell_wave_part(i,ind2(1:aux2))=1;
            clear ind1 ind2
        end
          
        if isempty(aux2)==1
            aux2=0;
        end
        
        PR(i,step)=aux2/size(table_u,1);
        
        %Second measure of PI        
        aux3=cumsum(sort(aux,'descend'))./sum(aux);        
        auc=sum(aux3*1/n_waves);        
%         PR_I(i,step)=auc-0.5;
        PR_I(i,step)=1-auc;

        clear aux aux2 aux3 auc

    end
    
    mean_PR(step)=mean(nonzeros(PR(:,step)),1);
    std_PR(step)=std(nonzeros(PR(:,step)),[],1);
    sem_PR(step)=std(nonzeros(PR(:,step)),[],1)/sqrt(length(nonzeros(PR(:,step))));
    
    mean_PR_I(step)=mean(nonzeros(PR_I(:,step)),1);
    std_PR_I(step)=std(nonzeros(PR_I(:,step)),[],1);
    
    new_mat2=new_mat;
    clear new_mat;
    new_mat=coarse_graining(new_mat2);
    
    clear spikes_per_wave aux aux2
end


%PR is ordered according to sorting_w. The rest of the quantities, such as
%MVL, is ordered in a "canonical" basis. In order to keep the PI of the
%locked cells, we need to identify the position of the locked cells in
%sorting_w, and keep only those. These positions will be given by
%index_sorting_locked. Sorting_w_locked is the sorting_w vector without the
%not_locked cells. The connection between these two is the following: If we
%take sorting_w, and keep only the components given by
%index_sorting_locked, we should get sorting_w_locked.

PR_mat(:,2)=PR(:,1);
PR_mat(:,1)=sorting_w;

count=0;
sorting_copy=sorting_w;
for j=1:length(not_locked)
    for i=1:N
        if sorting_copy(i)==not_locked(j)
            sorting_copy(i)=0;
            count=count+1;
            index_sorting_not_locked(count)=i; %Components of the sorted vector that are *not* locked
        end
    end
end
index_sorting_locked=1:N;
index_sorting_locked(index_sorting_not_locked)=[]; %Components of the sorted vector that are locked

sorting_w_not_locked=sorting_w(find(sorting_copy ==0)); %sorting_w without cells that are  locked

sorting_copy(sorting_copy ==0)=[];
sorting_w_locked=sorting_copy; %sorting_w without cells that are *not* locked

PR_locked=PR(:,1);
PR_I_locked=PR_I(:,1);
PR_not_locked=PR(:,1);
PR_I_not_locked=PR_I(:,1);

PR_locked=PR_locked(index_sorting_locked);
PR_I_locked=PR_I_locked(index_sorting_locked);
PR_not_locked=PR_not_locked(index_sorting_not_locked);
PR_I_not_locked=PR_I_not_locked(index_sorting_not_locked);


figure
errorbar(1:num_steps,mean_PR,sem_PR,'k-o','Linewidth',2);
ylabel('PI');
xlabel('Merging steps');
axis([0.5 6.5 0 1])
xticks([1 2 3 4 5 6])
set(gca,'fontsize',18)
box off

figure
errorbar(1:num_steps,mean_PR,std_PR,'k-o','Linewidth',2);
ylabel('PI');
xlabel('Merging steps');
axis([0.5 6.5 0 1])
xticks([1 2 3 4 5 6])
set(gca,'fontsize',18)
box off

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of PI

% tuning_curve=prob_phase_firing;
% not_locked=find(locking<=locking_sh_99');
% [~,b2]=sort(mean_p,'ascend');
[b,E]=discretize(PR(:,1),40);
cc= winter(40);
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
load(file_name_anat,'-mat');
pixel_size_new=1.18185;

figure
hold on
for i=1:N
    x=Anat.pos{1,i}(:,1);
    y=Anat.pos{1,i}(:,2);
    ov=Anat.overlap{1,i};
    
    index=find(sorting_w==i);
    if ismember(i,not_locked)==1
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,'r','filled');
    else
        scatter(x(find(ov<2))*pixel_size_new,y(find(ov<2))*pixel_size_new,1,cc(b(index),:),'filled');
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
colormap  winter(40)
co=colorbar('XTick',0:1,'XTickLabel',{'0.24','0.76'});
co.Label.String = 'PI';

clear b E

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures PI that I am not using

figure
scatter(MI(sorting_w_locked),PR_locked,50,'o','filled','k')
alpha 0.5
xlabel('MI (bits)')
ylabel('PI')
hold on
scatter(MI(sorting_w_not_locked),PR_not_locked,50,'o','filled','r') 
set(gca,'fontsize',16)

figure
scatter(MVL(sorting_w_locked),PR_locked,50,'o','filled','k')
alpha 0.5
xlabel('MVL')
ylabel('PI')
hold on
scatter(MVL(sorting_w_not_locked),PR_not_locked,50,'o','filled','r') 
set(gca,'fontsize',16)

figure
imagesc(cell_wave_part(sorting_w,:))
colormap magma(2)
set(gca,'fontsize',16)
ylabel('Neuron #');
xlabel('Wave #');
colorbar








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures PI_I

figure
scatter(MI(sorting_w_locked),PR_I_locked,50,'o','filled','k')
alpha 0.5
xlabel('MI (bits)')
ylabel('PI_A_r_e_a')
hold on
scatter(MI(sorting_w_not_locked),PR_I_not_locked,50,'o','filled','r') %CHECK
set(gca,'fontsize',16)

figure
scatter(MVL(sorting_w_locked),PR_I_locked,50,'o','filled','k')
alpha 0.5
xlabel('MVL')
ylabel('PI_A_r_e_a')
hold on
scatter(MVL(sorting_w_not_locked),PR_I_not_locked,50,'o','filled','r') %CHECK
set(gca,'fontsize',16)

figure
H=histogram(PR_I_locked,0:0.05:1);
x=H.BinEdges(2:end);
y=H.Values;
bar(x,y,'FaceColor',[137, 137, 255]./255,'Linewidth',1)
axis([0 1.0 -inf 100]);
ylabel('Counts');
ylabel('PI_A_r_e_a')
set(gca,'fontsize',18)
box off

figure
errorbar(1:num_steps,mean_PR_I,std_PR_I,'k-o','Linewidth',2);
ylabel('PI_A_r_e_a')
xlabel('Merging steps');
axis([0.5 6.5 0 1])
xticks([1 2 3 4 5 6])
set(gca,'fontsize',18)
box off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures of examples of cells with PI
aux=PR(:,1);

max_cell=412;
min_cell=97;
med_cell=271;

max_cell_s=find(sorting_w==412);
min_cell_s=find(sorting_w==97);
med_cell_s=find(sorting_w==271);

for w=1:size(table_u,1)
    spikes_per_wave_min(w)=sum(spikes(min_cell,table_u(w,1):table_u(w,2)));
end
for w=1:size(table_u,1)
    spikes_per_wave_max(w)=sum(spikes(max_cell,table_u(w,1):table_u(w,2)));
end
for w=1:size(table_u,1)
    spikes_per_wave_med(w)=sum(spikes(med_cell,table_u(w,1):table_u(w,2)));
end

cumsum_min=cumsum(sort(spikes_per_wave_min,'descend')./sum(spikes_per_wave_min));
cumsum_max=cumsum(sort(spikes_per_wave_max,'descend')./sum(spikes_per_wave_max));
cumsum_med=cumsum(sort(spikes_per_wave_med,'descend')./sum(spikes_per_wave_med));

cc=plasma(4);
figure
plot(cumsum_min,'-','Linewidth',2.5,'Color',cc(1,:));
hold on
plot(cumsum_med,'-','Linewidth',2.5,'Color',cc(2,:));
plot(cumsum_max,'-','Linewidth',2.5,'Color',cc(3,:));
ylabel({'Normalized cumulative';'sum of spikes'});
xlabel('Wave #');
set(gca,'fontsize',18)
axis([0.5 26 0 1.1])
xticks([1 25])
box off
legend({'PI=0.32';'PI=0.56';'PI=0.76'})
legend boxoff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures of examples of cells with PI
aux=PR(:,1);

max_cell=10;
min_cell=280;
med_cell=74;

max_cell_s=find(sorting_w==max_cell);
min_cell_s=find(sorting_w==min_cell);
med_cell_s=find(sorting_w==med_cell);

figure
aux=(abs(phase_f-mean_p(min_cell_s)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((min_cell_s),:)));
plot(zscore(dff((min_cell_s),:)) - 1.1*min(zscore(dff((min_cell_s),:))),'Color','k');
box off
yticks([])
axis off
hold on
axis([0 T 0 inf])
title(['MVL = ',num2str(MVL(min_cell_s)),' - PI = 0.24']);
scatter(timepoints,2*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.8
% h=area(1:T,(max_val_dff- 1.1*min(zscore(dff((max_cell_s),:)))).*wave_epochs,'LineStyle','-');
% h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
% h(1).FaceColor = col_area;
% h(1).FaceAlpha = 0.1;
clear aux vals rimepoints lox_local_minima max_val_dff

figure
aux=(abs(phase_f-mean_p(max_cell_s)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((max_cell_s),:)));
plot(zscore(dff((max_cell_s),:)) - 1.1*min(zscore(dff((max_cell_s),:))),'Color','k');
box off
yticks([])
axis off
hold on
axis([0 T 0 inf])
title(['MVL = ',num2str(MVL(max_cell_s)),' - PI = 0.76']);
scatter(timepoints,2*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.8
% h=area(1:T,(max_val_dff- 1.1*min(zscore(dff((max_cell_s),:)))).*wave_epochs,'LineStyle','-');
% h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
% h(1).FaceColor = col_area;
% h(1).FaceAlpha = 0.1;
clear aux vals rimepoints lox_local_minima max_val_dff

figure
aux=(abs(phase_f-mean_p(med_cell_s)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((med_cell_s),:)));
plot(zscore(dff((med_cell_s),:)) - 1.1*min(zscore(dff((med_cell_s),:))),'Color','k');
box off
yticks([])
axis off
hold on
axis([0 T 0 inf])
title(['MVL = ',num2str(MVL(med_cell_s)),' - PI = 0.56']);
scatter(timepoints,2*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.8
% h=area(1:T,(max_val_dff- 1.1*min(zscore(dff((max_cell_s),:)))).*wave_epochs,'LineStyle','-');
% h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
% h(1).FaceColor = col_area;
% h(1).FaceAlpha = 0.1;
clear aux vals rimepoints lox_local_minima max_val_dff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures raster + dff

aux=PR(:,1);

max_cell=3; %0.72
med_cell=26; %0.56
min_cell=258; %0.28

max_cell_s=find(sorting_w==max_cell);
min_cell_s=find(sorting_w==min_cell);
med_cell_s=find(sorting_w==med_cell);

for w=1:size(table_u,1)
    spikes_per_wave_min(w)=sum(spikes(min_cell_s,table_u(w,1):table_u(w,2)));
end
for w=1:size(table_u,1)
    spikes_per_wave_max(w)=sum(spikes(max_cell_s,table_u(w,1):table_u(w,2)));
end
for w=1:size(table_u,1)
    spikes_per_wave_med(w)=sum(spikes(med_cell_s,table_u(w,1):table_u(w,2)));
end

cumsum_min=cumsum(sort(spikes_per_wave_min,'descend')./sum(spikes_per_wave_min));
cumsum_max=cumsum(sort(spikes_per_wave_max,'descend')./sum(spikes_per_wave_max));
cumsum_med=cumsum(sort(spikes_per_wave_med,'descend')./sum(spikes_per_wave_med));

cc=plasma(4);
figure
plot(cumsum_min,'-','Linewidth',2.5,'Color',cc(1,:));
hold on
plot(cumsum_med,'-','Linewidth',2.5,'Color',cc(2,:));
plot(cumsum_max,'-','Linewidth',2.5,'Color',cc(3,:));
ylabel({'Normalized cumulative';'sum of spikes'});
xlabel('Wave #');
set(gca,'fontsize',18)
axis([0.5 26 0 1.1])
xticks([1 25])
box off
legend({strcat('PI=',num2str(aux(min_cell)));strcat('PI=',num2str(aux(med_cell)));strcat('PI=',num2str(aux(max_cell)))})
legend boxoff



aux=PR(:,1);
max_cell=10;
min_cell=280;
med_cell=74;

max_cell_s=find(sorting_w==max_cell);
min_cell_s=find(sorting_w==min_cell);
med_cell_s=find(sorting_w==med_cell);

[sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
ind_cell3=find(sorting_ascend==max_cell_s);
ind_cell2=find(sorting_ascend==med_cell_s);
ind_cell1=find(sorting_ascend==min_cell_s);
% 
% ind_cell3=find(sorting_w==cell_high_locking_3);
% ind_cell2=find(sorting_w==cell_high_locking_2);
% ind_cell1=find(sorting_w==cell_high_locking);

figure
subplot(4,1,1)
hold on
for i=1:size(spikes,1)
    %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_descend(i),:),5,'k','filled')
    scatter((1:size(spikes,2))./8,i*spikes(sorting_ascend(i),:),5,'k','filled')
    alpha 0.1
end
scatter((1:size(spikes,2))./8,ind_cell3.*spikes(max_cell_s,:),40,[183,212,219]/255,'filled')
scatter((1:size(spikes,2))./8,ind_cell2.*spikes(med_cell_s,:),40,[149,125,173]/255,'filled')
scatter((1:size(spikes,2))./8,ind_cell1.*spikes(min_cell_s,:),40,[244,209,181]/255,'filled')
axis([0 inf 2 inf])
ylabel('Neurons #');
yticks([100 400])
set(gca,'XColor', 'none')
set(gca,'fontsize',16);
subplot(4,1,2)
axis([-inf inf -2 18])
ylabel('zscore( DF/F )')
aux=(abs(phase_f-mean_p(max_cell_s)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[244,209,181]/255,'linewidth',2);
end
plot(zscore(dff((max_cell_s),:)),'Color','k','linewidth',1.5);
% set(gca,'XColor', 'none','YColor','none')
set(gca,'XColor', 'none')
yticks([0 8 16])
set(gca,'fontsize', 16);
subplot(4,1,3)
axis([-inf inf -2 18])
aux=(abs(phase_f-mean_p(med_cell_s)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[183,212,219]/255,'linewidth',2);
end
plot(zscore(dff((med_cell_s),:)),'Color','k','linewidth',1.5);
ylabel('zscore( DF/F )')
yticks([0 8 16])
set(gca,'XColor', 'none');
set(gca,'fontsize', 16);
clear aux vals timepoints
subplot(4,1,4)
axis([-inf inf -2 18])
aux=(abs(phase_f-mean_p(min_cell_s)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
hold on
for i=1:length(timepoints)
plot([timepoints(i),timepoints(i)],[-2,18],'-','color',[149,125,173]/255,'linewidth',2);
end
plot(zscore(dff((min_cell_s),:)),'Color','k','linewidth',1.5);
xticks([1,floor(T/2),T]);
xticklabels([0,floor(T/2)/8,T/8]);
ylabel('zscore( DF/F )');
yticks([0 8 16]);
% set(gca,'fontsize', 16);
% set(gca,'XColor', 'none','YColor','none')
% set(gca,'YColor','none')
set(gca,'fontsize',18)
xlabel('Time (s)');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures of examples of cells with MVL

col_dots=[0.9283,0.4730,0.3261];

cell_high_locking=18; %MVL=0.94
cell_high_locking_2=30; %MVL=0.91
cell_high_locking_3=474; %MVL=0.90
cell_med_locking=97; %MVL=0.62
cell_low_locking_not_locked=43; %MVL=0.1
cell_low_locking=217; %MVL=0.4
cell_low_locking_not_locked2=449; %MVL=0.15
cell_med_locking_2=86; %MVL=0.59

PR_cell_high_locking=PR_mat(find(PR_mat(:,1)==cell_high_locking),2); %0.6
PR_cell_high_locking_2=PR_mat(find(PR_mat(:,1)==cell_high_locking_2),2); %0.68
PR_cell_high_locking_3=PR_mat(find(PR_mat(:,1)==cell_high_locking_3),2); %0.56
PR_cell_med_locking=PR_mat(find(PR_mat(:,1)==cell_med_locking),2); %0.36
PR_cell_low_locking_not_locked=PR_mat(find(PR_mat(:,1)==cell_low_locking_not_locked),2); %0.64
PR_cell_low_locking=PR_mat(find(PR_mat(:,1)==cell_low_locking),2); %0.56
PR_cell_low_locking_not_locked2=PR_mat(find(PR_mat(:,1)==cell_low_locking_not_locked2),2); %0.4
PR_cell_med_locking_2=PR_mat(find(PR_mat(:,1)==cell_med_locking_2),2); %0.48


for w=1:size(table_u,1)
    spikes_per_wave_cell_high_locking(w)=sum(spikes(cell_high_locking,table_u(w,1):table_u(w,2)));
end
for w=1:size(table_u,1)
    spikes_per_wave_cell_high_locking_2(w)=sum(spikes(cell_high_locking_2,table_u(w,1):table_u(w,2)));
end
for w=1:size(table_u,1)
    spikes_per_wave_cell_high_locking_3(w)=sum(spikes(cell_high_locking_3,table_u(w,1):table_u(w,2)));
end
for w=1:size(table_u,1)
    spikes_per_wave_cell_med_locking(w)=sum(spikes(cell_med_locking,table_u(w,1):table_u(w,2)));
end
for w=1:size(table_u,1)
    spikes_per_wave_cell_low_locking_not_locked(w)=sum(spikes(cell_low_locking_not_locked,table_u(w,1):table_u(w,2)));
end
for w=1:size(table_u,1)
    spikes_per_wave_cell_low_locking(w)=sum(spikes(cell_low_locking,table_u(w,1):table_u(w,2)));
end
for w=1:size(table_u,1)
    spikes_per_wave_cell_low_locking_not_locked2(w)=sum(spikes(cell_low_locking_not_locked2,table_u(w,1):table_u(w,2)));
end
for w=1:size(table_u,1)
    spikes_per_wave_cell_med_locking_2(w)=sum(spikes(cell_med_locking_2,table_u(w,1):table_u(w,2)));
end

cumsum_spikes_per_wave_cell_high_locking=cumsum(sort(spikes_per_wave_cell_high_locking,'descend')./sum(spikes_per_wave_cell_high_locking));
cumsum_spikes_per_wave_cell_high_locking_2=cumsum(sort(spikes_per_wave_cell_high_locking_2,'descend')./sum(spikes_per_wave_cell_high_locking_2));
cumsum_spikes_per_wave_cell_high_locking_3=cumsum(sort(spikes_per_wave_cell_high_locking_3,'descend')./sum(spikes_per_wave_cell_high_locking_3));
cumsum_spikes_per_wave_cell_med_locking=cumsum(sort(spikes_per_wave_cell_med_locking,'descend')./sum(spikes_per_wave_cell_med_locking));
cumsum_spikes_per_wave_cell_low_locking_not_locked=cumsum(sort(spikes_per_wave_cell_low_locking_not_locked,'descend')./sum(spikes_per_wave_cell_low_locking_not_locked));
cumsum_spikes_per_wave_cell_low_locking=cumsum(sort(spikes_per_wave_cell_low_locking,'descend')./sum(spikes_per_wave_cell_low_locking));
cumsum_spikes_per_wave_cell_low_locking_not_locked2=cumsum(sort(spikes_per_wave_cell_low_locking_not_locked2,'descend')./sum(spikes_per_wave_cell_low_locking_not_locked2));
cumsum_spikes_per_wave_cell_med_locking_2=cumsum(sort(spikes_per_wave_cell_med_locking_2,'descend')./sum(spikes_per_wave_cell_med_locking_2));

cc=plasma(5);
figure
plot(cumsum_spikes_per_wave_cell_med_locking,'-','Linewidth',2.5,'Color',cc(1,:));
hold on
plot(cumsum_spikes_per_wave_cell_low_locking,'-','Linewidth',2.5,'Color',cc(2,:));
plot(cumsum_spikes_per_wave_cell_high_locking_3,'-','Linewidth',2.5,'Color',cc(3,:));
plot(cumsum_spikes_per_wave_cell_high_locking_2,'-','Linewidth',2.5,'Color',cc(4,:));
ylabel({'Normalized cumulative';'sum of spikes'});
xlabel('Wave #');
set(gca,'fontsize',18)
axis([0.5 26 0 1.1])
xticks([1 25])
box off
legend({'PI=0.36 - MVL=0.62';'PI=0.56 - MVL=0.4';'PI=0.56 - MVL=0.90';'PI=0.68 - MVL=0.91'})
legend boxoff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Examples of DFF of individual cells


% max_cell=412;
figure
aux=(abs(phase_f-mean_p(max_cell)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((max_cell),:)));
plot(zscore(dff((max_cell),:))- 1.1*min(zscore(dff((max_cell),:))),'Color',cc(3,:));
box off
yticks([])
axis off
hold on
scatter(timepoints,4*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.9
h=area(1:T,(max_val_dff- 1.7*min(zscore(dff((max_cell),:)))).*wave_epochs,'LineStyle','-');
h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
h(1).FaceColor = col_area;
h(1).FaceAlpha = 0.1;
hold on
clear aux vals rimepoints lox_local_minima max_val_dff



% min_cell=97;
figure
plot(zscore(dff((min_cell),:))- 1.1*min(zscore(dff((min_cell),:))),'Color',cc(1,:));
box off
axis([-inf inf -inf inf])
yticks([])
axis off
hold on
aux=(abs(phase_f-mean_p(min_cell)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((min_cell),:)));
scatter(timepoints,2*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.9
h=area(1:T,(max_val_dff- 1.9*min(zscore(dff((cell_high_locking_2),:)))).*wave_epochs,'LineStyle','-');
h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
h(1).FaceColor = col_area;
h(1).FaceAlpha = 0.1;
hold on
clear aux vals rimepoints lox_local_minima max_val_dff

% med_cell=271;
figure
plot(zscore(dff((med_cell),:))- 1.1*min(zscore(dff((med_cell),:))),'Color',cc(2,:));
box 
axis([-inf inf -inf inf])
yticks([])
axis off
hold on
aux=(abs(phase_f-mean_p(med_cell)));
[vals,timepoints]=findpeaks(-aux);
low_local_minima=find(aux(timepoints)<0.1);
timepoints=timepoints(low_local_minima);
max_val_dff=max(zscore(dff((med_cell),:)));
scatter(timepoints,3*ones(1,length(timepoints)),70,col_dots,'filled');
alpha 0.9
h=area(1:T,(max_val_dff- 1.9*min(zscore(dff((med_cell),:)))).*wave_epochs,'LineStyle','-');
h(1).EdgeColor = 'none';%[0.2 0.1 0.1];
h(1).FaceColor = col_area;
h(1).FaceAlpha = 0.1;
hold on
clear aux vals rimepoints lox_local_minima max_val_dff



%% Coupling to Population and Sensitivity - Using Tuning curve  % CONDITIONED ON WAVES

N_sh=200;
% spikes=spikes_d;
% phase=phase_r;  phase conditioned on waves
% spikes_d=[];
% spikes_d=spikes_r; spike matrix conditioned on waves

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
figure
plot(prob_k)

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
% figure
% scatter(inferred_ensemble,real_ensemble,'o','filled');
% axis([0 11 0 11])
% alpha 0.1
% xlabel('Inferred Ensemble - AUC')
% ylabel('Real Ensemble')
% set(gca,'fontsize',16)
% hold on

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


% H=hist3([inferred_ensemble',real_ensemble']);
% figure
% imagesc(H);
% colormap magma
% xticks([1 5 10])
% yticks([1 5 10])
% ylabel('');
% xlabel('Real ensemble');
% ylabel('Inferred ensemble')
% axis square
% set(gca,'fontsize',18)
% colorbar
% caxis([0 40])

% H=hist3([inferred_ensemble_sh';real_ensemble']);
% figure
% imagesc(H);
% colormap magma
% xticks([1 5 10])
% yticks([1 5 10])
% ylabel('');
% xlabel('Real ensemble');
% ylabel('Inferred ensemble - Shuffle')
% axis square
% set(gca,'fontsize',18)
% colorbar
% caxis([0 40])

figure
scatter(inferred_ensemble_sens_sh(:,3),real_ensemble,'o','filled');
axis([0 11 0 11])
alpha 0.1
xlabel('Inferred Ensemble - Sens')
ylabel('Real Ensemble')
set(gca,'fontsize',16)
hold on 
axis square
xticks([1 5 10])
yticks([1 5 10])

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Estimation of error in the prediction

correct=find(inferred_ensemble_sens==real_ensemble);
accuracy=length(correct)/N;
dist=abs(inferred_ensemble_sens-real_ensemble);

for n=1:N
   if dist(n)==9 dist(n)=1;
   elseif dist(n)==8 dist(n)=2;
   elseif dist(n)==7 dist(n)=3;
   elseif dist(n)==6 dist(n)=4;
   end       
end

MSE = sum(dist.*dist)/N;

for sh=1:200
    dist_sh(:,sh)=abs(inferred_ensemble_sens_sh(:,sh)-real_ensemble');
end

dist_sh_copy=dist_sh;

for sh=1:200    
    for n=1:N        
        if dist_sh_copy(n,sh)==9 dist_sh_copy(n,sh)=1;
        elseif dist_sh_copy(n,sh)==8 dist_sh_copy(n,sh)=2;
        elseif dist_sh_copy(n,sh)==7 dist_sh_copy(n,sh)=3;
        elseif dist_sh_copy(n,sh)==6 dist_sh_copy(n,sh)=4;
        end
    end
    MSE_sh(sh)=sum(dist_sh_copy(:,sh).*dist_sh_copy(:,sh))/N;
end

col=magma(7);
figure
bar([1,2],[MSE,mean(MSE_sh)],'FaceColor',col(4,:));
hold on
er= errorbar([1,2],[MSE,mean(MSE_sh)],[0,std(MSE_sh)]);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1.5;
ylabel('MSE');
xticks([1,2]);
xticklabels([{'Data'},{'Shuffle'}]);
set(gca,'fontsize',20)
box off
 
% figure
% errorbar(dist,mean(dist_sh,2),std(dist_sh,[],2),'*')

for sh=1:200   
    correct_sh=find(inferred_ensemble_sens_sh(:,sh)==real_ensemble');
    accuracy_sh(sh)=length(correct_sh)/N;
    clear correct_sh
end

figure
bar([1,2],[accuracy,mean(accuracy_sh)]);
hold on
er= errorbar([1,2],[accuracy,mean(accuracy_sh)],[0,std(accuracy_sh)]);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1.5;
ylabel('Accuracy');
xticks([1,2]);
xticklabels([{'Data'},{'Shuffle'}]);
set(gca,'fontsize',16)
axis([0 3 0 1])
box off

%% Population coupling - Using Pearson correlation - CONDITIONED ON WAVES

% zdff=zscore(dff')';
% 
% new_bin = 4;
% for i=1:length(phase_d)
%     sp_do(:,i)=sum(spikes_d(:,(i-1)*new_bin+1:i*new_bin),2);
% end

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
alpha=0.2
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


figure
imagesc(res_sh(sorting_w,:,1));
ylabel('Neurons #');
xlabel('Ensemble')
colormap plasma
xticks([1,5,10])
caxis([-0.2 0.5])
set(gca,'fontsize',18)
colorbar

thr=prctile(delta_coupling_sh',99);
alfa=find(delta_coupling>thr);
not_coupling=find(delta_coupling<=thr);
figure
scatter(thr,delta_coupling,'ko','filled')
alpha 0.5
hold on
axis([0 0.7 0 0.7])
hline=refline(1,0);
hline.LineStyle='--';
ylabel('Max coupling - Min coupling');
xlabel('Max coupling - Min coupling - 99th prctile');
set(gca,'fontsize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% single
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
figure
scatter(MI,delta_coupling,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255);
alpha 0.5
xlabel('MI (bits)');
ylabel('Max coupling - Min coupling');
set(gca,'fontsize',16)
hold on
scatter(MI(not_locked),delta_coupling(not_locked),40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor','r')
alpha 0.6
set(gca,'fontsize',18)
axis([-0.02 0.2 0 1])
box off

figure
scatter(MVL,delta_coupling,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255);
alpha 0.5
xlabel('MVL');
ylabel('Max coupling - Min coupling');
set(gca,'fontsize',16)
hold on
scatter(MVL(not_locked),delta_coupling(not_locked),40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor','r')
alpha 0.6
ylabel('Max coupling - Min coupling');
xlabel('MVL');
set(gca,'fontsize',18)
box off

figure
scatter(radius_f,delta_coupling,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255);
alpha 0.5
xlabel('Radius');
ylabel('Max coupling - Min coupling');
set(gca,'fontsize',16)
hold on
scatter(radius_f(not_locked),delta_coupling(not_locked),40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor','r')
alpha 0.6
set(gca,'fontsize',18)
box off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tuning curves for individual tuning curves

cell_high_locking=18; %MVL=0.94
cell_med_locking=97; %MVL=0.62
cell_low_locking_not_locked=43; %MVL=0.1
cell_low_locking=217; %MVL=0.4
cell_low_locking_not_locked2=449; %MVL=0.15
cell_med_locking_2=86; %MVL=0.59
cell_high_locking_2=30; %MVL=0.91
cell_high_locking_3=474; %MVL=0.90

figure
plot(1:10,res(cell_high_locking_3,:),'k','linewidth',3);
box off
hold on
errorbar(1:10,mean(res_sh(cell_high_locking_3,:,:),3),std(res_sh(cell_high_locking_3,:,:),[],3),'r','linewidth',2);
set(gca,'fontsize',16);
xlabel('Ensemble #');
ylabel('Coupling to ensemble activity');
axis([1 10 -0.3 0.5])
xticks([1 5 10])

figure
plot(1:10,res(cell_high_locking_2,:),'k','linewidth',3);
box off
hold on
errorbar(1:10,mean(res_sh(cell_high_locking_2,:,:),3),std(res_sh(cell_high_locking_2,:,:),[],3),'r','linewidth',2);
set(gca,'fontsize',16);
xlabel('Ensemble #');
ylabel('Coupling to ensemble activity');
axis([1 10 -0.3 0.5])
xticks([1 5 10])

figure
plot(1:10,res(cell_high_locking,:),'k','linewidth',3);
box off
hold on
errorbar(1:10,mean(res_sh(cell_high_locking,:,:),3),std(res_sh(cell_high_locking,:,:),[],3),'r','linewidth',2);
set(gca,'fontsize',16);
xlabel('Ensemble #');
ylabel('Coupling to ensemble activity');
axis([1 10 -0.3 0.5])
xticks([1 5 10])


figure
plot(1:10,res(cell_med_locking_2,:),'k','linewidth',3);
box off
hold on
errorbar(1:10,mean(res_sh(cell_med_locking_2,:,:),3),std(res_sh(cell_med_locking_2,:,:),[],3),'r','linewidth',2);
set(gca,'fontsize',16);
xlabel('Ensemble #');
ylabel('Coupling to ensemble activity');
axis([1 10 -0.3 0.5])
xticks([1 5 10])


figure
plot(1:10,res(cell_med_locking,:),'k','linewidth',3);
box off
hold on
errorbar(1:10,mean(res_sh(cell_med_locking,:,:),3),std(res_sh(cell_med_locking,:,:),[],3),'r','linewidth',2);
set(gca,'fontsize',16);
xlabel('Ensemble #');
ylabel('Coupling to ensemble activity');
axis([1 10 -0.3 0.5])

figure
plot(1:10,res(cell_low_locking,:),'k','linewidth',3);
box off
hold on
errorbar(1:10,mean(res_sh(cell_low_locking,:,:),3),std(res_sh(cell_low_locking,:,:),[],3),'r','linewidth',2);
set(gca,'fontsize',16);
xlabel('Ensemble #');
ylabel('Coupling to ensemble activity');
axis([1 10 -0.3 0.5])


figure
plot(1:10,res(cell_low_locking_not_locked,:),'k','linewidth',3);
box off
hold on
errorbar(1:10,mean(res_sh(cell_low_locking_not_locked,:,:),3),std(res_sh(cell_low_locking_not_locked,:,:),[],3),'r','linewidth',2);
set(gca,'fontsize',16);
xlabel('Ensemble #');
ylabel('Coupling to ensemble activity');
axis([1 10 -0.3 0.5])


figure
plot(1:10,res(cell_low_locking_not_locked2,:),'k','linewidth',3);
box off
hold on
errorbar(1:10,mean(res_sh(cell_low_locking_not_locked2,:,:),3),std(res_sh(cell_low_locking_not_locked2,:,:),[],3),'r','linewidth',2);
set(gca,'fontsize',16);
xlabel('Ensemble #');
ylabel('Coupling to ensemble activity');
axis([1 10 -0.3 0.5])
% 
% for i=1:length(not_locked)
%     
%     min_locking=not_locked(i);
%     figure
%     plot(1:10,res(min_locking,:),'k','linewidth',3);
%     box off
%     hold on
%     errorbar(1:10,mean(res_sh(min_locking,:,:),3),std(res_sh(min_locking,:,:),[],3),'r','linewidth',2);
%     set(gca,'fontsize',16);
%     xlabel('Ensemble #');
%     ylabel('Coupling to ensemble activity');
%     axis([1 10 -0.3 0.5])
%     
% end



% thr=prctile(delta_coupling_sh',99);
% alfa=find(delta_coupling>thr);
% figure
% scatter(thr,delta_coupling,'ko','filled')
% alpha 0.5
% hold on
% scatter(thr(not_locked),delta_coupling(not_locked),'ro','filled')
% axis([0 0.7 0 0.7])
% hline=refline(1,0);
% hline.LineStyle='--';
% ylabel('Max coupling - Min coupling');
% xlabel('Max coupling - Min coupling - 99th prctile');
% set(gca,'fontsize',16)
% 
% [~,large_delta_coupling]=max(delta_coupling);
% [~,min_delta_coupling]=min(delta_coupling);
% 
% figure
% plot(1:10,res(large_delta_coupling,:),'k','linewidth',3);
% box off
% hold on
% errorbar(1:10,mean(res_sh(large_locking,:,:),3),std(res_sh(large_locking,:,:),[],3),'r','linewidth',2);
% set(gca,'fontsize',16);
% xlabel('Ensemble #');
% ylabel('Coupling to ensemble activity');
% axis([1 10 -0.3 0.5])

figure
plot(1:10,res(min_delta_coupling,:),'k','linewidth',3);
box off
hold on
errorbar(1:10,mean(res_sh(large_locking,:,:),3),std(res_sh(large_locking,:,:),[],3),'r','linewidth',2);
set(gca,'fontsize',16);
xlabel('Ensemble #');
ylabel('Coupling to ensemble activity');
axis([1 10 -0.3 0.5])

% tot_cells=1:N;
% coupled=setdiff(tot_cells,not_coupling);
% locked=setdiff(tot_cells,not_locked);
% 
% prop_locked_coupled = length(intersect(coupled,locked))./N;
% prop_Nlocked_coupled = length(intersect(coupled,not_locked))./N;
% prop_locked_Ncoupled = length(intersect(not_coupling,locked))./N;
% prop_Nlocked_Ncoupled = length(intersect(not_coupling,not_locked))./N;
% 
% figure
% pie([prop_locked_coupled prop_locked_Ncoupled prop_Nlocked_coupled prop_Nlocked_Ncoupled])
% legend


% % % % % %% Cross correlation
% % % % % [~,sorting,~]=get_sorting(spikes_d);
% % % % % 
% % % % % maxlag=120;
% % % % % downsampling_factor=4;
% % % % % FRp = spikes_downsample(spikes_d,N,downsampling_factor);
% % % % % for i=1:N
% % % % %     FR(i,:)=full(fire_rate(FRp(i,:),3,'g')); %Smoothing 
% % % % % end
% % % % % FR=FRp;
% % % % % count=0;
% % % % % 
% % % % % %I choose 3 cellsbased on their locking 
% % % % % 
% % % % % [small_deltaMI,cell_smallMI]=min(MVL(not_locked));
% % % % % PR_cell_smallMI=PR_cell(non_inf(cell_smallMI));
% % % % % [medium_deltaMI,cell_mediumMI]=min(MVL(locked));
% % % % % PR_cell_mediumMI=PR_cell(infor(cell_mediumMI));
% % % % % [medium_largeMI,cell_largeMI]=max(MVL(locked));
% % % % % PR_cell_largeMI=PR_cell(infor(cell_largeMI));
% % % % % 
% % % % % cells_t=1:N;
% % % % % cells_t(not_locked)=[];
% % % % % % [cells_without_noninf,cell_min]=min(MI_TB(cells_t));
% % % % % % [cells_without_inf,cell_max]=max(MI_TB(cells_t));
% % % % % % cell_min=342;
% % % % % %I calculate cross correlations
% % % % % for i=1:N
% % % % %     for j=i+1:N
% % % % %         [val,time]=(xcorr(FR(i,:),FR(j,:),maxlag)); %Check whether I need to zscore
% % % % % 
% % % % %         [v,in]=max(zscore(val));     
% % % % %         
% % % % %         if time(in)>=0
% % % % %             Corr_mat(i,j)=v;
% % % % %             Corr_mat(j,i)=-v;
% % % % %             lag_mat(i,j)=time(in);
% % % % %             lag_mat(j,i)=-time(in);            
% % % % %         else
% % % % %             Corr_mat(i,j)=-v;
% % % % %             Corr_mat(j,i)=v;
% % % % %             lag_mat(i,j)=time(in);
% % % % %             lag_mat(j,i)=-time(in);
% % % % %         end        
% % % % %         clear val time
% % % % %     end
% % % % % end
% % % % % 
% % % % % [~,sorting_smallMI]=sort(Corr_mat(not_locked(cell_smallMI),:),'descend');
% % % % % [~,sorting_mediumMI]=sort(Corr_mat(locked(cell_mediumMI),:),'descend');
% % % % % [~,sorting_largeMI]=sort(Corr_mat(locked(cell_largeMI),:),'descend');
% % % % % 
% % % % % cell_min=not_locked(cell_smallMI);
% % % % % cell_med=locked(cell_mediumMI);
% % % % % cell_max=locked(cell_largeMI);
% % % % % 
% % % % % count=0;
% % % % % for i=[cell_min,cell_med,cell_max]%[non_inf(cell_smallMI),infor(cell_mediumMI),infor(cell_largeMI)]
% % % % %     count=count+1;
% % % % %     for j=1:N
% % % % %     xcorr_mat(j,:,count)=(zscore(xcorr(FR(i,:),FR(j,:),maxlag)));      
% % % % %     end
% % % % % end
% % % % % 
% % % % % figure
% % % % % imagesc(xcorr_mat(sorting,:,1));
% % % % % colormap magma
% % % % % caxis([-2 5])
% % % % % xticks([1 120 240])
% % % % % xticklabels([-60 0 60]);
% % % % % yticks([1 250 500])
% % % % % xlabel('Time (s)')
% % % % % ylabel('Neurons #')
% % % % % set(gca,'fontsize',26)
% % % % % axis square
% % % % % colorbar
% % % % % 
% % % % % figure
% % % % % imagesc(xcorr_mat(sorting,:,2));
% % % % % colormap magma
% % % % % caxis([-2 5])
% % % % % xticks([1 120 240])
% % % % % xticklabels([-60 0 60]);
% % % % % yticks([1 250 500])
% % % % % xlabel('Time (s)')
% % % % % set(gca,'fontsize',26)
% % % % % axis square
% % % % % colorbar
% % % % % 
% % % % % figure
% % % % % imagesc(xcorr_mat(sorting,:,3));
% % % % % colormap magma
% % % % % caxis([-2 5])
% % % % % xticks([1 120 240])
% % % % % xticklabels([-60 0 60]);
% % % % % yticks([1 250 500])
% % % % % xlabel('Time (s)')
% % % % % set(gca,'fontsize',26)
% % % % % axis square
% % % % % colorbar
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % count=0;
% % % % % for i=14%1:20%[non_inf(cell_smallMI),infor(cell_mediumMI),infor(cell_largeMI)]
% % % % %     count=count+1;
% % % % %     cell_nl=locked(i);
% % % % %     for j=1:N
% % % % %         xcorr_mat(j,:,count)=(zscore(xcorr(FR(cell_nl,:),FR(j,:),maxlag)));
% % % % %     end
% % % % %     
% % % % %     figure
% % % % %     imagesc(xcorr_mat(sorting,:,i));
% % % % % % title(MVL(locked(i)))
% % % % %     colormap magma
% % % % %     caxis([-2 5])
% % % % %     xticks([1 120 240])
% % % % %     xticklabels([-60 0 60]);
% % % % %     yticks([1 250 500])
% % % % %     xlabel('Time (s)')
% % % % %     ylabel('Neurons #')
% % % % %     set(gca,'fontsize',26)
% % % % %     axis square
% % % % %     colorbar
% % % % % end
% % % % % 
% % % % % 
% % % % % %% MI and PD
% % % % % 
% % % % % figure
% % % % % scatter(MI_TB(non_inf),PR_cell(non_inf),70,'filled');
% % % % % xlabel('MI of non informative cells');
% % % % % ylabel('PI of non informative cells');
% % % % % set(gca,'fontsize',16)
% % % % % 
% % % % % figure
% % % % % plot(dff(non_inf(cell_smallMI),:),'k','linewidth',2)
% % % % % box off
% % % % % figure
% % % % % plot(dff(infor(cell_mediumMI),:),'k','linewidth',2)
% % % % % box off
% % % % % xticks([])
% % % % % yticks([])
% % % % % figure
% % % % % plot(dff(infor(cell_largeMI),:),'k','linewidth',2)
% % % % % box off
% % % % % xticks([])
% % % % % yticks([])
% % % % % 
% % % % % 
% % % % % %% Distribution of PR per ensemble
% % % % % 
% % % % % figure
% % % % % scatter(PR(sorting(Ens(:,1))),(Ens(:,2)),'o','filled','markeredgecolor',[0.2902    0.2863    0.2863],'markerfacecolor',[ 0.6510    0.6510    0.6510])
% % % % % xlabel('PI')
% % % % % ylabel('Ensemble #')
% % % % % alpha 0.5
% % % % % axis([0 1 0 11])
% % % % % set(gca,'fontsize',16)
% % % % % yticks([1:10])
% % % % % 
% % % % % %% Locking restricted to waves 
% % % % % 
% % % % % 
% % % % % for i=1:N    
% % % % %     spikes_per_wave=[];
% % % % %     phase_per_wave=[];
% % % % %     
% % % % %     for w=1:size(table_u,1)
% % % % %         spikes_per_wave=[spikes_per_wave,(spikes_d(i,table_u(w,1):table_u(w,2)))];
% % % % %         phase_per_wave=[phase_per_wave,phase(table_u(w,1):table_u(w,2))'];
% % % % %     end
% % % % %     
% % % % %     p=phase_per_wave(find(spikes_per_wave))';
% % % % %     var_p_wave(i)=circ_var(p);
% % % % %     mean_p_wave(i)=circ_mean(p);
% % % % %     std_p_wave(i)=circ_std(p);
% % % % %     MVL_wave(i) = circ_r(p);
% % % % %     %     H=histogram(p,-pi:2*pi/40:pi,'Normalization','Probability');
% % % % %     %     prob_phase_firing(i,:)=H.Values;
% % % % %     clear p H
% % % % %     
% % % % % end
% % % % % 
% % % % % for i=1:N
% % % % %     disp(i)
% % % % %     for sh=1:250
% % % % %         
% % % % %         spikes_per_wave=[];
% % % % %         phase_per_wave=[];
% % % % %         
% % % % %         for w=1:size(table_u,1)
% % % % %             spikes_per_wave=[spikes_per_wave,(spikes_d(i,table_u(w,1):table_u(w,2)))];
% % % % %             phase_per_wave=[phase_per_wave,phase(table_u(w,1):table_u(w,2))'];
% % % % %         end
% % % % %         
% % % % %     
% % % % %         p_sh=phase_per_wave(floor((length(spikes_per_wave))*rand(1,length(find(spikes_per_wave)))+1))';
% % % % %         var_p__sh_wave(i,sh)=circ_var(p_sh);
% % % % %         std_p__sh_wave(i,sh)=circ_std(p_sh);
% % % % %         mean_p_sh_wave(i,sh)=circ_mean(p_sh);
% % % % %         MVL_sh_wave(i,sh) = circ_r(p_sh);
% % % % %         if sh==1
% % % % %             H=histogram(p_sh,-pi:2*pi/40:pi,'Normalization','Probability');
% % % % %             prob_phase_firing_sh(i,:)=H.Values;
% % % % %             clear H
% % % % %         end
% % % % % %         H=histogram(p_sh,-pi:2*pi/40:pi,'Normalization','Probability');
% % % % % %         prob_phase_firing_sh_tot(i,:,sh)=H.Values;
% % % % % %         clear H
% % % % %         
% % % % %         clear p_sh
% % % % %     end
% % % % % end
% % % % %         
% % % % % locking_wave=1-var_p_wave;
% % % % % locking_sh_wave=1-var_p__sh_wave;
% % % % % locking_sh_mean_wave=mean(MVL_sh_wave,2);
% % % % % locking_sh_99_wave=prctile(MVL_sh_wave,99,2);
% % % % % locking_sh_1_wave=prctile(MVL_sh_wave,1,2);
% % % % % 
% % % % % cells_mi_calc=find(locking_wave>locking_sh_99_wave');
% % % % % tot_cells=1:N;
% % % % % not_locked_wave=find(locking_wave<=locking_sh_99_wave');
% % % % % locked_wave=setdiff(tot_cells,not_locked_wave);
% % % % % 
% % % % % figure
% % % % % labels = {'Locked | Wave','Not locked | Wave'};
% % % % % pie([length(locked_wave)/N,length(not_locked_wave)/N],labels)
% % % % % colormap jet(4)
% % % % % set(gca,'fontsize',16)
% % % % % 
% % % % % [a,b]=sort(locking_wave,'ascend');
% % % % % figure
% % % % % hold on
% % % % % plot(locking_sh_99_wave(b),'.','MarkerSize',15,'Color','r');
% % % % % plot(locking_wave(b),'.','MarkerSize',15,'Color','k');
% % % % % alpha 0.6
% % % % % ylabel('Locking to phase');
% % % % % xlabel('Neuron #')
% % % % % set(gca,'fontsize',18)
% % % % % box off
% % % % % legend('Shuffle - 99th percentile', 'Data')
% % % % % %                 legend('Control - Mean', 'Control - 1th percentile', 'Control - 99th percentile', 'Data')
% % % % % legend boxoff
% % % % %                       
% % % % % 
% % % % % figure
% % % % % scatter(radius,MVL_wave,55,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255)
% % % % % alpha 0.7
% % % % % ylabel('Locking to phase');
% % % % % xlabel('Radius PC1-PC2');
% % % % % set(gca,'fontsize',18)
% % % % % axis([-0.02 0.2 0 1])
% % % % % box off
% % % % %                             
% % % % %                 
% % % % % table_cells(:,1)=radius;
% % % % % table_cells(:,2)=MVL_wave;
% % % % % table_cells(:,3)=sum(spikes_d,2);
% % % % % spikes_count_dis=discretize(table_cells(:,3),2);
% % % % % figure
% % % % % colormap viridis(2)
% % % % % scatter(radius,MVL_wave,55,spikes_count_dis,'filled')
% % % % % alpha 0.7
% % % % % ylabel('Locking to phase');
% % % % % xlabel('Radius PC1-PC2');
% % % % % set(gca,'fontsize',18)
% % % % % axis([-0.02 0.2 0 1])
% % % % % box off
% % % % % colorbar
% % % % % 
% % % % %                              
% % % % % 
% % % % % 


%% Coupling to Population and Sensitivity - Using Tuning curve  % *NOT* CONDITIONED ON WAVES

N_sh=200;
% spikes=spikes_d;
% phase=phase_r;  phase conditioned on waves
% spikes_d=[];
% spikes_d=spikes_r; spike matrix conditioned on waves

%Downsample and binarize
new_bin = 4;
num_bins=floor(size(spikes,2)/new_bin);
for i=1:num_bins %Length of phase conditioned on waves and previously downsampled to 4 bins (0.5 bin size) when computing MI
    sp_do(:,i)=sum(spikes(:,(i-1)*new_bin+1:i*new_bin),2);
end

for n=1:N
    sp_do(n,find(sp_do(n,:)>0))=1;
end

%Calculate coactivity
prob_k=rcoact_dist(sp_do);
figure
plot(prob_k)

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
% figure
% scatter(inferred_ensemble,real_ensemble,'o','filled');
% axis([0 11 0 11])
% alpha 0.1
% xlabel('Inferred Ensemble - AUC')
% ylabel('Real Ensemble')
% set(gca,'fontsize',16)
% hold on

figure
scatter(inferred_ensemble_sens,real_ensemble,'o','filled');
axis([0 11 0 11])
alpha 0.1
xlabel('Inferred Ensemble - Sens')
ylabel('Real Ensemble')
set(gca,'fontsize',16)
hold on 

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

% H=hist3([inferred_ensemble',real_ensemble']);
% figure
% imagesc(H);
% colormap magma
% xticks([1 5 10])
% yticks([1 5 10])
% ylabel('');
% xlabel('Real ensemble');
% ylabel('Inferred ensemble')
% axis square
% set(gca,'fontsize',18)
% colorbar
% caxis([0 40])

% H=hist3([inferred_ensemble_sh';real_ensemble']);
% figure
% imagesc(H);
% colormap magma
% xticks([1 5 10])
% yticks([1 5 10])
% ylabel('');
% xlabel('Real ensemble');
% ylabel('Inferred ensemble - Shuffle')
% axis square
% set(gca,'fontsize',18)
% colorbar
% caxis([0 40])

H=hist3([inferred_ensemble_sens_sh(:,1),real_ensemble']);
figure
imagesc(H);
colormap magma
xticks([1 5 10])
yticks([1 5 10])
ylabel('');
xlabel('Real ensemble');
ylabel('Inferred ensemble - Shuffle')
axis square
set(gca,'fontsize',18)
colorbar
caxis([0 40])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Estimation of error in the prediction

correct=find(inferred_ensemble_sens==real_ensemble);
accuracy=length(correct)/N;
dist=abs(inferred_ensemble_sens-real_ensemble);

for n=1:N
   if dist(n)==9 dist(n)=1;
   elseif dist(n)==8 dist(n)=2;
   elseif dist(n)==7 dist(n)=3;
   elseif dist(n)==6 dist(n)=4;
   end       
end

MSE = sum(dist.*dist)/N;

for sh=1:200
    dist_sh(:,sh)=abs(inferred_ensemble_sens_sh(:,sh)-real_ensemble');
end

dist_sh_copy=dist_sh;

for sh=1:200    
    for n=1:N        
        if dist_sh_copy(n,sh)==9 dist_sh_copy(n,sh)=1;
        elseif dist_sh_copy(n,sh)==8 dist_sh_copy(n,sh)=2;
        elseif dist_sh_copy(n,sh)==7 dist_sh_copy(n,sh)=3;
        elseif dist_sh_copy(n,sh)==6 dist_sh_copy(n,sh)=4;
        end
    end
    MSE_sh(sh)=sum(dist_sh_copy(:,sh).*dist_sh_copy(:,sh))/N;
end

col=magma(7);
figure
bar([1,2],[MSE,mean(MSE_sh)],'FaceColor',col(4,:));
hold on
er= errorbar([1,2],[MSE,mean(MSE_sh)],[0,std(MSE_sh)]);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1.5;
ylabel('MSE');
xticks([1,2]);
xticklabels([{'Data'},{'Shuffle'}]);
set(gca,'fontsize',20)
box off
 
% figure
% errorbar(dist,mean(dist_sh,2),std(dist_sh,[],2),'*')

for sh=1:200   
    correct_sh=find(inferred_ensemble_sens_sh(:,sh)==real_ensemble');
    accuracy_sh(sh)=length(correct_sh)/N;
    clear correct_sh
end

figure
bar([1,2],[accuracy,mean(accuracy_sh)]);
hold on
er= errorbar([1,2],[accuracy,mean(accuracy_sh)],[0,std(accuracy_sh)]);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1.5;
ylabel('Accuracy');
xticks([1,2]);
xticklabels([{'Data'},{'Shuffle'}]);
set(gca,'fontsize',16)
axis([0 3 0 1])
box off

