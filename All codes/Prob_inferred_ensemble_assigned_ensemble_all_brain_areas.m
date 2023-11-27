% 0.Load data

session_v1=load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_ensemble_Pearson_V1\locking_ensembles_Pearson_V1.mat');
session_pas=load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_ensemble_Pearson_PaS\locking_ensembles_Pearson_pas.mat');
load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_ensemble_Pearson_MEC\locking_ensembles_Pearson_mec.mat');

session_v1=session_v1.session;
session_pas=session_pas.session;

ses_pas=1:29;
ses_pas([19,21,23,25])=[];

% 1. Calculate probability of distance between inferred and real ensemble

%MEC
for i=1:size(session,2)
    prob_dist_mec(i,:)=histcounts(session{1,i}.dist,[0:6],'Normalization','Probability');
end

%V1
for i=1:size(session_v1,2)
    prob_dist_v1(i,:)=histcounts(session_v1{1,i}.dist,[0:6],'Normalization','Probability');
end

%PaS
count_pas=0;
for i=ses_pas
    count_pas=count_pas+1;
    prob_dist_pas(count_pas,:)=histcounts(session_pas{1,i}.dist,[0:6],'Normalization','Probability');
end

% 2. Repeats step 1 but for shuffled data

%MEC
for i=1:size(session,2)
    prob_dist_sh_mec(i,:)=histcounts(session{1,i}.dist_sh(:),[0:6],'Normalization','Probability');
end

%V1
for i=1:size(session_v1,2)
    prob_dist_sh_v1(i,:)=histcounts(session_v1{1,i}.dist_sh(:),[0:6],'Normalization','Probability');
end

%PaS
count_pas=0;
for i=ses_pas
    count_pas=count_pas+1;
    prob_dist_sh_pas(count_pas,:)=histcounts(session_pas{1,i}.dist_sh(:),[0:6],'Normalization','Probability');
end

% 3. Makes box plot for 2 brain areas


mat_PaS_V1=nan(25,3);
mat_PaS_V1(1:25,1)=prob_dist_pas(:,1);
mat_PaS_V1(1:25,2)=prob_dist_sh_pas(:,1);
mat_PaS_V1(1:19,3)=prob_dist_v1(:,1);
mat_PaS_V1(1:19,4)=prob_dist_sh_v1(:,1);

figure
boxplot(mat_PaS_V1,'Notch','off','Labels',{'PaS','PaS shuffle','VIS', 'VIS shuffle'})
ylabel({'Probability' ;'inferred ensemble = assigned ensemble'});
box off
set(gca,'fontsize',18)
yticks([0.1 0.2 0.3 0.4 0.5 0.6])

[p_V1,h_V1,stat_V1]=ranksum(mat_PaS_V1(1:19,3),mat_PaS_V1(1:19,4));
[p_PaS,h_PaS,stat_PaS]=ranksum(mat_PaS_V1(1:25,1),mat_PaS_V1(1:25,2));
[p_V1PaS,h_V1PaS,stat_V1PaS]=ranksum(mat_PaS_V1(1:19,1),mat_PaS_V1(1:25,3),'tail','left');

% 4. Makes box plot for 3 brain areas

figure
mat=nan(25,3);
mat(1:15,1)=prob_dist_mec(:,1);
mat(1:19,2)=prob_dist_v1(:,1);
mat(1:25,3)=prob_dist_pas(:,1);
figure
boxplot(mat,'Notch','off','Labels',{'MEC','V1','PaS'})
ylabel('Probability inferred = real');
box off
set(gca,'fontsize',18)

figure
bar([1,2,3],[mean_prob_mec(1),mean_prob_v1(1),mean_prob_pas(1)]);
hold on
errorbar([1,2,3],[mean_prob(1),mean_prob_v1(1),mean_prob_pas(1)],[sem_prob(1),sem_prob_v1(1),sem_prob_pas(1)],'k.','linewidth',2.5);
ylabel('Probability inferred = real');
box off
set(gca,'fontsize',18)
xticks=[1,2,3];
xticklabels({'MEC','V1','PaS'})


% 5. Computes the mean of the probabilities and plots them

mean_prob_mec=mean(prob_dist_mec);
sem_prob_mec=std(prob_dist_mec)./sqrt(15);

mean_prob_v1=mean(prob_dist_v1);
sem_prob_v1=std(prob_dist_v1)./sqrt(19);

mean_prob_pas=mean(prob_dist_pas);
sem_prob_pas=std(prob_dist_pas)./sqrt(25);

% mean_prob_sh_v1=mean(prob_dist_sh);
% sem_prob_sh_v1=std(prob_dist_sh)./sqrt(15);
% mean_prob_sh_pas=mean(prob_dist_sh_pas);
% sem_prob_sh_pas=std(prob_dist_sh_pas)./sqrt(15);


figure
errorbar([0:5],mean_prob_mec,sem_prob_mec,'linewidth',2.5);
hold on
errorbar([0:5],mean_prob_pas,sem_prob_pas,'linewidth',2.5);
errorbar([0:5],mean_prob_v1,sem_prob_v1,'linewidth',2.5);
xticks([0,1,2,3,4,5])
ylabel('Probability');
xlabel('Distance between inferred and real ensemble');
set(gca,'fontsize',18);
box off
legend('MEC','PaS','V1')
for i=1:6
    [p_mec_pas(i),h_mec_pas(i),stats_mec_pas(i)]=ranksum(prob_dist(:,i),prob_dist_pas(:,i));
end

for i=1:6
    [p_mec_v1(i),h_mec_v1(i),stats_mec_v1(i)]=ranksum(prob_dist(:,i),prob_dist_v1(:,i));
end

for i=1:6
    [p_pas_v1(i),h_pas_v1(i),stats_pas_v1(i)]=ranksum(prob_dist_pas(:,i),prob_dist_v1(:,i));
end


figure
plot([0:5],cumsum(mean_prob_mec),'linewidth',2.5);
hold on
plot([0:5],cumsum(mean_prob_pas),'linewidth',2.5);
plot([0:5],cumsum(mean_prob_v1),'linewidth',2.5);
xticks([0,1,2,3,4,5])
ylabel('Cumulative Probability');
xlabel('Distance between inferred and real ensemble');
set(gca,'fontsize',18);
box off



