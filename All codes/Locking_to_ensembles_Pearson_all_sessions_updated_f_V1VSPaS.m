clear all

session_pas=load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_ensemble_Pearson_PaS\locking_ensembles_Pearson_pas.mat');
session_v1=load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_ensemble_Pearson_V1\locking_ensembles_Pearson_V1.mat');


n_ses_v1=19;
n_ses_pas=29;

max_corr_v1=[];
corr_v1=[];
for i=1:n_ses_v1
    aux=session_v1.session{1,i}.res;
    max_corr_v1=[max_corr_v1;max(aux,[],2)];
    corr_v1=[corr_v1;aux(:)];
    cells_v1(i)=size(session_v1.session{1,i}.res,1);
    clear aux
end

ses_pas=1:n_ses_pas;
ses_pas([19,21,23,25])=[];

max_corr_pas=[];
corr_pas=[];
for i=ses_pas
    disp(i);
    aux=session_pas.session{1,i}.res;
    max_corr_pas=[max_corr_pas;max(aux,[],2)];
    corr_pas=[corr_pas;aux(:)];
    cells_pas(i)=size(session_pas.session{1,i}.res,1);

    clear aux
end


h_v1=histcounts(max_corr_v1,-1:0.05:1,'Normalization','cdf');
h_pas=histcounts(max_corr_pas,-1:0.05:1,'Normalization','cdf');

[h,p,ks2stat] = kstest2(max_corr_v1,max_corr_pas,'tail','smaller');
[h,p,ks2stat] = kstest2(max_corr_v1,max_corr_pas);


figure
plot(-1+0.05:0.05:1,h_v1,'linewidth',2)
hold on
plot(-1+0.05:0.05:1,h_pas,'linewidth',2)
box off
axis([-0.2 1 0 1.05]);
ylabel('Cumulative probability');
xlabel({'Maximum correlation between';'single-cell and ensemble activity'});
set(gca,'fontsize',16);
legend('V1','PaS');
legend boxoff

% 
% 
% h_v1=histcounts(corr_v1,-1:0.05:1,'Normalization','cdf');
% h_pas=histcounts(corr_pas,-1:0.05:1,'Normalization','cdf');
% 
% [h,p,ks2stat] = kstest2(max_corr_v1,max_corr_pas,'tail','smaller');
% 
% figure
% plot(-1+0.05:0.05:1,h_v1,'linewidth',2)
% hold on
% plot(-1+0.05:0.05:1,h_pas,'linewidth',2)
% box off
% axis([-0.2 1 0 1.05]);
% ylabel('Cumulative probability');
% xlabel({'Maximum correlation between';'single-cell and ensemble activity'});
% set(gca,'fontsize',16);
% legend('V1','PaS');
% legend boxoff




