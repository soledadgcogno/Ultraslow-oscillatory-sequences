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


%% V1 sessions


V1_ses=find(big_table(:,10)<0); %Sessions inthe in MEC
adults=find(big_table(:,3)>15); %Sessions inthe of adults
sessions=intersect(V1_ses,adults);
count=0;

ws=big_table(sessions,6);
index_nan=isnan(ws);
sessions(index_nan)=[];
session_number=big_table(sessions,4);
days=big_table(sessions,3);
repeated_days=find(diff(days)==0);
clear index_to_delete

index_to_delete=[];
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
V1_sessions=sessions;

%% Pearson locking to ensembles

n_sh_pearson=500;

clus=10;
count=0;
prctile_th=99;

for w=1:length(V1_sessions)
     row_w=V1_sessions(w);
    disp(w)
    
    count=count+1;
    
    if num2str(big_table(V1_sessions(w),2))=='1'
        mouse=92227;
    elseif num2str(big_table(V1_sessions(w),2))=='2'
        mouse=92229;
    elseif num2str(big_table(V1_sessions(w),2))=='3'
        mouse=60961;
    end
    
    count=count+1;
    day=big_table(V1_sessions(w),3);
    day_w=big_table(V1_sessions(w),3);
    s=big_table(V1_sessions(w),4);
    munit=big_table(V1_sessions(w),5);
    load([rec_data_path,strcat('recording_dates_',num2str(mouse),'.mat')]);
    if isfield(dates,'actual_day')==1
        day=find(dates.actual_day==day);
    end
    
    clus=10;
    disc_phase=10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
    load([rec_data_path,strcat('recording_dates_',num2str(mouse),'.mat')]);
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',num2str(mouse),'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

    load(file_name_spk,'-mat'); %Spike times
    spikes_d=full(spikes_d_s);
    [N,T]=size(spikes_d);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Unlike MEC there is no conditioning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% on waves
    
    spikes_r=spikes_d;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distribution of cells in ensembles
    
    clear Ens
    [~,sorting_w,~]=get_sorting(spikes_d); %Unlike MEC, we there is no smoothing kernel
%     [~,sorting_w,~]=get_sorting_smoothed(spikes_d,dt);
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Downsample spike matrix
    new_bin = 4;
    num_bins=floor(size(spikes_d,2)/new_bin);
    for i=1:num_bins %Length of phase conditioned on waves and previously downsampled to 4 bins (0.5 seconds bin size)
        sp_do(:,i)=sum(spikes_d(:,(i-1)*new_bin+1:i*new_bin),2);
    end
    
    for n=1:N
        sp_do(n,find(sp_do(n,:)>0))=1;
    end
    
    for n=1:N
%         disp(n)
        for ense=1:10
            cells_ens_aux=find(Ens(:,2)==ense);
            cells_ens=Ens(cells_ens_aux,1); %Cells that belong to ensemble ense
            cells=setdiff(cells_ens,n); %Removes from cells_ens cell *n* for which we are computing the tuning curve
            
            pop_sum=sum(sp_do(cells,:));
            [rho,pval]=corr(sp_do(n,:)',pop_sum');
            res(n,ense)=rho;
            
            for sh=1:n_sh_pearson
                temp=randperm(length(pop_sum));
                dff_sh=sp_do(n,temp); %shuffle spike train of cell "n"
                res_sh(n,ense,sh)=corr(dff_sh',pop_sum');
                clear temp
            end
        end
        
        delta_coupling(n)=max(res(n,:)) - min(res(n,:));
        aux=find(Ens(:,1)==n);
        real_ensemble(n)=Ens(aux,2);
        [~,inf_ensemble(n)]=max(res(n,:));
        
        dist(n)=abs(inf_ensemble(n)-real_ensemble(n));
        
        if dist(n)==9 dist(n)=1;
        elseif dist(n)==8 dist(n)=2;
        elseif dist(n)==7 dist(n)=3;
        elseif dist(n)==6 dist(n)=4;
        end
        
        for sh=1:n_sh_pearson
            delta_coupling_sh(n,sh)=max(res_sh(n,:,sh)) - min(res_sh(n,:,sh));
            [~,inf_ensemble_sh(n,sh)]=max(res_sh(n,:,sh));
            
            dist_sh(n,sh)=abs(inf_ensemble_sh(n,sh)-real_ensemble(n));
            if  dist_sh(n,sh)==9  dist_sh(n,sh)=1;
            elseif  dist_sh(n,sh)==8  dist_sh(n,sh)=2;
            elseif  dist_sh(n,sh)==7  dist_sh(n,sh)=3;
            elseif  dist_sh(n,sh)==6  dist_sh(n,sh)=4;
            end
        end   
        
        clear  pop_sum cells_ens_aux cells_ens cells 
          
    end
    
    %Store everything 
    
    session{1,w}.dist_sh=dist_sh;
    session{1,w}.dist=dist;
    session{1,w}.inf_ensemble_sh=inf_ensemble_sh;
    session{1,w}.inf_ensemble=inf_ensemble;
    session{1,w}.real_ensemble=real_ensemble;
    session{1,w}.delta_coupling=delta_coupling;
    session{1,w}.delta_coupling_sh=delta_coupling_sh;
    session{1,w}.res=res;
    session{1,w}.res_sh=res_sh;

    clear cells cells_di cells_d coeff dff exclude1 exclude2 latent_pca MI_TB phase phase_d phase_di score signal_dff snr SNR spikes spikes_d spikes_d_s p MVL
    clear not_locked locked locking locking_sh locking_sh_mean locking_sh_99 locking_sh_1 MVL_sh p_sh sp_do phase phase_d mean_p mean_p_sh prob_phase_firing
    clear std_p std_p__sh var_p var_p__sh calc_di coefft FRp MI MI_b MI_withbias phase_down phase_f phase_r radius_f scoret spikes_r spk_do ...
        table_u Anat d_within d_across delta_tissue delta_tissue_vec ens Ens ense_n  sorting_w r_i r_j PR new_mat new_mat2 spikes_sorted ...
        PR_locked PR_not_locked index_sorting_locked sorting_w_not_locked sorting_w_locked delta_phase_mean index_sorting_not_locked...
        Sens_sh AUC_sh TC2_sh TC3_sh Arg_sens2_sh Sens AUC TC2 TC3 Arg_sens2 p_1 p_0 mean_fr calc cell_wave_part prob_k sp_do ind ...
        cells cells_aux cells_ens_aux Ens sp_do res_sh res dist dist_sh  delta_coupling_sh inf_ensemble_sh inf_ensemble real_ensemble ...
        delta_coupling real_ensemble N n   
    
end

save('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_ensemble_Pearson_V1\locking_ensembles_Pearson_V1.mat','session')

  
  
%% Figures and test for significance

sessions=length(session);
% Test for significance
delta_coupling_conc=[];
count=0;

for w=1:sessions
    
    thr=prctile(session{1,w}.delta_coupling_sh',99);
    delta_coupling=session{1,w}.delta_coupling;
    
    % thr_locked=thr(locked);
    % delta_coupling_locked=delta_coupling(locked);
    
    prop_locked_to_ensembles(w)=length(find(delta_coupling>thr))/length(delta_coupling);
    
    
    dist=abs(session{1,w}.dist);
    dist_sh=abs(session{1,w}.dist_sh);
    
    prob_dist(w,:)=histcounts(dist,[0:6],'Normalization','Probability');
    prob_dist_sh(w,:)=histcounts(dist_sh,[0:6],'Normalization','Probability');
    
    delta_coupling_conc=[delta_coupling_conc,delta_coupling];
    
    for n=1:length(dist)
        count=count+1;
        index_count(count)=w;
        real_ens(count)=session{1,w}.real_ensemble(n);
        inf_ens(count)=session{1,w}.inf_ensemble(n);
        inf_ens_sh(count)=session{1,w}.inf_ensemble_sh(n,1);
    end
    
    clear dist dist_sh delta_coupling
    
end

for i=1:6
    [p(i),h(i),stats(i)]=ranksum(prob_dist(:,i),prob_dist_sh(:,i));
end

%Matrix after pooling cells
H=hist3([inf_ens',real_ens'],'Edges',{1:1:10 1:1:10});
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
caxis([0 250])
co.Label.String = 'Counts';

%Matrix after pooling cells - shuffle
H=hist3([inf_ens_sh',real_ens'],'Edges',{1:1:10 1:1:10});
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
caxis([0 250])
co.Label.String = 'Counts';


%For one session ; session=11
in_session11=find(index_count==11);

H=hist3([inf_ens(in_session11)',real_ens(in_session11)'],'Edges',{1:1:10 1:1:10});
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
caxis([0 20])
co.Label.String = 'Counts';

H=hist3([inf_ens_sh(in_session11)',real_ens(in_session11)'],'Edges',{1:1:10 1:1:10});
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
caxis([0 20])
co.Label.String = 'Counts';


%Probability of distance

mean_prob=mean(prob_dist);
sem_prob=std(prob_dist)./sqrt(15);
mean_prob_sh=mean(prob_dist_sh);
sem_prob_sh=std(prob_dist_sh)./sqrt(15);

figure
errorbar([0:5],mean_prob,sem_prob,'linewidth',2.5)
hold on
errorbar([0:5],mean_prob_sh,sem_prob_sh,'linewidth',2.5);
xticks([0,1,2,3,4,5])
ylabel('Probability');
xlabel('Distance between inferred and real ensemble');
set(gca,'fontsize',18);
box off
