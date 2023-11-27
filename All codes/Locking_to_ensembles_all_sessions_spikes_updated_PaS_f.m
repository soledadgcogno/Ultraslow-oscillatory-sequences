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
for w=5:length(repeated_days)
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

%% Calculates locking

N_sh_locking_ensemble=500;

clus=10;
count=0;
prctile_th=99;

for w=1:length(PaS_sessions)
    row_w=PaS_sessions(w);
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
    if isfield(dates,'actual_day')==1
        day=find(dates.actual_day==day);
    end
    
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Locking to ensemble
    
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
    TC2=zeros(N,N+1,10);
    for n=1:N %Loop on neurons
        p_1=find(sp_do(n,:)>0); %Frames in which the neuron fired
        p_0=find(sp_do(n,:)==0); %Frames in which the neuron *did not* fired
        mean_fr=(length(p_1)/size(sp_do,2)); %Mean firing rate of cell "n"
        
       % disp(n)
        for ense=1:10 %Loop on ensembles
            cells_ens_aux=find(Ens(:,2)==ense);
            cells_ens=Ens(cells_ens_aux,1); %Cells that belong to ensemble ense
            cells=setdiff(cells_ens,n); %Removes from cells_ens cell *n* for which we are computing the tuning curve
            a=sum(sp_do(cells,:)); %Number of cells that fired at each time point in ensemble "ense"
            prob_ki=rcoact_dist(sp_do(cells,:)); %Probability of coactivity for the cells in ensemble "ense"
            
            %Computing the sensitivity and the tuning curve
            for k=0:max(a)
                aux_k=find(a==k);
                if length(aux_k)>0
                    TC2(n,k+1,ense)=length(intersect(p_1,aux_k))./length(aux_k);
%                     TC3(n,k+1,ense)=(length(intersect(p_1,aux_k))./length(aux_k))/mean_fr;
                    %                 Arg_sens2(n,k+1,ense)= TC2(n,k+1,ense)* TC2(n,k+1,ense)/(length(p_1)/size(sp_do,2))*prob_ki(k+1);
                    Arg_sens2(n,k+1,ense)= TC2(n,k+1,ense)* TC2(n,k+1,ense)*prob_ki(k+1);                    
                else
                    TC2(n,k+1,ense)=0;
%                     TC3(n,k+1,ense)=0;
                    Arg_sens2(n,k+1,ense)=0;
                end                
                clear aux_k
            end
            
            %         Sens(n,ense) = sqrt(nansum(Arg_sens2(n,:,ense))- (length(p_1)/size(sp_do,2))*(length(p_1)/size(sp_do,2)));
            Sens(n,ense) = sqrt(nansum(Arg_sens2(n,:,ense))- (mean_fr*mean_fr));
%             AUC(n,ense)=sum(TC2(n,:,ense));
            
            clear aux2 aux
        end
        
        [c,d]=max(Sens(n,:)); %c is the maximum sensitivity. d is the ensemble that maximizes the sensitivity
%         [a,b]=max(AUC(n,:));  %a is the maximum AUC. b is the ensemble that maximizes the AUC
%         inferred_ensemble_allsessions{w,n}=b; %Inferred ensemble using the AUC of the tuning curve
        inferred_ensemble_sens_allsessions{w,n}=d; %Inferred ensemble using the sensitivity
        aux3=find(Ens(:,1)==n);
        real_ensemble_allsessions{w,n}=Ens(aux3,2);
        
        clear p_1 p_0 cells_ens_aux cells_ens cells a prob_ki
    end
    
    %Shuffle spike times of one cell at the time
    for n=1:N
        for sh=1:N_sh_locking_ensemble
            p_1=randperm(size(sp_do,2),length(find(sp_do(n,:)>0))); %Frames in which the neuron fired
            mean_fr=(length(p_1)/size(sp_do,2)); %Mean firing rate of cell "n"
           % disp(n)
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
%                 AUC_sh(n,ense,sh)=sum(TC2_sh(n,:,ense));
                clear aux2 aux
            end
            
            [c,d]=max(Sens_sh(n,:,sh));
%             [a,b]=max(AUC_sh(n,:,sh));
%             inferred_ensemble_sh_allsessions{w}(n,sh)=b;
            inferred_ensemble_sens_sh_allsessions{w}(n,sh)=d;
            
            clear p_1 p_0 cells_ens_aux cells_ens cells a prob_ki prob_k
        end
    end
    
    disp(w)
    
    clear cells cells_di cells_d coeff dff exclude1 exclude2 latent_pca MI_TB phase phase_d phase_di score signal_dff snr SNR spikes spikes_d spikes_d_s p MVL
    clear not_locked locked locking locking_sh locking_sh_mean locking_sh_99 locking_sh_1 MVL_sh p_sh sp_do phase phase_d mean_p mean_p_sh prob_phase_firing
    clear std_p std_p__sh var_p var_p__sh calc_di coefft FRp MI MI_b MI_withbias phase_down phase_f phase_r radius_f scoret spikes_r spk_do ...
          table_u Anat d_within d_across delta_tissue delta_tissue_vec ens Ens ense_n  sorting_w r_i r_j PR new_mat new_mat2 spikes_sorted ...
          PR_locked PR_not_locked index_sorting_locked sorting_w_not_locked sorting_w_locked delta_phase_mean index_sorting_not_locked...
          Sens_sh AUC_sh TC2_sh TC3_sh Arg_sens2_sh Sens AUC TC2 TC3 Arg_sens2 p_1 p_0 mean_fr calc cell_wave_part prob_k sp_do ind
end



%% Inferred ensembles


%Number of cells per session
for i=1:length(PaS_sessions)
    N_sessions(i)=size(inferred_ensemble_sens_sh_allsessions{1,i},1);
end

count=0;
for i=1:length(PaS_sessions)
    clear dist dist_sh
    for n=1:N_sessions(i)
        count=count+1;
        session_index(count)=i;
        ind_session(count)=i;
        inferred_ens(count)=inferred_ensemble_sens_allsessions{i,n};
        real_ens(count)=real_ensemble_allsessions{i,n};
        inferred_ens_sh(count,:)=inferred_ensemble_sens_sh_allsessions{1,i}(n,:);
        
        null_corrected=inferred_ens_sh(count,:)-real_ens(count);
        
        dist(n)=abs(inferred_ens(count)-real_ens(count));
        if dist(n)==9 dist(n)=1;
            elseif dist(n)==8 dist(n)=2;
            elseif dist(n)==7 dist(n)=3;
            elseif dist(n)==6 dist(n)=4;
        end
        
        for sh=1:500
            dist_sh(n,sh)=abs(inferred_ensemble_sens_sh_allsessions{1,i}(n,sh)-real_ens(count)');
            if dist_sh(n,sh)==9 dist_sh(n,sh)=1;
                elseif dist_sh(n,sh)==8 dist_sh(n,sh)=2;
                elseif dist_sh(n,sh)==7 dist_sh(n,sh)=3;
                elseif dist_sh(n,sh)==6 dist_sh(n,sh)=4;
            end
        end


    end
    
    session{1,i}.dist=dist;
    session{1,i}.dist_sh=dist_sh;

end


%Probability of distance

for i=1:length(PaS_sessions)
    prob_dist(i,:)=histcounts(session{1,i}.dist,[0:6],'Normalization','Probability');
    prob_dist_sh(i,:)=histcounts(session{1,i}.dist_sh,[0:6],'Normalization','Probability');
end

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



%Number of cells per session
for i=1:length(waves)
    N_sessions(i)=length(sorting_all_sessions{i});
end

count=0;
for i=1:length(waves)
    for n=1:N_sessions(i)
        count=count+1;
        ind_session(count)=i;
        inferred_ens(count)=inferred_ensemble_sens_allsessions{i,n};
        real_ens(count)=real_ensemble_allsessions{i,n};
        inferred_ens_sh(count,:)=inferred_ensemble_sens_sh_allsessions{1,i}(n,:);
        
        null_corrected=inferred_ens_sh(count,:)-real_ens(count);
    end
end


figure
scatter(inferred_ens,real_ens,'o','filled');
axis([0 11 0 11])
alpha 0.1
xlabel('Inferred Ensemble - Sens')
ylabel('Real Ensemble')
set(gca,'fontsize',16)
hold on 
axis square
xticks([1 5 10])
yticks([1 5 10])

H=hist3([inferred_ens',real_ens'],'Edges',{1:1:10 1:1:10});
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


figure
scatter(inferred_ens_sh(:,3),real_ens,'o','filled');
axis([0 11 0 11])
alpha 0.1
xlabel('Inferred Ensemble - Sens')
ylabel('Real Ensemble')
set(gca,'fontsize',16)
hold on 
axis square
xticks([1 5 10])
yticks([1 5 10])

H=hist3([inferred_ens_sh(:,3),real_ens']);
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
caxis([0 250])
co.Label.String = 'Counts';

% Now I correct for the distances

correct=find(inferred_ens==real_ens);
accuracy=length(correct)/sum(N_sessions)   ;
dist=abs(inferred_ens-real_ens);

for n=1:sum(N_sessions)   
   if dist(n)==9 dist(n)=1;
   elseif dist(n)==8 dist(n)=2;
   elseif dist(n)==7 dist(n)=3;
   elseif dist(n)==6 dist(n)=4;
   end       
end

MSE = sum(dist.*dist)/sum(N_sessions);

for sh=1:200
    dist_sh(:,sh)=abs(inferred_ens_sh(:,sh)-real_ens');
end

dist_sh_copy=dist_sh;

for sh=1:200    
    for n=1:sum(N_sessions)        
        if dist_sh_copy(n,sh)==9 dist_sh_copy(n,sh)=1;
        elseif dist_sh_copy(n,sh)==8 dist_sh_copy(n,sh)=2;
        elseif dist_sh_copy(n,sh)==7 dist_sh_copy(n,sh)=3;
        elseif dist_sh_copy(n,sh)==6 dist_sh_copy(n,sh)=4;
        end
    end
    MSE_sh(sh)=sum(dist_sh_copy(:,sh).*dist_sh_copy(:,sh))/sum(N_sessions);
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

for i=1:length(waves)
    idx=find(ind_session==i);
    
    sh_pooled=dist_sh_copy(idx,:);
end

%single cell
for i=1:sum(N_sessions)
   thr_5(i)=prctile(dist_sh_copy(i,:),5) ;
end

figure
histogram(dist_sh_copy(3,:))