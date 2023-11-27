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

N_sh=1000;
clus=10;
count=0;
prctile_th=99;

Radius_f_allsessions_vec=[];
MVL_allsessions_vec=[];
MI_allsessions_vec=[];
PR_all_sessions_step1=[];
PR_locked_all_sessions=[];

for w=19:length(PaS_sessions)
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
    file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];    
    load(file_name_anat,'-mat'); %Anatomical information
    %     load(file_name_dff,'-mat'); %DFF
    
    
    if isfield(dates,'actual_day')==1
        day=find(dates.actual_day==day);
    end
    
%     file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
%     file_name_dff=[dpath ['DFF_120ms_Do_SNRH','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

    load(file_name_spk,'-mat'); %Spike times
    spikes_d=full(spikes_d_s);
    [N,T]=size(spikes_d);
    %     dff=signal_dff;    
    
    
           
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Condition on having waves
    
    dt=big_table(row_w,8);
    
    if dt>0
        for i=1:N
            FRp(i,:)=full(smoothdata(spikes_d(i,:),'gaussian',1*dt)); %smooth using as kernel the dt chosen for each session
        end
        
        [coefft,scoret,~] = pca(FRp');
        phase_f=(atan2(smooth(scoret(:,2),floor(1*dt)),smooth(scoret(:,1),floor(1*dt))));
        radius_f=sqrt(coefft(:,2).*coefft(:,2)+coefft(:,1).*coefft(:,1));
    else
        FRp=spikes_d;
        
        [coefft,scoret,~] = pca(FRp');
        phase_f=(atan2(scoret(:,2),scoret(:,1)));
        radius_f=sqrt(coefft(:,2).*coefft(:,2)+coefft(:,1).*coefft(:,1));
    end
    
    %     for i=1:N
    %         FRp(i,:)=full(fire_rate(spikes_d(i,:),1*dt,'g')); %smooth using as kernel the dt chosen for each session
    %     end


    %We do not identify waves because there are no waves
    
%     num_clus_discr=10;
%     make_fig=0;
%     [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes_d);
%     
%     spikes_r=[]; %Reduced spike matrix; only contains wave epochs
%     phase_r=[]; %Reduced phase; only contains wave epochs
%     
%    
%     for wa=1:size(table_u,1)               
%         spikes_r=[spikes_r,spikes_d(:,table_u(wa,1):table_u(wa,2))];
%         phase_r=[phase_r;phase_f(table_u(wa,1):table_u(wa,2))];      
%     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    spikes=spikes_d;
    phase=phase_f;
%     phase=phase_r;
%     spikes_d=[];
%     spikes_d=spikes_r;
     phase_r=phase_f;   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MVL
    
    for i=1:N
        p=phase(find(spikes_d(i,:)));
        var_p(i)=circ_var(p);
        mean_p(i)=circ_mean(p);
        std_p(i)=circ_std(p);
        MVL(i) = circ_r(p);
        prob_phase_firing(i,:)=histcounts(p,-pi:2*pi/40:pi,'Normalization','Probability');
        clear p H
    end
    
    mean_p_all_sessions{w}=mean_p;
    std_p_all_sessions{w}=std_p;

    for i=1:N
        %disp(i)
        for sh=1:N_sh
            p_sh=phase(randperm(length(phase_r),length(find(spikes_d(i,:)))));
            var_p__sh(i,sh)=circ_var(p_sh);
            std_p__sh(i,sh)=circ_std(p_sh);
            mean_p_sh(i,sh)=circ_mean(p_sh);
            MVL_sh(i,sh) = circ_r(p_sh);
            
            clear H p_sh
        end
    end
    
    mean_p_sh_all_sessions{w}=mean_p_sh;
    std_p_sh_all_sessions{w}=std_p__sh;
    
    locking=MVL;
    locking_sh=MVL_sh;
    locking_sh_mean=mean(MVL_sh,2);
    locking_sh_99=prctile(MVL_sh,99,2);
    locking_sh_1=prctile(MVL_sh,1,2);
    
    MVL_allsessions{w}=MVL;
    MVL_allsessions_vec=[MVL_allsessions_vec,MVL];

    
    not_locked=find(locking<=locking_sh_99');
    locked=setdiff(1:N,not_locked);
    
    non_inf(w)=length(not_locked)/N;
    infor(w)=length(locked)/N;
    
    locked_all_sessions{w}=locked;
    not_locked_all_sessions{w}=not_locked;
    
    Radius_f_allsessions{w}=radius_f;
    Radius_f_allsessions_vec=[Radius_f_allsessions_vec;radius_f];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mutual information and Radius
    
    disc_phase=10; %Number of bins for the phase
%     disc_spk=4; %Number of bins for the spikes
    new_bin = 4;    
    
    %Downsample signals
    phase_down=downsample(phase,new_bin);
%     phase_di=discretize(phase_down,-pi:2*pi/disc_phase:pi);
    
    for i=1:length(phase_down)
        spk_do(:,i)=sum(spikes_d(:,(i-1)*new_bin+1:i*new_bin),2);
    end
    
    cells=1:N;
    for ind=1:N
      %  disp(ind)
        i=cells(ind);
        
        %Preprocess signal
        
        calc=spk_do(ind,:);
%         calc_di=discretize(calc,disc_spk);
        
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
    
    MI_allsessions{w}=MI;
    MI_allsessions_vec=[MI_allsessions_vec,MI];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distribution of cells in ensembles
    
    clear Ens
    % [~,sorting,~]=get_sorting(spikes_d);
    [~,sorting_w,~]=get_sorting_smoothed(spikes,dt);
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distance within and across ensemble, need to come after distribution of cells in ensembles
    
    pixel_size_new=1.18185;
    pixel_size_old=1.78211;
    
    if w>10
        pixel_size=pixel_size_old;
    else
        pixel_size=pixel_size_new;       
    end
    
    d_within_all_sessions=[];
    d_across_all_sessions=[];
    
    %Delta position on the tissue
    d_within=[];
    d_across=[];
    c_w=0;%counter within
    c_a=0;%counter across
    
    count=0;
    for i=1:N
        for j=i+1:N
            count=count+1;
            r_i=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
            r_j=[Anat.med{1,j}(1)*pixel_size,Anat.med{1,j}(2)*pixel_size];
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
%             delta_phase_mean_vec(count)=angdiff(pref_phase(i),pref_phase(j));
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
        
        d_within_ens(i,1,w)=mean(aux_within);
        d_within_ens(i,2,w)=std(aux_within);
        d_within_ens(i,3,w)=std(aux_within)./sqrt(length(aux_within));
        
        d_across_ens(i,1,w)=mean(aux_across);
        d_across_ens(i,2,w)=std(aux_across);
        d_across_ens(i,3,w)=std(aux_across)./sqrt(length(aux_across));

        
        d_within_all_sessions=[d_within_all_sessions,d_within];
        d_across_all_sessions=[d_across_all_sessions,d_across];
        
         clear within across aux_within aux_across
    end
      
    %%%%%%% Participation index is not computed for PaS
    
% % % %     spikes_sorted=spikes(sorting_w,:);
% % % %     new_mat=spikes_sorted;
% % % %     
% % % %     num_steps=6;
% % % %     step=0;
% % % %     mean_PR(1)=0;
% % % %     std_PR(1)=0.5;
% % % %     mean_PR_I(1)=0;
% % % %     std_PR_I(1)=0.5;
% % % %     n_waves=size(table_u,1);
% % % %     % for i=1:N
% % % %     %     new_mat(i,:)=spikes_d(sorting_angle(i),:);
% % % %     % end
% % % %     
% % % %     cell_wave_part=zeros(N,size(table_u,1));
% % % %     
% % % %     for j=1:num_steps
% % % %         
% % % %         step=step+1;
% % % %         for i=1:size(new_mat,1) %Loop on elements of mat
% % % %             calc=(new_mat(i,:)); %spike train of row i
% % % %             
% % % %             for wa=1:size(table_u,1) %Spikes of each cell (or row of new_mat) per wave
% % % %                 spikes_per_wave(i,wa)=sum(calc(table_u(wa,1):table_u(wa,2)));
% % % %             end
% % % %             
% % % %             aux=spikes_per_wave(i,:);
% % % %             aux2=find((cumsum(sort(aux,'descend'))./sum(aux))>0.90,1);%Number of waves needed to account for 95% of the spikes
% % % %             
% % % %             if j==1 %This indicayes the waves in which a neuron fired, based on the criterion of the waves for which 95% of spiking is captured
% % % %                 [ind1,ind2]=sort(aux,'descend');
% % % %                 cell_wave_part(i,ind2(1:aux2))=1;
% % % %                 clear ind1 ind2
% % % %             end
% % % %             
% % % %             if isempty(aux2)==1
% % % %                 aux2=0;
% % % %             end
% % % %             
% % % %             PR(i,step)=aux2/size(table_u,1);
% % % %             
% % % %             %Second measure of PI
% % % %             aux3=cumsum(sort(aux,'descend'))./sum(aux);
% % % %             auc=sum(aux3*1/n_waves);
% % % %             %         PR_I(i,step)=auc-0.5;
% % % %             PR_I(i,step)=1-auc;
% % % %             
% % % %             clear aux aux2 aux3 auc
% % % %             
% % % %         end
% % % %         
% % % %         mean_PR_all_sessions(step,w)=mean(nonzeros(PR(:,step)),1);
% % % %         std_PR_all_sessions(step,w)=std(nonzeros(PR(:,step)),[],1);
% % % %         sem_PR_all_sessions(step,w)=std(nonzeros(PR(:,step)),[],1)/sqrt(length(nonzeros(PR(:,step))));
% % % %  
% % % %         mean_PR_I_all_sessions(step,w)=mean(nonzeros(PR_I(:,step)),1);
% % % %         std_PR_I_all_sessions(step,w)=std(nonzeros(PR_I(:,step)),[],1);
% % % %         
% % % %         new_mat2=new_mat;
% % % %         clear new_mat;
% % % %         new_mat=coarse_graining(new_mat2);
% % % %         
% % % %         clear spikes_per_wave aux aux2
% % % %     end
% % % %     
% % % %     PR_all_sessions{w}=PR;
% % % %     PR_all_sessions_step1=[PR_all_sessions_step1;PR(:,1)];
% % % %     
% % % %     
% % % %     
% % % %     %PR is ordered according to sorting_w. The rest of the quantities, such as
% % % %     %MVL, is ordered in a "canonical" basis. In order to keep the PI of the
% % % %     %locked cells, we need to identify the position of the locked cells in
% % % %     %sorting_w, and keep only those. These positions will be given by
% % % %     %index_sorting_locked. Sorting_w_locked is the sorting_w vector without the
% % % %     %not_locked cells. The connection between these two is the following: If we
% % % %     %take sorting_w, and keep only the components given by
% % % %     %index_sorting_locked, we should get sorting_w_locked.
% % % %     
% % % %     count=0;
% % % %     sorting_copy=sorting_w;
% % % %     for j=1:length(not_locked)
% % % %         for i=1:N
% % % %             if sorting_copy(i)==not_locked(j)
% % % %                 sorting_copy(i)=0;
% % % %                 count=count+1;
% % % %                 index_sorting_not_locked(count)=i; %Components of the sorted vector that are *not* locked
% % % %             end
% % % %         end
% % % %     end
% % % %     index_sorting_locked=1:N;
% % % %     index_sorting_locked(index_sorting_not_locked)=[]; %Components of the sorted vector that are locked
% % % %     
% % % %     sorting_w_not_locked=sorting_w(find(sorting_copy ==0)); %sorting_w without cells that are  locked
% % % %     
% % % %     sorting_copy(sorting_copy ==0)=[];
% % % %     sorting_w_locked=sorting_copy; %sorting_w without cells that are *not* locked
% % % %     
% % % %     PR_locked=PR(:,1);
% % % %     PR_not_locked=PR(:,1);
% % % %     
% % % %     PR_locked=PR_locked(index_sorting_locked);
% % % %     PR_not_locked=PR_not_locked(index_sorting_not_locked);
% % % %     
% % % %     PR_locked_all_sessions=[PR_locked_all_sessions;PR_locked];
% % % %     sorting_all_sessions{w}=sorting_w;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Locking to ensemble
% % % %     N_sh_locking_ensemble=200;
% % % %     
% % % %     %Downsample and binarize
% % % %     new_bin = 4;
% % % %     num_bins=floor(size(spikes_d,2)/new_bin);
% % % %     for i=1:num_bins %Length of phase conditioned on waves and previously downsampled to 4 bins (0.5 bin size) when computing MI
% % % %         sp_do(:,i)=sum(spikes_d(:,(i-1)*new_bin+1:i*new_bin),2);
% % % %     end    
% % % %     for n=1:N
% % % %         sp_do(n,find(sp_do(n,:)>0))=1;
% % % %     end
% % % %     
% % % %     %Calculate coactivity
% % % %     prob_k=rcoact_dist(sp_do);
% % % %     
% % % %     for n=1:N %Loop on neurons
% % % %         p_1=find(sp_do(n,:)>0); %Frames in which the neuron fired
% % % %         p_0=find(sp_do(n,:)==0); %Frames in which the neuron *did not* fired
% % % %         mean_fr=(length(p_1)/size(sp_do,2)); %Mean firing rate of cell "n"
% % % %         
% % % %        % disp(n)
% % % %         for ense=1:10 %Loop on ensembles
% % % %             cells_ens_aux=find(Ens(:,2)==ense);
% % % %             cells_ens=Ens(cells_ens_aux,1); %Cells that belong to ensemble ense
% % % %             cells=setdiff(cells_ens,n); %Removes from cells_ens cell *n* for which we are computing the tuning curve
% % % %             a=sum(sp_do(cells,:)); %Number of cells that fired at each time point in ensemble "ense"
% % % %             prob_ki=rcoact_dist(sp_do(cells,:)); %Probability of coactivity for the cells in ensemble "ense"
% % % %             
% % % %             %Computing the sensitivity and the tuning curve
% % % %             for k=0:N
% % % %                 aux_k=find(a==k);
% % % %                 if length(aux_k)>0
% % % %                     TC2(n,k+1,ense)=length(intersect(p_1,aux_k))./length(aux_k);
% % % %                     TC3(n,k+1,ense)=(length(intersect(p_1,aux_k))./length(aux_k))/mean_fr;
% % % %                     %                 Arg_sens2(n,k+1,ense)= TC2(n,k+1,ense)* TC2(n,k+1,ense)/(length(p_1)/size(sp_do,2))*prob_ki(k+1);
% % % %                     Arg_sens2(n,k+1,ense)= TC2(n,k+1,ense)* TC2(n,k+1,ense)*prob_ki(k+1);
% % % %                     
% % % %                 else
% % % %                     TC2(n,k+1,ense)=0;
% % % %                     TC3(n,k+1,ense)=0;
% % % %                     Arg_sens2(n,k+1,ense)=0;
% % % %                 end
% % % %                 
% % % %                 clear aux_k
% % % %             end
% % % %             
% % % %             %         Sens(n,ense) = sqrt(nansum(Arg_sens2(n,:,ense))- (length(p_1)/size(sp_do,2))*(length(p_1)/size(sp_do,2)));
% % % %             Sens(n,ense) = sqrt(nansum(Arg_sens2(n,:,ense))- (mean_fr*mean_fr));
% % % %             AUC(n,ense)=sum(TC2(n,:,ense));
% % % %             
% % % %             clear aux2 aux
% % % %         end
% % % %         
% % % %         [c,d]=max(Sens(n,:)); %c is the maximum sensitivity. d is the ensemble that maximizes the sensitivity
% % % %         [a,b]=max(AUC(n,:));  %a is the maximum AUC. b is the ensemble that maximizes the AUC
% % % %         inferred_ensemble_allsessions{w,n}=b; %Inferred ensemble using the AUC of the tuning curve
% % % %         inferred_ensemble_sens_allsessions{w,n}=d; %Inferred ensemble using the sensitivity
% % % %         aux3=find(Ens(:,1)==n);
% % % %         real_ensemble_allsessions{w,n}=Ens(aux3,2);
% % % %         
% % % %         clear p_1 p_0 cells_ens_aux cells_ens cells a prob_ki
% % % %     end
% % % %     
% % % %     %Shuffle spike times of one cell at the time
% % % %     for n=1:N
% % % %         for sh=1:N_sh_locking_ensemble
% % % %             p_1=randperm(size(sp_do,2),length(find(sp_do(n,:)>0))); %Frames in which the neuron fired
% % % %             mean_fr=(length(p_1)/size(sp_do,2)); %Mean firing rate of cell "n"
% % % %            % disp(n)
% % % %             for ense=1:10
% % % %                 cells_ens_aux=find(Ens(:,2)==ense);
% % % %                 cells_ens=Ens(cells_ens_aux,1);
% % % %                 cells=setdiff(cells_ens,n);
% % % %                 a=sum(sp_do(cells,:));  %Number of cells that fired at each time point in ensemble "ense"
% % % %                 prob_ki=rcoact_dist(sp_do(cells,:));
% % % %                 
% % % %                 %Computing the sensitivity and the tuning curve
% % % %                 for k=0:N
% % % %                     aux_k=find(a==k);
% % % %                     if length(aux_k)>0
% % % %                         TC2_sh(n,k+1,ense)=length(intersect(p_1,aux_k))./length(aux_k);
% % % %                         TC3_sh(n,k+1,ense)=(length(intersect(p_1,aux_k))./length(aux_k))/mean_fr;
% % % %                         %                     Arg_sens2_sh(n,k+1,ense)= TC2_sh(n,k+1,ense)* TC2_sh(n,k+1,ense)/(length(p_1)/size(sp_do,2))*prob_ki(k+1);
% % % %                         Arg_sens2_sh(n,k+1,ense)= TC2_sh(n,k+1,ense)* TC2_sh(n,k+1,ense)*prob_ki(k+1);
% % % %                     else
% % % %                         TC2_sh(n,k+1,ense)=0;
% % % %                         TC3_sh(n,k+1,ense)=0;
% % % %                         Arg_sens2_sh(n,k+1,ense)=0;
% % % %                     end
% % % %                     clear aux_k
% % % %                 end
% % % %                 Sens_sh(n,ense,sh) = sqrt(nansum(Arg_sens2_sh(n,:,ense))- (length(p_1)/size(sp_do,2))*(length(p_1)/size(sp_do,2)));
% % % %                 AUC_sh(n,ense,sh)=sum(TC2_sh(n,:,ense));
% % % %                 clear aux2 aux
% % % %             end
% % % %             
% % % %             [c,d]=max(Sens_sh(n,:,sh));
% % % %             [a,b]=max(AUC_sh(n,:,sh));
% % % %             inferred_ensemble_sh_allsessions{w}(n,sh)=b;
% % % %             inferred_ensemble_sens_sh_allsessions{w}(n,sh)=d;
% % % %             
% % % %             clear p_1 p_0 cells_ens_aux cells_ens cells a prob_ki prob_k
% % % %         end
% % % %     end
    
    clear cells cells_di cells_d coeff dff exclude1 exclude2 latent_pca MI_TB phase phase_d phase_di score signal_dff snr SNR spikes spikes_d spikes_d_s p MVL
    clear not_locked locked locking locking_sh locking_sh_mean locking_sh_99 locking_sh_1 MVL_sh p_sh sp_do phase phase_d mean_p mean_p_sh prob_phase_firing
    clear std_p std_p__sh var_p var_p__sh calc_di coefft FRp MI MI_b MI_withbias phase_down phase_f phase_r radius_f scoret spikes_r spk_do ...
          table_u Anat d_within d_across delta_tissue delta_tissue_vec ens Ens ense_n  sorting_w r_i r_j PR new_mat new_mat2 spikes_sorted ...
          PR_locked PR_not_locked index_sorting_locked sorting_w_not_locked sorting_w_locked delta_phase_mean index_sorting_not_locked...
          Sens_sh AUC_sh TC2_sh TC3_sh Arg_sens2_sh Sens AUC TC2 TC3 Arg_sens2 p_1 p_0 mean_fr calc cell_wave_part prob_k sp_do ind
end

%% Figures fraction of cells locked to waves

figure
bar(1,mean(infor)*100,'FaceColor', [170 170 170]/255);
hold on
bar(2,mean(non_inf)*100,'FaceColor', [170 170 170]/255);
scatter(ones(1,length(infor)),infor*100,'k');
hold on
scatter(2*ones(1,length(non_inf)),non_inf*100,'k');
axis([0 3 0 100])
xticks([1,2]);
xticklabels({'Locked','Not locked'})
ylabel('Cells %')
set(gca,'fontsize',16);
box off

figure
b=bar([1,2],[mean(infor*100),mean(non_inf*100)],0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];%[ 0.4703    0.1098    0.4286]*1.1;
b.CData(2,:) = [1 0 0];%'r';%[0.0704    0.7457    0.7258];
hold on
er=errorbar([1,2],[mean(infor*100),mean(non_inf*100)],[std(infor*100)/(sqrt(15)),std(non_inf*100)/(sqrt(15))]);
er.LineStyle='none';
er.Color='k';
er.LineWidth=2.5;
axis([0.5,2.5,0,101])
ylabel('Cells %')
xticks([1,2])
xticklabels({'Locked', 'Not locked'});
box off
set(gca,'fontsize',26)
yticks([0 50 100])
[H_locking,P_locking,ci_locking,stat_locking]=ttest((infor*100-non_inf*100));


%% Figures of MVL, Radius and MI pooling all cells

MI_not_locked=[];
MVL_not_locked=[];
radius_not_locked=[];

for i=1:length(waves)
    mi_per_session=MI_allsessions{i};
    radius_per_session=Radius_f_allsessions{i};
    MVL_per_session=MVL_allsessions{i};
    
    MI_not_locked=[MI_not_locked,mi_per_session(not_locked_all_sessions{i})];
    radius_not_locked=[radius_not_locked;radius_per_session(not_locked_all_sessions{i})];
    MVL_not_locked=[MVL_not_locked,MVL_per_session(not_locked_all_sessions{i})];

    clear mi_per_session radius_per_session MVL_per_session
end

figure
scatter(MI_allsessions_vec,MVL_allsessions_vec,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255)
alpha 0.5
ylabel('Locking to phase');
xlabel('MI corrected for bias (bits)');
set(gca,'fontsize',18)
axis([-0.02 0.35 0 1])
box off
hold on
scatter(MI_not_locked,MVL_not_locked,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor','r')
alpha 0.6
linear_model_MI_MVL = fitlm(MI_allsessions_vec,MVL_allsessions_vec);
[R,P] = corrcoef(MI_allsessions_vec,MVL_allsessions_vec);


figure
hold on
scatter(Radius_f_allsessions_vec,MI_allsessions_vec,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255)
alpha 0.5
ylabel('MI corrected for bias (bits)');
xlabel('Radius PC1-PC2')
set(gca,'fontsize',20)
axis([0 0.4 0 0.4])
scatter(radius_not_locked,MI_not_locked,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor','r')
alpha 0.6
linear_model_MI_Rad = fitlm(Radius_f_allsessions_vec,MI_allsessions_vec);
[R,P] = corrcoef(Radius_f_allsessions_vec,MI_allsessions_vec);

figure
hold on
scatter(Radius_f_allsessions_vec,MVL_allsessions_vec,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor',[128,128,128]./255)
alpha 0.5
ylabel('Locking to phase');
xlabel('Radius PC1-PC2')
set(gca,'fontsize',20)
axis([0 0.35 0 1])
scatter(radius_not_locked,MVL_not_locked,40,'filled','MarkerEdgeColor',[32,32,32]./255,'MarkerFaceColor','r')
alpha 0.6
linear_model_MVL_Rad = fitlm(Radius_f_allsessions_vec,MVL_allsessions_vec);
[R,P] = corrcoef(Radius_f_allsessions_vec,MVL_allsessions_vec);
%% PI of all cells, no filtering for locked ones

%Histogram
figure
H=histogram(PR_all_sessions_step1,0:0.05:1);
x=H.BinEdges(2:end);
y=H.Values;
bar(x,y,'FaceColor',[137, 137, 255]./255,'Linewidth',1)
ylim([0 800]);
ylabel('Counts');
xlabel('PI');
set(gca,'fontsize',18)
box off

%Pooling cells

PR_allsessions_1=[];
PR_allsessions_2=[];
PR_allsessions_3=[];
PR_allsessions_4=[];
PR_allsessions_5=[];
PR_allsessions_6=[];

for i=1:length(waves)
    PR_allsessions_1=[PR_allsessions_1;nonzeros(PR_all_sessions{1,i}(:,1))];
    PR_allsessions_2=[PR_allsessions_2;nonzeros(PR_all_sessions{1,i}(:,2))];
    PR_allsessions_3=[PR_allsessions_3;nonzeros(PR_all_sessions{1,i}(:,3))];
    PR_allsessions_4=[PR_allsessions_4;nonzeros(PR_all_sessions{1,i}(:,4))];
    PR_allsessions_5=[PR_allsessions_5;nonzeros(PR_all_sessions{1,i}(:,5))];
    PR_allsessions_6=[PR_allsessions_6;nonzeros(PR_all_sessions{1,i}(:,6))];        
end

mean_PR(1)=mean(PR_allsessions_1);
mean_PR(2)=mean(PR_allsessions_2);
mean_PR(3)=mean(PR_allsessions_3);
mean_PR(4)=mean(PR_allsessions_4);
mean_PR(5)=mean(PR_allsessions_5);
mean_PR(6)=mean(PR_allsessions_6);


std_PR(1)=std(PR_allsessions_1);
std_PR(2)=std(PR_allsessions_2);
std_PR(3)=std(PR_allsessions_3);
std_PR(4)=std(PR_allsessions_4);
std_PR(5)=std(PR_allsessions_5);
std_PR(6)=std(PR_allsessions_6);

sem_PR(1)=std(PR_allsessions_1)/sqrt(length(PR_allsessions_1));
sem_PR(2)=std(PR_allsessions_2)/sqrt(length(PR_allsessions_2));
sem_PR(3)=std(PR_allsessions_3)/sqrt(length(PR_allsessions_3));
sem_PR(4)=std(PR_allsessions_4)/sqrt(length(PR_allsessions_4));
sem_PR(5)=std(PR_allsessions_5)/sqrt(length(PR_allsessions_5));
sem_PR(6)=std(PR_allsessions_6)/sqrt(length(PR_allsessions_6));
        
figure
errorbar(1:6,mean_PR,sem_PR,'k-o','Linewidth',2);
ylabel('PI');
xlabel('Merging steps');
axis([0.5 6.5 0 1])
xticks([1 2 3 4 5 6])
set(gca,'fontsize',18)
box off

% Per session
for i=1:length(waves)
    PR_persessions_1(i)=mean(nonzeros(PR_all_sessions{1,i}(:,1)));
    PR_persessions_2(i)=mean(nonzeros(PR_all_sessions{1,i}(:,2)));
    PR_persessions_3(i)=mean(nonzeros(PR_all_sessions{1,i}(:,3)));
    PR_persessions_4(i)=mean(nonzeros(PR_all_sessions{1,i}(:,4)));
    PR_persessions_5(i)=mean(nonzeros(PR_all_sessions{1,i}(:,5)));
    PR_persessions_6(i)=mean(nonzeros(PR_all_sessions{1,i}(:,6)));    
end

PR_persessions_mean(1)=mean(PR_persessions_1);
PR_persessions_mean(2)=mean(PR_persessions_2);
PR_persessions_mean(3)=mean(PR_persessions_3);
PR_persessions_mean(4)=mean(PR_persessions_4);
PR_persessions_mean(5)=mean(PR_persessions_5);
PR_persessions_mean(6)=mean(PR_persessions_6);

PR_persessions_sem(1)=std(PR_persessions_1)/sqrt(length(PR_persessions_1));
PR_persessions_sem(2)=std(PR_persessions_2)/sqrt(length(PR_persessions_2));
PR_persessions_sem(3)=std(PR_persessions_3)/sqrt(length(PR_persessions_3));
PR_persessions_sem(4)=std(PR_persessions_4)/sqrt(length(PR_persessions_4));
PR_persessions_sem(5)=std(PR_persessions_5)/sqrt(length(PR_persessions_5));
PR_persessions_sem(6)=std(PR_persessions_6)/sqrt(length(PR_persessions_6));

figure
errorbar(1:6,PR_persessions_mean,PR_persessions_sem,'k-o','Linewidth',2);
ylabel('PI');
xlabel('Merging steps');
axis([0.5 6.5 0 1])
xticks([1 2 3 4 5 6])
set(gca,'fontsize',18)
box off

[p_PI,tbl_PI,stats_PI] = anova1([PR_persessions_1',PR_persessions_2',PR_persessions_3',PR_persessions_4',PR_persessions_5',PR_persessions_6']);
%% Distance within and across ensembles

%Per session

for w=1:length(waves)
mean_within(w)=mean(d_within_ens(:,1,w));
mean_across(w)=mean(d_across_ens(:,1,w));

end

total_mean_within=mean(mean_within);
total_sem_within=std(mean_within)/sqrt(length(waves));
total_mean_across=mean(mean_across);
total_sem_across=std(mean_across)/sqrt(length(waves));

figure
b=bar([1,2],[total_mean_within,total_mean_across],0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [0.2422    0.1504    0.6603]*1.1;
b.CData(2,:) = [0.0704    0.7457    0.7258];
hold on
er=errorbar([1,2],[total_mean_within,total_mean_across],[total_sem_within,total_sem_across]);
er.LineStyle='none';
er.Color='k';
er.LineWidth=1.5;
axis([0.5,2.5,0,350])
ylabel('Distance on tissue (\mum)')
xticks([1,2])
xticklabels({'Within ensembles', 'Across ensembles'});
box off
set(gca,'fontsize',16)
[h_dist,p_dist,ci_dist,stats_dist] = ttest2(mean_within,mean_across);
[h_dist,p_dist,stats_dist] = ranksum(mean_within,mean_across);


%Pooling all distances
figure
b=bar([1,2],[mean(d_within_all_sessions),mean(d_across_all_sessions)],0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [0.2422    0.1504    0.6603]*1.1;
b.CData(2,:) = [0.0704    0.7457    0.7258];
hold on
er=errorbar([1,2],[mean(d_within_all_sessions),mean(d_across_all_sessions)],[std(d_within_all_sessions),std(d_across_all_sessions)]);
er.LineStyle='none';
er.Color='k';
er.LineWidth=1.5;
axis([0.5,2.5,0,350])
ylabel('Distance on tissue (\mum)')
xticks([1,2])
xticklabels({'Within ensembles', 'Across ensembles'});
box off
set(gca,'fontsize',16)


%% Continuity of mean phase

figure
for w=1:length(waves)
    mean_p=mean_p_all_sessions{w};
    mean_p_locked=mean_p(locked_all_sessions{w});
    std_p=std_p_all_sessions{w};
    std_p_locked=std_p(locked_all_sessions{w});
    locked=locked_all_sessions{w};
    [~,sorting_mean_angle]=sort(mean_p_locked,'ascend');
    N_n=size(locked,2);
    
    
    subplot(4,4,w)
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
    ylabel({'Preferred phase';'(rad)'})
    xlabel('Neuron #')
    set(gca,'fontsize',18)
    box off
    
    clear meap_p mean_p_locked std_p_locked std_p locked locked sorting_mean_angle
end


 %Quantification only with locked cells
 
for w=1:length(waves)
    locked=locked_all_sessions{w};
    mean_p=mean_p_all_sessions{w};
    mean_p_locked=mean_p(locked);
    Y=histcounts(mean_p_locked,-pi:2*pi/10:pi,'Normalization','Probability');
%     N=length(locked);
        
    p_max=(1/10)*ones(1,10);
    H_max= -sum(p_max .* log(p_max));
    
    p_real= Y;
    H_real= -sum(p_real .* log(p_real));
    
    H_ratio(w) = H_real/H_max;    
    
    clear mean_p Y locked mean_p mean_p_locked
end

for w=1:length(waves)
    locked=locked_all_sessions{w};
    for sh=1:N_sh
        mean_p_sh=mean_p_sh_all_sessions{w}(:,sh);
        mean_p_sh_locked=mean_p_sh(locked);

        Y=histcounts(mean_p_sh,-pi:2*pi/10:pi,'Normalization','Probability');
%         N=length(sorting_all_sessions{w});
    
        p_max=(1/10)*ones(1,10);
        H_max= -sum(p_max .* log(p_max));
    
        p_real= Y;
        H_real= -sum(nonzeros(p_real) .* log(nonzeros(p_real)));
    
        H_ratio_sh(w,sh) = H_real/H_max;    

        clear mean_p_sh Y mean_p_sh_locked
    end
    clear locked
end

figure
errorbar(H_ratio,mean(H_ratio_sh,2),std(H_ratio_sh,[],2)/sqrt(N_sh),'*','Linewidth',1.5);
axis([0.2 1.0 0.2 1]);
r=refline(1,0);
r.LineStyle='--';
r.Color='k';
r.LineWidth=2;
set(gca,'fontsize',18);
ylabel('H_r_a_t_i_o - Shuffle');
xlabel('H_r_a_t_i_o - Data');
xticks([0.2 0.5 1])
yticks([0.2 0.5 1]);
box off

figure
errorbar(H_ratio,mean(H_ratio_sh,2),std(H_ratio_sh,[],2),'*','Linewidth',1.5);
axis([0.85 1.0 0.2 1]);
r=refline(1,0);
r.LineStyle='--';
r.Color='k';
r.LineWidth=2;
set(gca,'fontsize',18);
ylabel('H_r_a_t_i_o - Shuffle');
xlabel('H_r_a_t_i_o - Data');
xticks([0.85 1])
yticks([0.2 0.5 1]);
box off

figure
hold on
bar([1,2],[mean(H_ratio),mean(mean(H_ratio_sh,2))])
aux=mean(H_ratio_sh,2);
errorbar([1,2],[mean(H_ratio),mean(mean(H_ratio_sh,2))],[std(H_ratio)/sqrt(15),...
    std(mean(H_ratio_sh,2))/sqrt(15)],'.k','linewidth',2);
ylabel('H_r_a_t_i_o');
xticks([1,2]);
xticklabels({'Data','Shuffle'});
set(gca,'fontsize',16);

[p_h,t_h,stats]=ranksum(H_ratio,mean(H_ratio_sh,2));
% hold on
% for i=1:length(waves)
% scatter([[1,2]],[H_ratio(i),aux(i)],'k');
% end

 %Quantification with all cells
%  
% for w=1:length(waves)
%     mean_p=mean_p_all_sessions{w};
%     Y=histogram(mean_p,-pi:2*pi/10:pi);
%     N=length(sorting_all_sessions{w});
%     
%     
%     p_max=(1/10)*ones(1,10);
%     H_max= -sum(p_max .* log(p_max));
%     
%     p_real= Y.Values./N;
%     H_real= -sum(p_real .* log(p_real));
%     
%     H_ratio(w) = H_real/H_max;    
%     
%     clear mean_p Y
% end
% 
% for w=1:length(waves)
%     for sh=1:N_sh
%     mean_p_sh=mean_p_sh_all_sessions{w}(:,sh);
%     Y=histogram(mean_p_sh,-pi:2*pi/10:pi);
%     N=length(sorting_all_sessions{w});
%     
%     p_max=(1/10)*ones(1,10);
%     H_max= -sum(p_max .* log(p_max));
%     
%     p_real= Y.Values./(N);
%     H_real= -sum(nonzeros(p_real) .* log(nonzeros(p_real)));
%     
%     H_ratio_sh(w,sh) = H_real/H_max;    
% 
%     clear mean_p_sh Y 
%     end
% end
% 
% figure
% errorbar(H_ratio,mean(H_ratio_sh,2),std(H_ratio_sh,[],2),'*','Linewidth',1.5);
% axis([0.89 1.01 0.2 1]);
% r=refline(1,0);
% r.LineStyle='--';
% r.Color='k';
% r.LineWidth=2;
% set(gca,'fontsize',18);
% ylabel('H_r_a_t_i_o - Shuffle');
% xlabel('H_r_a_t_i_o - Data');
% xticks([0.9 1])
% yticks([0.2 0.6 1]);
% 
% figure
% bar([1,2],[mean(H_ratio),mean(mean(H_ratio_sh,2))])
%     aux=mean(H_ratio_sh,2);
% hold on
% for i=1:13
% scatter([[1,2]],[H_ratio(i),aux(i)],'k');
% end


