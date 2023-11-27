%% Frequency of single cell oscillations - All sessions

clear all
close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
[big_table, waves] = get_big_table();


% Calculates locking

N_sh=200;
make_fig=0;
num_clus_discr=10;
clus=10;
count_s=0;
for w=1:length(waves)
    row_w=waves(w);
    disp(w)
    count_s=count_s+1;
    mouse=['L',num2str(big_table(row_w,1)),'M',num2str(big_table(row_w,2))];
    day=big_table(row_w,3);
    s=big_table(row_w,4);
    munit=big_table(row_w,5);
    clus=10;
    disc_phase=10;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    file_name_dff=[dpath ['DFF_120ms_Do_SNRH','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

    load(file_name_spk,'-mat'); %Spike times
    spikes=full(spikes_d_s);
    [N,~]=size(spikes);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Condition on having waves
    dt=big_table(row_w,8);
    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);
    spikes_w=[];
    for i=1:size(table_u,1)
        spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
    end
    T=size(spikes_w,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ensembles
% %     % Sorting - Full
% %     if dt==0
% %         [~,sorting,~]=get_sorting(spikes_w);
% %     else
% %         [~,sorting,~]=get_sorting_smoothed(spikes_w,dt);
% %     end
% % 
% %     cells_per_ens=floor(N/10);
% %     for e=1:10
% %         ens(e,:)=sorting((e-1)*cells_per_ens+1:e*cells_per_ens);
% %         ense_n(e,:)=e*ones(1,cells_per_ens);
% %     end
% %     Ens(:,1)=reshape(ens,cells_per_ens*10,1);
% %     Ens(:,2)=reshape(ense_n,cells_per_ens*10,1);
% % 
% %     if (N-10*cells_per_ens)>0
% %         Ens(10*cells_per_ens+1:N,1)=sorting(10*cells_per_ens+1:N);
% %         Ens(10*cells_per_ens+1:N,2)=ones(length(sorting(10*cells_per_ens+1:N)),1)*10;
% %     end
% % 
% %     Ens1=Ens;
% %     clear Ens ense ense_n;
% % 
% %     % Sorting - 1st half
% %     if dt==0
% %         [~,sorting,~]=get_sorting(spikes_w);
% %     else
% %         [~,sorting,~]=get_sorting_smoothed(spikes_w(:,1:floor(T/2)),dt);
% %     end
% % 
% %     for e=1:10
% %         ens(e,:)=sorting((e-1)*cells_per_ens+1:e*cells_per_ens);
% %         ense_n(e,:)=e*ones(1,cells_per_ens);
% %     end
% %     Ens(:,1)=reshape(ens,cells_per_ens*10,1);
% %     Ens(:,2)=reshape(ense_n,cells_per_ens*10,1);
% % 
% %     if (N-10*cells_per_ens)>0
% %         Ens(10*cells_per_ens+1:N,1)=sorting(10*cells_per_ens+1:N);
% %         Ens(10*cells_per_ens+1:N,2)=ones(length(sorting(10*cells_per_ens+1:N)),1)*10;
% %     end
% % 
% %     Ens2=Ens;
% %     clear Ens ense ense_n;
% % 
% %     % Sorting - 2nd half
% %     if dt==0
% %         [~,sorting,~]=get_sorting(spikes_w);
% %     else
% %         [~,sorting,~]=get_sorting_smoothed(spikes_w(:,floor(T/2):end),dt);
% %     end
% % 
% %     cells_per_ens=floor(N/10);
% %     for e=1:10
% %         ens(e,:)=sorting((e-1)*cells_per_ens+1:e*cells_per_ens);
% %         ense_n(e,:)=e*ones(1,cells_per_ens);
% %     end
% %     Ens(:,1)=reshape(ens,cells_per_ens*10,1);
% %     Ens(:,2)=reshape(ense_n,cells_per_ens*10,1);
% % 
% %     if (N-10*cells_per_ens)>0
% %         Ens(10*cells_per_ens+1:N,1)=sorting(10*cells_per_ens+1:N);
% %         Ens(10*cells_per_ens+1:N,2)=ones(length(sorting(10*cells_per_ens+1:N)),1)*10;
% %     end
% % 
% %     Ens3=Ens;
% %     clear Ens ense ense_n;

    %Sortings

    [~,sorting_f,~]=get_sorting_smoothed(spikes_w,floor(dt));
    [~,sorting_1,~]=get_sorting_smoothed(spikes_w(:,1:floor(T/2)),floor(dt));
    [~,sorting_2,~]=get_sorting_smoothed(spikes_w(:,floor(T/2):end),floor(dt));

    count=0;
    for i=1:N
        for j=i+1:N
            count=count+1;
            aux_i=find(sorting_f==i);
            aux_j=find(sorting_f==j);
            dist_f(count)=abs(aux_i-aux_j); if (dist_f(count)>floor(N/2)) dist_f(count)=N-dist_f(count); end

            aux_i=find(sorting_1==i);
            aux_j=find(sorting_1==j);
            dist_1(count)=abs(aux_i-aux_j); if (dist_1(count)>floor(N/2)) dist_1(count)=N-dist_1(count); end

            aux_i=find(sorting_2==i);
            aux_j=find(sorting_2==j);
            dist_2(count)=abs(aux_i-aux_j); if (dist_2(count)>floor(N/2)) dist_2(count)=N-dist_2(count); end
        end
    end

    for sh=1:N_sh
        count=0;
        sorting_f_sh=sorting_f(randperm(N));
        sorting_1_sh=sorting_1(randperm(N));
        sorting_2_sh=sorting_2(randperm(N));
        for i=1:N
            for j=i+1:N
                count=count+1;
                aux_i=find(sorting_f_sh==i);
                aux_j=find(sorting_f_sh==j);
                dist_f_sh(count)=abs(aux_i-aux_j); if (dist_f_sh(count)>floor(N/2)) dist_f_sh(count)=N-dist_f_sh(count); end

                aux_i=find(sorting_1_sh==i);
                aux_j=find(sorting_1_sh==j);
                dist_1_sh(count)=abs(aux_i-aux_j); if (dist_1_sh(count)>floor(N/2)) dist_1_sh(count)=N-dist_1_sh(count); end

                aux_i=find(sorting_2_sh==i);
                aux_j=find(sorting_2_sh==j);
                dist_2_sh(count)=abs(aux_i-aux_j); if (dist_2_sh(count)>floor(N/2)) dist_2_sh(count)=N-dist_2_sh(count); end
            end
        end

            corr_shuffle(count_s,sh,1)=corr(dist_f_sh',dist_1_sh');
            corr_shuffle(count_s,sh,2)=corr(dist_f_sh',dist_2_sh');
            corr_shuffle(count_s,sh,3)=corr(dist_1_sh',dist_2_sh');
    end

    corr_data(count_s,1)=corr(dist_f',dist_1');
    corr_data(count_s,2)=corr(dist_f',dist_2');
    corr_data(count_s,3)=corr(dist_1',dist_2');


    clear spikes spikes_d_s spikes_w table_u sorting sorting_1 sorting_2 sorting_1_sh sorting_2_sh sorting_f_sh dist_1 dist_2 dist_f ...
        dist_1_sh dist_2_sh dist_f_sh cells_d sorting_f
 

end

mat=zeros(15,3);
for i=1:15
    for s=1:3
        thr(i,s)=prctile(corr_shuffle(i,:,s),95);
        if(corr_data(i,s)>thr(i,s)) mat(i,s)=1; end
    end    
end


figure
subplot(1,3,1)
histogram(corr_data(:,1),-0.1:0.1:1);
hold on
histogram(thr(:,1),-0.1:0.1:1);
axis([0 1 0 inf])
subplot(1,3,2)
histogram(corr_data(:,2),-0.1:0.1:1);
hold on
histogram(thr(:,2),-0.1:0.1:1);
axis([0 1 0 inf])
subplot(1,3,3)
histogram(corr_data(:,3),-0.1:0.1:1);
hold on
histogram(thr(:,3),-0.1:0.1:1);
axis([0 1 0 inf])

% 
% figure
% histogram

% figure
% histogram2(dist_1,dist_2)
% histogram2(dist_1,dist_2,'DisplayStyle','tile','ShowEmptyBins','on')
