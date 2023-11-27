clear all
close all

[big_table, waves] = get_big_table();

%% Calculates the COM for different temporal resolutions and without shuffling

N_sh=100;
count=0;

for w=1:length(waves)
    row_w=waves(w);
    disp(w);
    
    count=count+1;
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);
    munit=big_table(waves(w),5);
    
    dt=floor(big_table(waves(w),8));
    if isinteger(dt)
    else
        dt=floor(dt);
    end
    
%     dt_Ent=floor(big_table(waves(w),11));

    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_spk,'-mat');
    spikes=full(spikes_d_s);
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_anat,'-mat'); %Anatomical information


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distance in the tissue

    pixel_size_new=1.18185;
    pixel_size_old=1.78211;

    if w>10
        pixel_size=pixel_size_old;
    else
        pixel_size=pixel_size_new;
    end
    
    %anatomical position of each cell
    for i=1:size(spikes,1)
        r_i(i,:)=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
%         dist_origin(i)=norm(r_i);
    end

    %spikes from all sequences
    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);        
    
     %Original timescale  
    figure
    spikes_w=[];
    for i=1:size(table_u,1)
        spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
        spikes_seq=spikes(:,table_u(i,1):table_u(i,2));
        
        for t=1:size(spikes_seq,2); com(t,:)=sum(spikes_seq(:,t).*r_i)/sum(spikes_seq(:,t)); end
        for t=2:size(com,1); d_com_129ms(t-1,:)=norm(com(t,:)-com(1,:)); end %calculates distance of COM

        cc=copper(size(spikes_seq,2));
        subplot(ceil(sqrt(size(table_u,1))),ceil(sqrt(size(table_u,1))),i); hold on;...
            for t=1:size(spikes_seq,2); plot(com(t,1),com(t,2),'o','Color',cc(t,:),'MarkerFaceColor',cc(t,:)); end; ... 
            axis([0 700 0 700]); title(['Sequence = ',num2str(i)]); 

        clear spikes_seq cc com
    end

     %Downsamples in time with 1 s  
    figure
    spikes_w=[];
    for i=1:size(table_u,1)
        spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
        spikes_seq=spikes(:,table_u(i,1):table_u(i,2));
        [~,T] = size(spikes_seq);
        time_w=floor(T/8);      %Number of time points
         win_t=8;       
        spk_down=zeros(N,time_w);
        for l=1:time_w
            spk_down(:,l)=(sum(spikes_seq(:,(l-1)*win_t +1:(l-1)*win_t +win_t),2));
        end

        for t=1:size(spk_down,2); com(t,:)=sum(spk_down(:,t).*r_i)/sum(spk_down(:,t)); end %calculates COM
        for t=2:size(com,1); d_com_1s(t-1,:)=norm(com(t,:)-com(1,:)); end %calculates distance of COM

        cc=copper(size(spk_down,2));
        subplot(ceil(sqrt(size(table_u,1))),ceil(sqrt(size(table_u,1))),i); hold on; for t=1:size(spk_down,2); plot(com(t,1),com(t,2),...
                'o','Color',cc(t,:),'MarkerFaceColor',cc(t,:)); end; ...
            axis([0 700 0 700]); title(['Sequence = ',num2str(i)]); set(gca,'fontsize',18);
 
        clear spikes_seq cc com spk_down time_w
    end
    
        
    %Downsamples in time with characteristic time scale  
    figure
    spikes_w=[];
    for i=1:size(table_u,1)
        spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
        spikes_seq=spikes(:,table_u(i,1):table_u(i,2));
        [~,T] = size(spikes_seq);
        time_w=floor(T/dt);      %Number of time points
         win_t=dt;       
        spk_down=zeros(N,time_w);
        for l=1:time_w
            spk_down(:,l)=(sum(spikes_seq(:,(l-1)*win_t +1:(l-1)*win_t +win_t),2));
        end

        for t=1:size(spk_down,2); com(t,:)=sum(spk_down(:,t).*r_i)/sum(spk_down(:,t)); end
        for t=2:size(com,1); d_com_dt(t-1,:)=norm(com(t,:)-com(1,:)); end %calculates distance of COM

        cc=copper(size(spk_down,2));
        subplot(ceil(sqrt(size(table_u,1))),ceil(sqrt(size(table_u,1))),i); hold on; for t=1:size(spk_down,2); plot(com(t,1),com(t,2),...
                'o','Color',cc(t,:),'MarkerFaceColor',cc(t,:)); end; ...
            axis([0 700 0 700]); title(['Sequence = ',num2str(i)]); set(gca,'fontsize',18);
 
        clear spikes_seq cc com spk_down time_w
    end

    %Shuffle anatomical position of each cell
    for sh=1:N_sh
        sh_ord=randperm(N);
        r_i_sh=r_i(sh_ord,:);

        for i=1%:size(table_u,1)
            spikes_seq=spikes(:,table_u(i,1):table_u(i,2));
            [~,T] = size(spikes_seq);
            time_w=floor(T/8);      %Number of time points
            win_t=8;
            spk_down=zeros(N,time_w);
            for l=1:time_w
                spk_down(:,l)=(sum(spikes_seq(:,(l-1)*win_t +1:(l-1)*win_t +win_t),2));
            end

            for t=1:size(spk_down,2); com_sh(t,:)=sum(spk_down(:,t).*r_i_sh)/sum(spk_down(:,t)); end %calculates COM
            for t=2:size(com_sh,1); d_com_1s_sh(t-1,sh)=norm(com_sh(t,:)-com_sh(t-1,:)); end %calculates distance of COM

%             cc=copper(size(spk_down,2));
%             subplot(ceil(sqrt(size(table_u,1))),ceil(sqrt(size(table_u,1))),i); hold on; for t=1:size(spk_down,2); plot(com(t,1),com(t,2),...
%                     'o','Color',cc(t,:),'MarkerFaceColor',cc(t,:)); end; ...
%                 axis([0 700 0 700]); title(['Sequence = ',num2str(i)]); set(gca,'fontsize',18);

            clear spikes_seq cc com spk_down time_w
        end

    end
    

   
    
    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2 r_i table_u N T spikes
    
end

%% Calculates the COM for only 1 second and also for shuffling 

N_sh=100;
count=0;
load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\locking_all_sessions.mat');


for w=8%1:length(waves)
    row_w=waves(w);
    disp(w);
    
    count=count+1;
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);
    munit=big_table(waves(w),5);
    mean_p=mean_p_all_sessions{w};

    
    dt=floor(big_table(waves(w),8));
    if isinteger(dt)
    else
        dt=floor(dt);
    end
    
%     dt_Ent=floor(big_table(waves(w),11));

    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_spk,'-mat');
    spikes=full(spikes_d_s);
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_anat,'-mat'); %Anatomical information


    for i=1:size(spikes,1)
        spikes(i,:)=(full(fire_rate(spikes(i,:),8,'g')));  % 9 seconds for L8M2 Day 19 ; 15 seconds for L9M4 Day 17
    end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Distance in the tissue

    pixel_size_new=1.18185;
    pixel_size_old=1.78211;

    if w>10
        pixel_size=pixel_size_old;
    else
        pixel_size=pixel_size_new;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Anatomical position of each cell
    for i=1:size(spikes,1)
        r_i(i,:)=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
        dist_origin(i)=norm(r_i(i,:));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Shuffled anatomical position of each cell
    N=size(spikes,1);
    for sh=1:N_sh
        sh_ord=randperm(N);
        r_i_sh(:,:,sh)=r_i(sh_ord,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation of traveling wave
        
     [~,idx_p]=sort(mean_p,'ascend');
     [~,idx_dist]=sort(dist_origin,'ascend');

     [~,sorting_descend,~] = get_sorting(spikes);

     for i=1:N
         r_i_tw(sorting_descend(i),:)=r_i(idx_dist(i),:);
     end

     pp=discretize(mean_p,[-pi:2*pi/30:pi]);
     cc=jet(30);

%      figure
%      subplot(1,2,1)
%      hold on
%      for n=1:N
%          scatter(r_i(n,1),r_i(n,2),[],cc(pp(n),:),"filled");
%      end
%      subplot(1,2,2)
%      hold on
%      for n=1:N
%          scatter(r_i_tw(n,1),r_i_tw(n,2),[],cc(pp(n),:),"filled");
%      end

   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spikes from all sequences
    [table_u,N,~]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);        
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Downsamples in time with 1 s  


    time_bin_size=16;
    for i=19%1:size(table_u,1)

        spikes_seq=spikes(:,table_u(i,1):table_u(i,2));
        [~,T] = size(spikes_seq);
        time_w=floor(T/time_bin_size);      %Total number of time points in the sequence
        win_t=time_bin_size;
        spk_down=zeros(N,time_w);
        spk_down_t=zeros(N,time_w);

        for l=1:time_w
            spk_down(:,l)=(sum(spikes_seq(:,(l-1)*win_t +1:(l-1)*win_t +win_t),2));
            spk_down_t(:,l)=spk_down(:,l);
            spk_down_t(spk_down(:,l)<mean(spk_down(:,l)),l)=0;
        end

        seq_m{i}=spk_down;

        spk_down=spk_down_t;
        for t=1:size(spk_down,2); com(t,:)=sum(spk_down(:,t).*r_i)/sum(spk_down(:,t)); end %calculates COM
        for t=2:size(com,1); d_com_1s(t-1,:)=norm(com(t,:)-com(1,:)); end %calculates distance of COM

        for t=1:size(spk_down,2); com_tw(t,:)=sum(spk_down(:,t).*r_i_tw)/sum(spk_down(:,t)); end %calculates COM
        for t=2:size(com,1); d_com_1s_tw(t-1,:)=norm(com_tw(t,:)-com_tw(1,:)); end %calculates distance of COM

        for sh=1:N_sh
            for t=1:size(spk_down,2); com_sh(t,:,sh)=sum(spk_down(:,t).*r_i_sh(:,:,sh))/sum(spk_down(:,t)); end %calculates COM
            for t=2:size(com,1); d_com_1s_sh(t-1,sh)=norm(com_sh(t,:,sh)-com_sh(t-1,:,sh)); end %calculates distance of COM
        end


        figure
        subplot(1,3,1)
        plot(com(:,1),com(:,2));
        axis([0 600 0 600])
        subplot(1,3,2)
        plot(com_tw(:,1),com_tw(:,2));
        axis([0 600 0 600])
        subplot(1,3,3)
        plot(com_sh(:,1,1),com_sh(:,2,1));
        axis([0 600 0 600])


        score=zeros(1,size(com,1)-1);
        for t=2:size(com,1); if (d_com_1s(t-1)<prctile(d_com_1s_sh(t-1,:),1)); ...
                    score(t-1)=1; end; end

        max_score(i)=max(score);
        clear spikes_seq cc com spk_down time_w score d_com_1s com_sh d_com_1s_sh com_tw
    end


    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2 r_i table_u N T spikes mean_p

end
