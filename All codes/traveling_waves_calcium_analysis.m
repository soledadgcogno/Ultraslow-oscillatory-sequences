% Code for snippets of traveling wave as a function of time

clear all
close all

[big_table, waves] = get_big_table();


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


%     for i=1:size(spikes,1)
%         spikes(i,:)=(full(fire_rate(spikes(i,:),8,'g')));  % 9 seconds for L8M2 Day 19 ; 15 seconds for L9M4 Day 17
%     end

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

    [sorting_ascend,sorting_descend,sorting_0] = get_sorting(spikes);

    for i=1:N
%                 r_i_tw(idx_p(i),:)=r_i(idx_dist(i),:);
                r_i_tw(sorting_descend(i),:)=r_i(idx_dist(i),:);
    end
% % % %    for i=1:size(spikes,1)
% % % %         r_i(i,:)=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
% % % %         dist_origin(i)=norm(r_i(i,:));
% % % %     end

    clear delta
    count=0;
    for i=400%:N
        for j=i+1:N
            count=count+1;
            delta(count,1)=j-i;
            delta(count,2)=norm(r_i_tw(sorting_descend(i),:)-r_i_tw(sorting_descend(j),:));
        end
    end

%     figure
%     scatter(delta(:,1),delta(:,2))

%     clear r_i_tw
%     newlist=r_i;
%     r_i_tw(sortind_descend(1),:)=r_i(idx_dist(1),:); 
%     newlist(idx_dist(1),:)=[];
%     
%     for i=2:N %CHECK AND PLOT DT VS DEPOS
%         aux=pdist2(newlist, r_i_tw(sortind_descend(i-1),:));
%         [~,aux2]=min(aux);
%         r_i_tw(sortind_descend(i),:)=newlist(aux2,:); 
%         newlist(aux2,:)=[];
%     end
% 
%     clear delta2
%     count=0;
%     for i=1:N
%         for j=i+1:N
%             count=count+1;
%             delta2(count,1)=j-i;
%             delta2(count,2)=norm(r_i_tw(sortind_descend(i),:)-r_i_tw(sortind_descend(j),:));
%         end
%     end
%     
%     figure
%     scatter(delta2(:,1),delta2(:,2))

   
    pp=discretize(mean_p,[-pi:2*pi/30:pi]);
    cc=jet(30);

    figure
    subplot(1,2,1)
    hold on
    for n=1:N
        scatter(r_i(n,1),r_i(n,2),[],cc(pp(n),:),"filled");
    end    
    subplot(1,2,2)
    hold on
    for n=1:N
        scatter(r_i_tw(n,1),r_i_tw(n,2),[],cc(pp(n),:),"filled");
    end 
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spikes from all sequences
    [table_u,N,~]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);        
    
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binning of position

    bin_size=12; %in um
    n_bins=600/bin_size;

    posx=discretize(r_i(:,1),0:bin_size:900); 
    posy=discretize(r_i(:,2),0:bin_size:900);
   
    C1=[];
    C2=[];
    for l=1:n_bins
        C1 = [C1;(l*ones(n_bins,1))];
        C2 = [C2;(1:n_bins)'];
    end
    C=[C1,C2];

    for n=1:N
        a=find(posx(n)==C(:,1));
        b=find(posy(n)==C(:,2));
        
        idx(n)=intersect(a,b);

        clear a b
    end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binning of shuffled
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% positions

%     bin_size=18;%18
    n_bins=600/bin_size;

    for sh=1:N_sh
        posx_sh=discretize(r_i_sh(:,1,sh),0:bin_size:900);
        posy_sh=discretize(r_i_sh(:,2,sh),0:bin_size:900);

        for n=1:N
            a=find(posx_sh(n)==C(:,1));
            b=find(posy_sh(n)==C(:,2));
            idx_sh(n,sh)=intersect(a,b);
            clear a b
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binning of positions for
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% traveling wave

%     bin_size=18;%18
    n_bins=600/bin_size;

    posx_tw=discretize(r_i_tw(:,1),0:bin_size:900);
    posy_tw=discretize(r_i_tw(:,2),0:bin_size:900);
 

    for n=1:N
        a=find(posx_tw(n)==C(:,1));
        b=find(posy_tw(n)==C(:,2));        
        idx_tw(n)=intersect(a,b);
        clear a b
    end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Downsample in time with 1 s  

    num_panels=8;
    time_bin_size=16;

    for i=19% 1:size(table_u,1) %Loop on all sequences
        
        spikes_seq=spikes(:,table_u(i,1):table_u(i,2));
        [~,T] = size(spikes_seq);
        time_w=floor(T/time_bin_size);      %Total number of time points in the sequence
        win_t=time_bin_size;

        spk_down=zeros(N,time_w);
        spk_down_t=zeros(N,time_w);
        spk_pos=zeros(n_bins,n_bins,time_w);
        spk_pos_tw=zeros(n_bins,n_bins,time_w);
        spk_pos_sh=zeros(n_bins,n_bins,time_w);

        for l=1:time_w
            spk_down(:,l)=(sum(spikes_seq(:,(l-1)*win_t +1:(l-1)*win_t +win_t),2));
            spk_down_t(:,l)=spk_down(:,l);
            spk_down_t(spk_down(:,l)<mean(spk_down(:,l)),l)=0;


            for m=1:size(C,1)
                if sum(idx==m)>0
                spk_pos(C(m,1),C(m,2),l)=sum(spk_down_t(idx==m,l));
                end

                if sum(idx_tw==m)>0
                spk_pos_tw(C(m,1),C(m,2),l)=sum(spk_down_t(idx_tw==m,l));
                end

                if sum(idx_tw==m)>0
                spk_pos_sh(C(m,1),C(m,2),l)=sum(spk_down_t(idx_sh(:,1)==m,l));
                end
            end     
        end

        % Creates the figure
        aux=floor(time_w/num_panels);
        timepoints_f=aux*[1:num_panels];

        timepoints_f(1)=[];
        timepoints_f(end)=[];
        num_panels=length(timepoints_f);


        spk_down_b=spk_down_t;
        spk_down_b(spk_down_b>0)=1;

        figure
        hold on
        for p=1:num_panels
            subplot(3,num_panels,p); 
            hold on
            for n=1:N
                scatter(1:time_w,n*spk_down_b(sorting_ascend(n),:),1,'k.');
            end
            axis([1 time_w 1 N])
            hold on
            for n=1:N 
                scatter(timepoints_f(p),n*spk_down_b(sorting_ascend(n),timepoints_f(p)),20,'r.');
            end
            axis square

            for p=1:num_panels
                subplot(4,num_panels,p+num_panels); imagesc(spk_pos(:,:,timepoints_f(p))); axis square
            end  

            for p=1:num_panels
                subplot(4,num_panels,p+2*num_panels); imagesc(spk_pos_sh(:,:,timepoints_f(p))); axis square
            end   

            for p=1:num_panels
                subplot(4,num_panels,p+3*num_panels); imagesc(spk_pos_tw(:,:,timepoints_f(p))); axis square
            end

        end

        
       

        seq_m{i}=spk_down;

      
        clear spikes_seq cc com spk_down time_w score d_com_1s com_sh d_com_1s_sh com_tw
    end


   
    
    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2 r_i table_u N T spikes mean_p
    
end

%% Code for correlations inside the space domain 

clear all
close all

[big_table, waves] = get_big_table();


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

    for i=1:N
        r_i_tw(idx_p(i),:)=r_i(idx_dist(i),:);
    end

   
    pp=discretize(mean_p,[-pi:2*pi/30:pi]);
    cc=jet(30);

    figure
    subplot(1,2,1)
    hold on
    for n=1:N
        scatter(r_i(n,1),r_i(n,2),[],cc(pp(n),:),"filled");
    end    
    subplot(1,2,2)
    hold on
    for n=1:N
        scatter(r_i_tw(n,1),r_i_tw(n,2),[],cc(pp(n),:),"filled");
    end 
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spikes from all sequences
    [table_u,N,~]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);        
    
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binning of positions

    bin_size=100;%18
    n_bins=900/bin_size;

    posx=discretize(r_i(:,1),0:bin_size:900);
    posy=discretize(r_i(:,2),0:bin_size:900);
   
    C1=[];
    C2=[];
    for l=1:n_bins
        C1 = [C1;(l*ones(n_bins,1))];
        C2 = [C2;(1:n_bins)'];
    end
    C=[C1,C2];

    for n=1:N
        a=find(posx(n)==C(:,1));
        b=find(posy(n)==C(:,2));
        
        idx(n)=intersect(a,b);

        clear a b
    end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binning of shuffled
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% positions

    bin_size=100;%18
    n_bins=900/bin_size;

    for sh=1:N_sh
        posx_sh=discretize(r_i_sh(:,1,sh),0:bin_size:900);
        posy_sh=discretize(r_i_sh(:,2,sh),0:bin_size:900);

        for n=1:N
            a=find(posx_sh(n)==C(:,1));
            b=find(posy_sh(n)==C(:,2));
            idx_sh(n,sh)=intersect(a,b);
            clear a b
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binning of positions for
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% traveling wave

    bin_size=100;%18
    n_bins=900/bin_size;

    posx_tw=discretize(r_i_tw(:,1),0:bin_size:900);
    posy_tw=discretize(r_i_tw(:,2),0:bin_size:900);
 

    for n=1:N
        a=find(posx_tw(n)==C(:,1));
        b=find(posy_tw(n)==C(:,2));        
        idx_tw(n)=intersect(a,b);
        clear a b
    end


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Downsamples in time with 1 s  

    for i=1:size(table_u,1)
        
        spikes_seq=spikes(:,table_u(i,1):table_u(i,2));
        [~,T] = size(spikes_seq);
        time_w=floor(T/8);      %Number of time points
        win_t=8;
        spk_down=zeros(N,time_w);
        spk_pos=zeros(n_bins,n_bins,time_w);
        for l=1:time_w
            spk_down(:,l)=(sum(spikes_seq(:,(l-1)*win_t +1:(l-1)*win_t +win_t),2));

%             for m=1:size(C,1)
%                 if sum(idx==m)>0
%                 spk_pos(C(m,1),C(m,2),l)=sum(spk_down(idx==m,l));
%                 end
%             end

        end

        correlations=[];
        for m=1:size(C,1)
            if sum(idx==m)>1
                Pcorr=corr(spk_down(idx==m,:)');
                aux2=triu(Pcorr,1);
                aux3=aux2(:);
                aux3(aux3==0)=[];
                correlations=[correlations;aux3];
            end
        end
        cdf=histcounts(correlations,[-1:0.1:1],'Normalization','cdf');

       correlations_tw=[];
        for m=1:size(C,1)
            if sum(idx_tw==m)>1
                Pcorr_tw=corr(spk_down(idx_tw==m,:)');
                aux2=triu(Pcorr_tw,1);
                aux3=aux2(:);
                aux3(aux3==0)=[];
                correlations_tw=[correlations_tw;aux3];
            end
        end
        cdf_tw=histcounts(correlations_tw,[-1:0.1:1],'Normalization','cdf');

        figure
        hold on
        for sh=1:N_sh
            correlations_sh=[];
            for m=1:size(C,1)
                if sum(idx_sh(:,sh)==m)>1
                    Pcorr_sh=corr(spk_down(idx_sh(:,sh)==m,:)');
                    aux2=triu(Pcorr_sh,1);
                    aux3=aux2(:);
                    aux3(aux3==0)=[];%CORREGIR
                    correlations_sh=[correlations_sh;aux3];
                end
            end
            correlations_sh_full{sh}=correlations_sh;
            cdf_sh(sh,:)=histcounts(correlations_sh,[-1:0.1:1],'Normalization','cdf');

            plot(cdf_sh(sh,:),'color',[0.9290 0.6940 0.1250])
        end
        plot(cdf,'k','linewidth',3); 
        plot(cdf_tw,'b','linewidth',3);

        
%         % Creates the figure
%         num_panels=6;
%         aux=floor(time_w/num_panels);
%         timepoints_f=aux*[1:num_panels];
% 
%         figure
%         for p=1:num_panels
%             subplot(1,num_panels,p); imagesc(spk_pos(:,:,timepoints_f(p))); axis square
%         end

        seq_m{i}=spk_down;

      
        clear spikes_seq cc com spk_down time_w score d_com_1s com_sh d_com_1s_sh com_tw
    end


   
    
    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2 r_i table_u N T spikes mean_p
    
end

%% Code for correlations of "spatial ensembles" 

clear all
close all

[big_table, waves] = get_big_table();


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

    for i=1:N
        r_i_tw(idx_p(i),:)=r_i(idx_dist(i),:);
    end

   
    pp=discretize(mean_p,[-pi:2*pi/30:pi]);
    cc=jet(30);

    figure
    subplot(1,2,1)
    hold on
    for n=1:N
        scatter(r_i(n,1),r_i(n,2),[],cc(pp(n),:),"filled");
    end    
    subplot(1,2,2)
    hold on
    for n=1:N
        scatter(r_i_tw(n,1),r_i_tw(n,2),[],cc(pp(n),:),"filled");
    end 
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spikes from all sequences
    [table_u,N,~]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);        
    
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binning of positions

    bin_size=100;%18
    n_bins=900/bin_size;

    posx=discretize(r_i(:,1),0:bin_size:900);
    posy=discretize(r_i(:,2),0:bin_size:900);
   
    C1=[];
    C2=[];
    for l=1:n_bins
        C1 = [C1;(l*ones(n_bins,1))];
        C2 = [C2;(1:n_bins)'];
    end
    C=[C1,C2];

    for n=1:N
        a=find(posx(n)==C(:,1));
        b=find(posy(n)==C(:,2));
        
        idx(n)=intersect(a,b);

        clear a b
    end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binning of shuffled
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% positions

    bin_size=100;%18
    n_bins=900/bin_size;

    for sh=1:N_sh
        posx_sh=discretize(r_i_sh(:,1,sh),0:bin_size:900);
        posy_sh=discretize(r_i_sh(:,2,sh),0:bin_size:900);

        for n=1:N
            a=find(posx_sh(n)==C(:,1));
            b=find(posy_sh(n)==C(:,2));
            idx_sh(n,sh)=intersect(a,b);
            clear a b
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binning of positions for
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% traveling wave

    bin_size=100;%18
    n_bins=900/bin_size;

    posx_tw=discretize(r_i_tw(:,1),0:bin_size:900);
    posy_tw=discretize(r_i_tw(:,2),0:bin_size:900);
 

    for n=1:N
        a=find(posx_tw(n)==C(:,1));
        b=find(posy_tw(n)==C(:,2));        
        idx_tw(n)=intersect(a,b);
        clear a b
    end


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Downsample in time with 1 s  

    for i=1:size(table_u,1)
        
        spikes_seq=spikes(:,table_u(i,1):table_u(i,2));
        [~,T] = size(spikes_seq);
        time_w=floor(T/8);      %Number of time points
        win_t=8;
        spk_down=zeros(N,time_w);
        spk_pos=zeros(n_bins,n_bins,time_w);
        for l=1:time_w
            spk_down(:,l)=(sum(spikes_seq(:,(l-1)*win_t +1:(l-1)*win_t +win_t),2));
        end

        sp_ensemble=zeros(size(C,1),time_w);
        sp_ensemble_tw=zeros(size(C,1),time_w);

        for m=1:size(C,1)
            if sum(idx==m)>1; sp_ensemble(m,:)=mean(spk_down(idx==m,:)); end
            if sum(idx_tw==m)>1; sp_ensemble_tw(m,:)=mean(spk_down(idx_tw==m,:)); end
        end
        
        table_dist_lag=nan(n_bins*n_bins,n_bins*n_bins,2);
        table_dist_lag_tw=nan(n_bins*n_bins,n_bins*n_bins,2);

        for m=1:size(C,1)
            for n=m+1:size(C,1)

                if (sum(sp_ensemble(m,:))>0 && sum(sp_ensemble(n,:))>0)
                    [r,lags]=xcorr(sp_ensemble(m,:),sp_ensemble(n,:));
                    [~,lagm]=max(r);
                    table_dist_lag(m,n,1)=lags(lagm)*(win_t)/7.73;
                    table_dist_lag(m,n,2)=pdist([C(m,:);C(n,:)]); %USE RIGHT UNITS
                end

                if (sum(sp_ensemble_tw(m,:))>0 && sum(sp_ensemble_tw(n,:))>0)
                    [r,lags]=xcorr(sp_ensemble_tw(m,:),sp_ensemble_tw(n,:));
                    [~,lagm]=max(r);
                    table_dist_lag_tw(m,n,1)=lags(lagm)*(win_t)/7.73;
                    table_dist_lag_tw(m,n,2)=pdist([C(m,:);C(n,:)]); %USE RIGHT UNITS
                end

            end
        end

        figure
        scatter(table_dist_lag(:,:,1),table_dist_lag(:,:,2));
        
        figure
        scatter(table_dist_lag_tw(:,:,1),table_dist_lag_tw(:,:,2));


        seq_m{i}=spk_down;

      
        clear spikes_seq cc com spk_down time_w score d_com_1s com_sh d_com_1s_sh com_tw
    end


   
    
    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2 r_i table_u N T spikes mean_p
    
end

