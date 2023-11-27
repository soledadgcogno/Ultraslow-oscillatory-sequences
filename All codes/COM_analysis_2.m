%% Calculates the COM for data, shuffling and TW using a bin size of time_bin_size
clear all 
close all
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

% Params
pixel_size_new=1.18185;
pixel_size_old=1.78211;
num_clus_discr=10;
N_sh=500;
sf=7.73;
count=0;
data=load('C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Final scripts\Outputs\locking_all_sessions_280821\locking_all_sessions.mat');

[big_table, waves] = get_big_table();
time_bin_size_vec=[8,16,40];
pooled_dis_sh=zeros(500*N_sh,length(time_bin_size_vec));
pooled_dis=zeros(500,length(time_bin_size_vec));
pooled_prctile_90=zeros(500,length(time_bin_size_vec));
pooled_prctile_95=zeros(500,length(time_bin_size_vec));
pooled_prctile_5=zeros(500,length(time_bin_size_vec));
pooled_prctile_10=zeros(500,length(time_bin_size_vec));
pooled_dis_sh_1=zeros(500*N_sh,length(time_bin_size_vec));
pooled_dis_1=zeros(500,length(time_bin_size_vec));
pooled_prctile_90_1=zeros(500,length(time_bin_size_vec));
pooled_prctile_95_1=zeros(500,length(time_bin_size_vec));
pooled_prctile_5_1=zeros(500,length(time_bin_size_vec));
pooled_prctile_10_1=zeros(500,length(time_bin_size_vec));
n_seq=zeros(1,length(waves));

disp('Session: '); disp(' ');
for w=1:length(waves)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data
    row_w=data.waves(w);
    disp(w);    
    count=count+1;
    mouse=['L',num2str(data.big_table(data.waves(w),1)),'M',num2str(data.big_table(data.waves(w),2))];
    day=data.big_table(data.waves(w),3);
    s=data.big_table(data.waves(w),4);
    munit=data.big_table(data.waves(w),5);
    mean_p=data.mean_p_all_sessions{w};    
    dt=floor(data.big_table(data.waves(w),8));
    if isinteger(dt)
    else
        dt=floor(dt);
    end
    
%     dt_Ent=floor(big_table(waves(w),11));
    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    file_name_spk=[data.dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_spk,'-mat');
    spikes=full(spikes_d_s); clear spikes_d_s;
    load([data.dpath,strcat('recording_dates_',mouse,'.mat')]);
    file_name_anat=[data.dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_anat,'-mat'); %Anatomical information

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Anatomical position of each cell
    if w>10; pixel_size=pixel_size_old; else; pixel_size=pixel_size_new; end    
    for i=1:size(spikes,1)
        r_i(i,:)=[Anat.med{1,i}(1)*pixel_size,Anat.med{1,i}(2)*pixel_size];
        dist_origin(i)=norm(r_i(i,:));
    end
    % Shuffled anatomical position of each cell
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
    n_seq(w)=size(table_u,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Downsample in time  
    count_t=0;
    for time_bin_size=time_bin_size_vec
        bin_size=time_bin_size/sf;
        count_t=count_t+1;
        for i=1:size(table_u,1) %Loop in sequences
            spikes_seq=spikes(:,table_u(i,1):table_u(i,2));
            [~,T] = size(spikes_seq);
            time_w=floor(T/time_bin_size);      %Total number of time points in the sequence
            win_t=time_bin_size;
            spk_down=zeros(N,time_w);
            %         spk_down_t=zeros(N,time_w);

            %Denoising of cells in each time bin
            for l=1:time_w
                spk_down(:,l)=(sum(spikes_seq(:,(l-1)*win_t +1:(l-1)*win_t +win_t),2));
                %             spk_down_t(:,l)=spk_down(:,l);
                %             spk_down_t(spk_down(:,l)<mean(spk_down(:,l)),l)=0;
            end
            %         seq_m{i}=spk_down;

            %         spk_down=spk_down_t;
            for t=1:size(spk_down,2); com(t,:)=sum(spk_down(:,t).*r_i)/sum(spk_down(:,t)); end %calculates COM
            for t=2:size(com,1); d_com_1s(t-1,:)=norm(com(t,:)-com(1,:)); end %calculates distance of COM from COM at t=1
            for t=2:size(com,1); d_com_1s_c(t-1,:)=norm(com(t,:)-com(t-1,:)); end %calculates distance of COM at t from COM at t-1
%             flowx=diff(com(:,1))/bin_size;  flowy=diff(com(:,2))/bin_size; %calculates flow in x and y directions
%             flow=sqrt(flowx.*flowx + flowy.*flowy); %calculates flow

            for t=1:size(spk_down,2); com_tw(t,:)=sum(spk_down(:,t).*r_i_tw)/sum(spk_down(:,t)); end %calculates COM
            for t=2:size(com,1); d_com_1s_tw(t-1,:)=norm(com_tw(t,:)-com_tw(1,:)); end %calculates distance of COM from COM at t=1
            for t=2:size(com,1); d_com_1s_tw_c(t-1,:)=norm(com_tw(t,:)-com_tw(t-1,:)); end %calculates distance of COM at t from COM at t-1
%             flowx_tw=diff(com_tw(:,1))/bin_size;  flowy_tw=diff(com_tw(:,2))/bin_size; %calculates flow in x and y directions
%             flow_tw=sqrt(flowx_tw.*flowx_tw + flowy_tw.*flowy_tw); %calculates flow

            for sh=1:N_sh
                for t=1:size(spk_down,2); com_sh(t,:,sh)=sum(spk_down(:,t).*r_i_sh(:,:,sh))/sum(spk_down(:,t)); end %calculates COM
                for t=2:size(com,1); d_com_1s_sh(t-1,sh)=norm(com_sh(t,:,sh)-com_sh(1,:,sh)); end %calculates distance of COM from COM at t=1
                for t=2:size(com,1); d_com_1s_sh_c(t-1,sh)=norm(com_sh(t,:,sh)-com_sh(t-1,:,sh)); end %calculates distance of COM from COM at t-1
%                 flowx_sh=diff(com_sh(:,1,sh))/bin_size;  flowy_sh=diff(com_sh(:,2,sh))/bin_size; %calculates flow in x and y directions
%                 flow_sh(:,sh)=sqrt(flowx_sh.*flowx_sh + flowy_sh.*flowy_sh); %calculates flow
            end

       
            %         figure
            %         subplot(1,3,1)
            %         plot(com(:,1),com(:,2));
            %         axis([0 600 0 600])
            %         subplot(1,3,2)
            %         plot(com_tw(:,1),com_tw(:,2));
            %         axis([0 600 0 600])
            %         subplot(1,3,3)
            %         plot(com_sh(:,1,1),com_sh(:,2,1));
            %         axis([0 600 0 600])

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quantification of displacemente of the COM from time point to time point
            %         %Data
            %         score=zeros(1,size(com,1)-1);
            %         for t=2:size(com,1); if (d_com_1s_c(t-1)<prctile(d_com_1s_sh_c(t-1,:),1)); ...
            %                     score(t-1)=1; end; end
            %         max_score(i)=max(score);
            %         frac_score(i)=sum(score)/length(score);
            %
            %         %TW
            %         score_tw=zeros(1,size(com,1)-1);
            %         for t=2:size(com,1); if (d_com_1s_tw_c(t-1)<prctile(d_com_1s_sh_c(t-1,:),1)); ...
            %                     score_tw(t-1)=1; end; end
            %         max_score_tw(i)=max(score_tw);
            %         frac_score_tw(i)=sum(score_tw)/length(score_tw);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flow quantification
            %Data
%             score_flow=zeros(1,size(flow,1)-1);
%             for t=1:size(flow,1); if (flow(t)<prctile(flow_sh(t,:),1)); ...
%                         score_flow(t)=1; end; end
%             max_score_flow(i)=max(score_flow);
%             frac_score_flow(i)=sum(score_flow)/length(score_flow);
% 


            if w==8 && i==19
                figure
                [h,bins]=histcounts(sum(d_com_1s_sh_c,1),375:50:1025);
                bar(bins(1:end-1)+25,h)
                hold on
                xline(sum(d_com_1s_c),'r-','linewidth',5);
                xline(prctile(sum(d_com_1s_sh_c,1),95),'g-','linewidth',5);
                xline(prctile(sum(d_com_1s_sh_c,1),5),'b-','linewidth',5);
                ylabel('Counts');
                xlabel('Cumulative distance travelled (\mum)')
                set(gca,'fontsize',16);
                box off
                axis([400 1000 0 150])
                xticks([400 500 600 700 800 900 1000])
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quantification of total displacemente of the COM
            tot_dis(i,count_t,w)=sum(d_com_1s_c); % number_of_bins x 1 %This is the one i want to plot
            tot_dis_sh(i,:,count_t)=sum(d_com_1s_sh_c,1); % number_of_bins x Nsh
            prctile_90_sh(i,count_t,w)=prctile(tot_dis_sh(i,:,count_t),90);
            prctile_10_sh(i,count_t,w)=prctile(tot_dis_sh(i,:,count_t),10);     
            prctile_99_sh(i,count_t,w)=prctile(tot_dis_sh(i,:,count_t),99);
            prctile_95_sh(i,count_t,w)=prctile(tot_dis_sh(i,:,count_t),95);%This is the one i want to quantify - after Edvard's suggestion
            prctile_5_sh(i,count_t,w)=prctile(tot_dis_sh(i,:,count_t),5); %This is the one i want to plot
            prctile_50_sh(i,count_t,w)=prctile(tot_dis_sh(i,:,count_t),50); %This is the one i want to plot - after Edvard's suggestion
            prctile_1_sh(i,count_t,w)=prctile(tot_dis_sh(i,:,count_t),1);   
            tot_dis_tw(i,count_t,w)=sum(d_com_1s_tw_c);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quantification of displacemente of the COM with respect of t=1
            tot_dis_1(i,count_t,w)=d_com_1s(end);
            tot_dis_sh_1(i,:,count_t)=d_com_1s_sh(end,:);
            prctile_90_sh_1(i,count_t,w)=prctile(tot_dis_sh_1(i,:,count_t),90);
            prctile_95_sh_1(i,count_t,w)=prctile(tot_dis_sh_1(i,:,count_t),95);
            prctile_5_sh_1(i,count_t,w)=prctile(tot_dis_sh_1(i,:,count_t),5);
            prctile_10_sh_1(i,count_t,w)=prctile(tot_dis_sh_1(i,:,count_t),10);    
            prctile_99_sh_1(i,count_t,w)=prctile(tot_dis_sh_1(i,:,count_t),99);
            prctile_1_sh_1(i,count_t,w)=prctile(tot_dis_sh_1(i,:,count_t),1);    
            tot_dis_tw_1(i,count_t,w)=sum(d_com_1s_tw);

            clear spikes_seq cc com spk_down time_w score d_com_1s d_com_1s_c com_sh d_com_1s_sh com_tw spk_down_t d_com_1s_tw flowx_sh flowy_sh flow_sh flow flow_tw flowx flowy ...
                flowx_tw flowy_tw d_com_1s_c d_com_1s_tw_c d_com_1s_sh_c score_flow score_tw_flow aux 
        end
         clear aux; aux=tot_dis_sh(:,:,count_t);
         tot_dis_sh_m(1:length(aux(:)),count_t)=aux(:);

         clear aux1; aux1=tot_dis_sh_1(:,:,count_t);
         tot_dis_sh_m1(1:length(aux1(:)),count_t)=aux1(:);
    end

    pooled_dis=[pooled_dis;tot_dis(:,:,w)]; 
    pooled_dis_1=[pooled_dis_1;tot_dis_1(:,:,w)]; 


    pooled_dis_sh=[pooled_dis_sh;tot_dis_sh_m];     
    pooled_prctile_90=[pooled_prctile_90;prctile_90_sh(:,:,w)];
    pooled_prctile_10=[pooled_prctile_10;prctile_10_sh(:,:,w)];
    pooled_prctile_95=[pooled_prctile_95;prctile_95_sh(:,:,w)];
    pooled_prctile_5=[pooled_prctile_95;prctile_5_sh(:,:,w)];


    pooled_dis_sh_1=[pooled_dis_sh_1;tot_dis_sh_m1];     
    pooled_prctile_90_1=[pooled_prctile_90_1;prctile_90_sh_1(:,:,w)];
    pooled_prctile_10_1=[pooled_prctile_10_1;prctile_10_sh_1(:,:,w)];

    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2 r_i N T spikes mean_p frac_score_tw_flow frac_score_flow ...
        max_score_flow max_score_tw_flow mean_p pp prctile_sh r_i r_i_sh r_i_tw sh_ord sorting_descend spikes cells_d tot_dis_sh_m aux dist_origin ...
        idx_p idx_dist
end

%% New figures

exp_data=[];
shuf_data=[];

for w=1:15
    exp_data=[exp_data;tot_dis(:,:,w)];
    shuf_data=[shuf_data;prctile_5_sh(:,:,w)];    
end

to_delete=find(sum(exp_data,2)==0);
exp_data(to_delete,:)=[];
shuf_data(to_delete,:)=[];

column=3;
ed=[10:50:1200];

[dist_data,edges]=histcounts(exp_data(:,column),ed,'Normalization','probability');
binwidth=edges(2)-edges(1);
[dist_sh,edges]=histcounts(shuf_data(:,column),ed,'Normalization','probability');


fig=figure;
plot(edges(1:end-1)+binwidth/2,cumsum(dist_data),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(dist_sh),'linewidth',2.5);
box off; axis square;
legend({'Data','Shuffle'});
set(gca,'fontsize',16);
ylabel('Cumulative normalized frequency');
xlabel('Cumulative distance travelled (\mum)');

for i=1:3
    aux(i)=length(find(exp_data(:,i)<shuf_data(:,i)));
end

% Using the median, as Edvard suggested
exp_data=[];
shuf_data=[];
for w=1:15
    exp_data=[exp_data;tot_dis(:,:,w)];
    shuf_data=[shuf_data;prctile_50_sh(:,:,w)];  
end

[dist_data,edges]=histcounts(exp_data(:,column),ed,'Normalization','probability');
binwidth=edges(2)-edges(1);
[dist_sh,edges]=histcounts(shuf_data(:,column),ed,'Normalization','probability');

fig=figure;
plot(edges(1:end-1)+binwidth/2,cumsum(dist_sh),'linewidth',2.5);
hold on
plot(edges(1:end-1)+binwidth/2,cumsum(dist_data),'linewidth',2.5);
box off; axis square;
legend({'Shuffle','Data'});
set(gca,'fontsize',16);
ylabel('Cumulative normalized frequency');
xlabel('Cumulative distance travelled (\mum)');

%Quantification for 95th percentile
exp_data=[];
shuf_data=[];
for w=1:15
    exp_data=[exp_data;tot_dis(:,:,w)];
    shuf_data=[shuf_data;prctile_95_sh(:,:,w)];    
end

to_delete=find(sum(exp_data,2)==0);
exp_data(to_delete,:)=[];
shuf_data(to_delete,:)=[];

for i=1:3
    aux(i)=length(find(exp_data(:,i)>shuf_data(:,i)));
end

%% Old figures that we didn't end up using
% %% Figures
% 
% %Histograms for total displacemente of the COM
% pooled_dis(sum(pooled_dis,2)==0,:)=[];
% pooled_dis_sh(sum(pooled_dis_sh,2)==0,:)=[];
% for i=1:count_t %Loop on number of time bin sizes
%     data_sh=pooled_dis_sh(:,i);
%     data=pooled_dis(:,i);
%     figure
%     yyaxis right
%     histogram(data(:),[0:max(data_sh)/100:max(data_sh)]);
%     ylabel('Counts - Data');
%     yyaxis left
%     histogram(data_sh(:),[0:max(data_sh)/100:max(data_sh)]);
%     xlabel('COM cumulative distance travelled (\mum)');
%     ylabel('Counts - Shuffle');
%     set(gca,'fontsize',20,'XColor','k');
% %     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\travelled_distance_histogram_',num2str(i),'.png'));
% %     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\travelled_distance_histogram_',num2str(i),'.fig'));
% %     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\travelled_distance_histogram_',num2str(i),'.svg'));
% 
%     [h_tot(i),p_tot(i),ks2stat_tot(i)] = kstest2(data_sh,data);
%     [pW_tot(i),hW_tot(i),statsW_tot(i)] = ranksum(data_sh,data);
% 
%     clear data_sh data
% end
% % 
% % %Histograms for COM(end)-COM(1)
% % pooled_dis_sh_1(sum(pooled_dis_sh_1,2)==0,:)=[];
% % pooled_dis_1(sum(pooled_dis_1,2)==0,:)=[];
% % clear data_sh data
% % for i=1:count_t
% %     data_sh=pooled_dis_sh_1(:,i);
% %     data=pooled_dis_1(:,i);
% %     figure
% %     yyaxis right
% %     histogram(data(:),[0:max(data_sh)/100:max(data_sh)]);
% %     ylabel('Counts - Data');
% %     yyaxis left
% %     histogram(data_sh(:),[0:max(data_sh)/100:max(data_sh)]);
% %     xlabel('COM net displacement (\mum)');
% %     ylabel('Counts - Shuffle');
% %     set(gca,'fontsize',20,'XColor','k');
% % %     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\net_displacement_histogram_',num2str(i),'.png'));
% % %     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\net_displacement_histogram_',num2str(i),'.fig'));
% % %     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\net_displacement_histogram_',num2str(i),'.svg'));
% %     [h_net(i),p_net(i),ks2stat_net(i)] = kstest2(data_sh,data);
% %     [pW_net(i),hW_net(i),statsW_net(i)] = ranksum(data_sh,data);
% % 
% %     clear data_sh data
% % end
% 
% %Figure of total displacemente of the COM for example session
% cc=lines(count_t);
% clear aux; aux=tot_dis(:,:,8); aux(sum(aux,2)==0,:)=[];
% figure
% for i=count_t%1:count_t
%     plot(aux(:,i),'color',cc(1,:),'linewidth',2)
%     hold on
%     plot(prctile_90_sh(1:size(aux,1),i,8),'--','color',cc(1,:),'linewidth',2);
%     plot(prctile_10_sh(1:size(aux,1),i,8),'-.','color',cc(1,:),'linewidth',2);
% end
% lgd=legend({['Data - ',num2str(round(time_bin_size_vec(1)/sf,1)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile',...
%     ['Data - ',num2str(floor(time_bin_size_vec(2)/sf)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile',...
%     ['Data - ',num2str(floor(time_bin_size_vec(3)/sf)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile'});
% lgd.FontSize=8;
% ylabel({'COM cumulative distance';'travelled (\mum)'});
% xlabel('Cycle #');
% set(gca,'fontsize',20,'YColor','k','XColor','k');
% axis([0 25 500 inf]);
% xticks([1 12 24]);
% box off
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\travelled_distance_example_session.png'));
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\travelled_distance_example_session.fig'));
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\travelled_distance_example_session.svg'));
% 
% %Figure of COM(end)-COM(1) for example session
% % cc=lines(count_t);
% % clear aux; aux=tot_dis_1(:,:,8); aux(sum(aux,2)==0,:)=[];
% % figure
% % for i=1:count_t
% %     plot(aux(:,i),'color',cc(i,:),'linewidth',2)
% %     hold on
% %     plot(prctile_90_sh_1(1:size(aux,1),i,8),'--','color',cc(i,:),'linewidth',2);
% %     plot(prctile_10_sh_1(1:size(aux,1),i,8),'-.','color',cc(i,:),'linewidth',2);
% % end
% % lgd=legend({['Data - ',num2str(round(time_bin_size_vec(1)/sf,1)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile',...
% %     ['Data - ',num2str(floor(time_bin_size_vec(2)/sf)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile',...
% %     ['Data - ',num2str(floor(time_bin_size_vec(3)/sf)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile'});
% % lgd.FontSize=8;
% % ylabel({'COM net displacement (\mum)'});
% % xlabel('Cycle #');
% % set(gca,'fontsize',20,'YColor','k','XColor','k');
% % axis([0 25 0 inf]);
% % xticks([1 12 24]);
% % box off
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\net_displacement_example_session.png'));
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\net_displacement_example_session.fig'));
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\net_displacement_example_session.svg'));
% 
% % 
% % %Figure of total displacemente of the COM for example session - 99th
% % %percentile
% % cc=lines(count_t);
% % clear aux; aux=tot_dis(:,:,8); aux(sum(aux,2)==0,:)=[];
% % figure
% % for i=1:count_t
% %     plot(aux(:,i),'color',cc(i,:),'linewidth',2)
% %     hold on
% %     plot(prctile_90_sh(1:size(aux,1),i,8),'--','color',cc(i,:),'linewidth',2);
% %     plot(prctile_10_sh(1:size(aux,1),i,8),'-.','color',cc(i,:),'linewidth',2);
% % end
% % lgd=legend({['Data - ',num2str(round(time_bin_size_vec(1)/sf,1)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile',...
% %     ['Data - ',num2str(floor(time_bin_size_vec(2)/sf)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile',...
% %     ['Data - ',num2str(floor(time_bin_size_vec(3)/sf)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile'});
% % lgd.FontSize=8;
% % ylabel({'COM total travelled distance (\mum)'});
% % xlabel('Cycle #');
% % set(gca,'fontsize',20,'YColor','k','XColor','k');
% % axis([0 25 0 inf]);
% % xticks([1 12 24]);
% % box off
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\99thprctile-travelled_distance_example_session.png'));
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\99thprctile-travelled_distance_example_session.fig'));
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\99thprctile-travelled_distance_example_session.svg'));
% % 
% % %Figure of COM(end)-COM(1) for example session - 99th
% % %percentile
% % cc=lines(count_t);
% % clear aux; aux=tot_dis_1(:,:,8); aux(sum(aux,2)==0,:)=[];
% % figure
% % for i=1:count_t
% %     plot(aux(:,i),'color',cc(i,:),'linewidth',2)
% %     hold on
% %     plot(prctile_90_sh_1(1:size(aux,1),i,8),'--','color',cc(i,:),'linewidth',2);
% %     plot(prctile_10_sh_1(1:size(aux,1),i,8),'-.','color',cc(i,:),'linewidth',2);
% % end
% % lgd=legend({['Data - ',num2str(round(time_bin_size_vec(1)/sf,1)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile',...
% %     ['Data - ',num2str(floor(time_bin_size_vec(2)/sf)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile',...
% %     ['Data - ',num2str(floor(time_bin_size_vec(3)/sf)),' s'],'Shuffle 90^t^h percentile','Shuffle 10^t^h percentile'});
% % lgd.FontSize=8;
% % ylabel({'COM net displacement (\mum)'});
% % xlabel('Cycle #');
% % set(gca,'fontsize',20,'YColor','k','XColor','k');
% % axis([0 25 0 inf]);
% % xticks([1 12 24]);
% % box off
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\99thprctile-net_displacement_example_session.png'));
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\99thprctile-net_displacement_example_session.fig'));
% % saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\Traveling Waves\COM\99thprctile-net_displacement_example_session.svg'));
% 
% %% Quantifications of the figures from the previous cell
% 
% %Quantifications with 10 and 90 prctile
% count_seq=0;
% for w=1:length(waves)
%     clear aux; aux=tot_dis(:,:,w); aux(sum(aux,2)==0,:)=[];
%     for s=1:n_seq(w)
%         count_seq=count_seq+1;
%         for i=1:count_t
%             if(aux(s,i)>prctile_10_sh(s,i,w) && aux(s,i)<prctile_90_sh(s,i,w) ) sig_tot(count_seq,i)= 0;
% 
%             else sig_tot(count_seq,i)= 1;
%             end
%         end
%     end
% end
% 
% for i=1:3
% fraction_tot(i)=sum(sig_tot(:,i))/size(sig_tot,1);
% pout_tot(i)=myBinomTest(sum(sig_tot(:,i)),size(sig_tot,1),0.05,'one');
% end
% 
% 
% count_seq=0;
% for w=1:length(waves)
%     clear aux; aux=tot_dis_1(:,:,w); aux(sum(aux,2)==0,:)=[];
%     for s=1:n_seq(w)
%         count_seq=count_seq+1;
%         for i=1:count_t
%             if(aux(s,i)>prctile_10_sh_1(s,i,w) && aux(s,i)<prctile_90_sh_1(s,i,w) ) sig_net(count_seq,i)= 0;
%             else sig_net(count_seq,i)= 1;
%             end
%         end
%     end
% end
% 
% 
% %Quantifications with 5 prctile
% count_seq=0;
% for w=1:length(waves)
%     clear aux; aux=tot_dis(:,:,w); aux(sum(aux,2)==0,:)=[];
%     for s=1:n_seq(w)
%         count_seq=count_seq+1;
%         for i=1:count_t
%             if(aux(s,i)<prctile_5_sh(s,i,w) ) sig_tot(count_seq,i)= 1;
% 
%             else sig_tot(count_seq,i)= 0;
%             end
%         end
%     end
% end
% 
% for i=1:3
%     outcomes(i)=sum(sig_tot(:,i));
%     fraction_tot(i)=sum(sig_tot(:,i))/size(sig_tot,1);
% end
% 
% 
% 
% %Quantifications with 95 prctile
% count_seq=0;
% for w=1:length(waves)
%     clear aux; aux=tot_dis(:,:,w); aux(sum(aux,2)==0,:)=[];
%     for s=1:n_seq(w)
%         count_seq=count_seq+1;
%         for i=1:count_t
%             if(aux(s,i)<prctile_95_sh(s,i,w) ) sig_tot(count_seq,i)= 0;
% 
%             else sig_tot(count_seq,i)= 1;
%             end
%         end
%     end
% end
% 
% for i=1:3
% fraction_tot(i)=sum(sig_tot(:,i))/size(sig_tot,1);
% pout_tot(i)=myBinomTest(sum(sig_tot(:,i)),size(sig_tot,1),0.05,'one');
% end
% 
% 
% %Quantifications with 1 and 99 prctile
% sig_tot=zeros(421,3);
% count_seq=0;
% for w=1:length(waves)
%     clear aux; aux=tot_dis(:,:,w); aux(sum(aux,2)==0,:)=[];
%     for s=1:n_seq(w)
%         count_seq=count_seq+1;
%         for i=1:count_t
%             if(aux(s,i)>prctile_1_sh(s,i,w) && aux(s,i)<prctile_99_sh(s,i,w) ) sig_tot(count_seq,i)= 0;
%             else sig_tot(count_seq,i)= 1;
%             end
%         end
%     end
% end
% 
% sig_net=zeros(421,3);
% count_seq=0;
% for w=1:length(waves)
%     clear aux; aux=tot_dis_1(:,:,w); aux(sum(aux,2)==0,:)=[];
%     for s=1:n_seq(w)
%         count_seq=count_seq+1;
%         for i=1:count_t
%             if(aux(s,i)>prctile_1_sh_1(s,i,w) && aux(s,i)<prctile_99_sh_1(s,i,w) ) sig_net(count_seq,i)= 0;
%             else sig_net(count_seq,i)= 1;
%             end
%         end
%     end
% end