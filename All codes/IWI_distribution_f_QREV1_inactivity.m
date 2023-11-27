%Version 1 with fixed window length

clear all
close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

[big_table, waves] = get_big_table();

% IWI distribution
clus=10;
IWI_ses_pooled=[];
sf=7.73; 
win_length=5; %in seconds
flipped_sorting(1,15)=0;
count=0;
for w=1:length(waves)
    row_w=waves(w);
    disp(w);
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);
    munit=big_table(waves(w),5);
    
    dt=floor(big_table(waves(w),8));
    if isinteger(dt)
    else
        dt=floor(dt);
    end
    
    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_spk,'-mat');
    spikes=full(spikes_d_s);
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);

    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);        

    %Separation into ensembles
    [sorting_ascend,sortind_descend,sorting_0] = get_sorting(spikes) ;

    %Checks that the wave has an upwards direction
    clear FRpre
    if dt==0
        FRpre = spikes_downsample(spikes(sortind_descend,:),clus,1);
    else
        FRpre = spikes_downsample(spikes(sortind_descend,:),clus,dt);
    end

    clear signal; [~,signal]=max(FRpre); %Calculate signal by keeping the cluster with largest activity
    clear signal_s; signal_s(1)=signal(1);
    signal_s(2:length(find(diff(signal)~=0))+1) = signal(find(diff(signal)~=0) + 1);
    clear aux; aux=diff(signal_s);
    clear pos neg; pos=length(find(aux>6)); neg=length(find(aux<-6));
    if pos>neg; flipped_sorting(w)=1; sortind_descend=flip(sortind_descend); end

    figure
    spy(spikes(sortind_descend,:));
    pbaspect([8,1,1])
%     CHEQUEAR QUE EL SORTING SIEMPRE VAYA EN LA MISMA DIRECCION
    FRpi = spikes_downsample(spikes(sortind_descend,:),clus,1);
    

    %Prepares new spike matrix by keeping time bins with sequences only 
    
    clear IWI_ses
    for i=2:size(table_u,1)
        IWI_ses(i-1)=(table_u(i,1)-table_u(i-1,2))./sf;

        if (IWI_ses(i-1)>30)
            count=count+1;
            rate_before(count)=mean(mean(spikes(:,table_u(i,1)-ceil(win_length*sf) : table_u(i,1) )));
            rate_after(count)=mean(mean(spikes(:,table_u(i,1) : table_u(i,1)+ceil(win_length*sf))));
            percentage_change(count)=100*(rate_after(count)-rate_before(count))/rate_before(count);

            for ens=1:10
                rate_before_ens(ens,count)=(mean(FRpi(ens,table_u(i,1)-ceil(win_length*sf) : table_u(i,1) )));
                rate_after_ens(ens,count)=(mean(FRpi(ens,table_u(i,1) : table_u(i,1)+ceil(win_length*sf))));
            end
        end
    end

    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2 IWI_ses FRpi sorting_0 sortind_descend sorting_ascend
    
end

%% Version 2 with number of dt

clear all
close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

[big_table, waves] = get_big_table();

% IWI distribution
clus=10;
IWI_ses_pooled=[];
sf=7.73; 
win_length=5; %in number of dts
flipped_sorting(1,15)=0;
count=0;
n_dts=1;
bin_size=1/sf;

for w=1:length(waves)
    row_w=waves(w);
    disp(w);
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);
    munit=big_table(waves(w),5);
    
    dt=floor(big_table(waves(w),8));
    if isinteger(dt)
    else
        dt=floor(dt);
    end
    
    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_spk,'-mat');
    spikes=full(spikes_d_s);
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);

    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);        

    %Separation into ensembles
    [sorting_ascend,sortind_descend,sorting_0] = get_sorting(spikes) ;

    %Checks that the wave has an upwards direction
    clear FRpre
    if dt==0
        FRpre = spikes_downsample(spikes(sortind_descend,:),clus,1);
    else
        FRpre = spikes_downsample(spikes(sortind_descend,:),clus,dt);
    end

    clear signal; [~,signal]=max(FRpre); %Calculate signal by keeping the cluster with largest activity
    clear signal_s; signal_s(1)=signal(1);
    signal_s(2:length(find(diff(signal)~=0))+1) = signal(find(diff(signal)~=0) + 1);
    clear aux; aux=diff(signal_s);
    clear pos neg; pos=length(find(aux>6)); neg=length(find(aux<-6));
    if pos>neg; flipped_sorting(w)=1; sortind_descend=flip(sortind_descend); end

%     figure
%     spy(spikes(sortind_descend,:));
%     pbaspect([8,1,1])
%     CHEQUEAR QUE EL SORTING SIEMPRE VAYA EN LA MISMA DIRECCION
    FRpi = spikes_downsample(spikes(sortind_descend,:),clus,1);
    

    %Prepares new spike matrix by keeping time bins with sequences only 
    
    clear IWI_ses
    for i=2:size(table_u,1)
        IWI_ses(i-1)=(table_u(i,1)-table_u(i-1,2))./sf;

        if (IWI_ses(i-1)>35)
            count=count+1;
            rate_before(count)=mean(mean(spikes(:,table_u(i,1)-ceil(n_dts*dt) : table_u(i,1) )));
            rate_after(count)=mean(mean(spikes(:,table_u(i,1) : table_u(i,1)+ceil(n_dts*dt))));
            percentage_change(count)=100*(rate_after(count)-rate_before(count))/rate_before(count);
            change(count)= rate_after(count)./bin_size - rate_before(count)./bin_size;

            for ens=1:10
                rate_before_ens(ens,count)=(mean(FRpi(ens,table_u(i,1)-ceil(win_length*sf) : table_u(i,1) )));
                rate_after_ens(ens,count)=(mean(FRpi(ens,table_u(i,1) : table_u(i,1)+ceil(win_length*sf))));
                change_ens(ens,count)= rate_after_ens(ens,count)./bin_size - rate_before_ens(ens,count)./bin_size;
            end
        end
    end

    clear adj adj_sh freq_data_up freq_data_down_sh freq_data_up_sh freq_data_down freq_data_up mat_sh spikes sorting signal signal signal_s ...
        sorting mat adj_sh_th aux_mat aux_mat_sh signal_sh signal_sh_s spikes_d_s subset_2_sh freq_sum_up freq_sum_down coeff2 ...
        cells_d angle22 aux2 FRp thr_mat vals mat_sh table_u spikes_w aux_mat_2 IWI_ses FRpi sorting_0 sortind_descend sorting_ascend
    
end

% Figures and quantifications
%For all the network
bin_size=1/sf;
figure
boxplot([rate_before'./bin_size,rate_after'/bin_size]);
xticks([1,2])
xticklabels({'Before sequence onset', 'After sequence onset' });
ylabel('Mean calcium event rate (Hz)');
set(gca,'fontsize',18)
ylim([0 0.5]);
title('All neurons')
% title(['Window length: ',num2str(win_length),' s']);
set(gcf,'units','points','position',[500,300,500,400])
box off
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\RL\CalciumRate\CalciumRate_AllNeurons','_',num2str(n_dts),'dt.png'));
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\RL\CalciumRate\CalciumRate_AllNeurons','_',num2str(n_dts),'dt.fig'));
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\RL\CalciumRate\CalciumRate_AllNeurons','_',num2str(n_dts),'dt.svg'));
close all


figure
boxplot(change);
xticklabels({'All neurons' });
ylabel({'Change in mean';'calcium event rate (Hz)'});
set(gca,'fontsize',16)
box off
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\firing rate before and after sequence onset\CHANGECalciumRate_AllNeurons','_',num2str(n_dts),'dt.png'));
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\firing rate before and after sequence onset\CHANGECalciumRate_AllNeurons','_',num2str(n_dts),'dt.fig'));
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\firing rate before and after sequence onset\CHANGECalciumRate_AllNeurons','_',num2str(n_dts),'dt.svg'));
close all


[p_all,h_all,stat_all]=signrank(change);


[p,h,stat]=signrank(rate_before,rate_after);
disp([])
disp(['p value = ',num2str(p)])
disp(['Z = ',num2str(stat.zval)])
disp(['counts = ',num2str(count)])

%ensembles

for i=[1,5,10]
    figure
    boxplot([rate_before_ens(i,:)'./bin_size,rate_after_ens(i,:)'./bin_size])
    xticks([1,2])
    xticklabels({'Before sequence onset', 'After sequence onset' });
    ylabel('Mean calcium event rate (Hz)');
    set(gca,'fontsize',16)
    ylim([-0.01 1]);
    % title(['Window length: ',num2str(win_length),' s']);
    set(gcf,'units','points','position',[500,300,500,400])
    box off
    title(['Ensemble = ',num2str(i)]);
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\RL\CalciumRate\CalciumRate_Ensemble',num2str(i),'_',num2str(n_dts),'dt.png'));
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\RL\CalciumRate\CalciumRate_Ensemble',num2str(i),'_',num2str(n_dts),'dt.fig'));
    saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\RL\CalciumRate\CalciumRate_Ensemble',num2str(i),'_',num2str(n_dts),'dt.svg'));
    close all

    [p_ens(i),h_ens(i),stat_ens]=signrank(rate_before_ens(i,:),rate_after_ens(i,:));
    disp(['p value = ',num2str(p)])
    disp(['Z = ',num2str(stat.zval)])
    disp(['counts = ',num2str(count)])
end

%All ensembles
figure
boxplot(change_ens')
ylabel({'Change in mean';'calcium event rate (Hz)'});
xlabel('Ensemble #');
set(gca,'fontsize',16)
box off

saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\firing rate before and after sequence onset\CHANGECalciumRate_Ensembles','_',num2str(n_dts),'dt.png'));
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\firing rate before and after sequence onset\CHANGECalciumRate_Ensembles','_',num2str(n_dts),'dt.fig'));
saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\firing rate before and after sequence onset\CHANGECalciumRate_Ensembles','_',num2str(n_dts),'dt.svg'));


for i=1:10
    [p_ens(i),h_ens(i),stat_ens]=signrank(change_ens(i,:));
    [p_ens_right(i),h_ens_right(i),stat_ens_right]=signrank(change_ens(i,:),[],'tail','right');
    [p_ens_left(i),h_ens_left(i),stat_ens_left]=signrank(change_ens(i,:),[],'tail','left');
end