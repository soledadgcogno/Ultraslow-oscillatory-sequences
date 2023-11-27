%% Autocorrelations
clear all
close all

dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
plot_all=1;
[big_table, waves] = get_big_table();

for w=1:length(waves)
    row_w=waves(w);
    disp(w);
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);
    munit=big_table(waves(w),5);
    load([dpath,strcat('recording_dates_',mouse,'.mat')]);
    file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

    load(file_name,'-mat');
    spikes=full(spikes_d_s);
    downsampling_factor=4;
    fs_120 = 7.73;
    fs_30 = 30.95;
    fs = fs_120/downsampling_factor;
    maxlag=floor(560*fs_120);%240;%250 - 520;

    [N,T]=size(spikes);
    FRp = spikes_downsample(spikes,N,downsampling_factor);

    FR=spikes;
    count=0;
 
    currentFolder = sprintf('C:\\Users\\xscogno\\Dropbox\\Paper Waves\\Figures Waves Paper\\Revision - 1st round\\Extended 3\\session%d',w);
    mkdir(currentFolder)

    for i=1:N
        [alfa(i,:),lags]=xcorr(FR(i,:),FR(i,:),maxlag,'coeff');
        [val(i,:),~]=zscore(alfa(i,:)); %Check whether I need to zscore
        signal=alfa(i,:);
        clear Powerspec2 Powerfreq2
        %Calculate spectogram using pwelch method
        [Powerspec2,Powerfreq2] = doPwelch(signal,fs_120,2*4096);
        [Powerspec2z,Powerfreq2z] = doPwelch(zscore(signal),fs_120,2*4096);
        pwelch_fft(i,:)=Powerspec2;
        pwelch_fft_z(i,:)=Powerspec2z;
        pwelch_fft_section(i,:)=Powerspec2(5:100);
        %Check if the cell is oscillatory using the pwelch
        %method
        [peak_freq(i),peak_psd(i),quality(i)]=check_peak_quality_3c_f(Powerspec2,Powerfreq2,0.04);
        %         period_fft(i)=1/peak_fft(i);
        period_ps(i)=1/peak_freq(i);
        if quality(i)==1
            count=count+1;
            cells_osc(count)=i;
        else
        end
        clear  Powerspec2 P2 P1  Y

        %1 shuffle - circular
        rand_num_circ=floor((T-120/(1/fs_120))*rand) + floor(120/(1/fs_120)); % I shifht by at least 2 minutes
        spk_circ_shuffle=circshift(FR(i,:),rand_num_circ);
        alfa_circshuffle(i,:)=xcorr(spk_circ_shuffle,spk_circ_shuffle,maxlag,'coeff');
        alfa_circshuffle_z(i,:)=zscore(alfa_circshuffle(i,:));
        [Powerspec2_circshuffle(i,:),~] = doPwelch(alfa_circshuffle(i,:),fs_120,2*4096);
        [Powerspec2_circshuffle_z(i,:),~] = doPwelch(zscore(alfa_circshuffle(i,:)),fs_120,2*4096);
        
        %1 shuffle - total
        spk_tot_shuffle=shuffle(FR(i,:));
        alfa_totshuffle(i,:)=xcorr(spk_tot_shuffle,spk_tot_shuffle,maxlag,'coeff');
        alfa_totshuffle_z(i,:)=zscore(alfa_totshuffle(i,:));
        [Powerspec2_totshuffle(i,:),~] = doPwelch(alfa_totshuffle(i,:),fs_120,2*4096);
        [Powerspec2_totshuffle_z(i,:),~] = doPwelch(zscore(alfa_totshuffle(i,:)),fs_120,2*4096);
        
        clear spk_tot_shuffle
    end

    [~,sorting_ps2]=sort(peak_psd,'descend');
    %     [val_sorting_freq2,sorting_freq2]=sort(peak_freq,'descend');
    %     [sorting_ascend,sortind_descend,sorting_0] = get_sorting(spikes);

    %Stacked autocorrelations
    figure
    imagesc(val(sorting_ps2,:))
    caxis([0 0.5])
    colormap cividis;
    set(gca,'fontsize',18)
    xticks([1 maxlag maxlag*2])
    xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
    xlabel('Time (s)');
    ylabel('Neurons #');
    yticks([100 400])
    axis square
    colorbar
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\preliminary\session',num2str(w),'\DataAutoFull_session',num2str(w),'.png'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\DataAutoFull_session',num2str(w),'.png'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\DataAutoFull_session',num2str(w),'.svg'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\DataAutoFull_session',num2str(w),'.fig'));

    %Mean over cells
%     figure
%     sgtitle('Mean over cells')
%     subplot(1,2,1)
%     plot(mean(alfa,1),'k','linewidth',2.5);
%     xticks([1 maxlag maxlag*2])
%     xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
%     xlabel('Time lag (s)');
%     ylabel('Autocorrelation');
%     axis square
%     axis([-inf inf 0 0.1])
%     set(gca,'fontsize',10);
%     subplot(1,2,2)
%     plot(Powerfreq2(1:65),mean(pwelch_fft(:,1:65),1),'k','linewidth',2.5);
%     ylabel('PSD (mV^2/Hz)');
%     % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
%     % title(MVL(1))
%     % ylabel('PSD (mV^2/Hz)');
%     % ylabel('PSD ');
%     xlabel('Frequency (Hz)');
%     set(gca,'fontsize',10,'YColor','k','XColor','k');
%     axis square
% %     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\preliminary\session',num2str(w),'\MeanPop_session',num2str(w),'.png'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\MeanPop_session',num2str(w),'.png'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\MeanPop_session',num2str(w),'.png'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\MeanPop_session',num2str(w),'.png'));

    %Mean over cells - Calculating zscore
    lr=[255, 204, 203]/255;
    lg=[169,169,169]/255;
    figure
    sgtitle('Mean over cells - zscore')
    hold on
    errorbar(1:length(lags),mean(alfa_totshuffle_z,1),std(alfa_totshuffle_z,1)/sqrt(N),'color',lr,'linewidth',0.5);
    plot(1:length(lags),mean(alfa_totshuffle_z,1),'r','linewidth',2.5);
    xticks([1 maxlag maxlag*2])
    xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
    xlabel('Time lag (s)');
    ylabel('Z-scored autocorrelation');
    axis square
    axis([-inf inf -0.5 1])
    set(gca,'fontsize',16,'YColor','k','XColor','k');
    hold on
    errorbar(1:length(lags),mean(val,1),std(val,1)./sqrt(N),'color',lg,'linewidth',0.5);
    plot(1:length(lags),mean(val,1),'k','linewidth',2.5);
    xticks([1 maxlag maxlag*2])
    xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
    xlabel('Time lag (s)');
    axis square
    axis([-inf inf -0.5 1])
    set(gca,'fontsize',16,'YColor','k','XColor','k');
    set(gcf,'renderer','painters');
    hold on
    box off
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\preliminary\session',num2str(w),'\MeanAutoZ_session',num2str(w),'.png'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\MeanAutoZ_session',num2str(w),'.png'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\MeanAutoZ_session',num2str(w),'.svg'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\MeanAutoZ_session',num2str(w),'.fig'));
    saveas(gcf,strcat(currentFolder,'\MeanAutoZ_painter.png'));
    saveas(gcf,strcat(currentFolder,'\MeanAutoZ_painter.svg'));
    saveas(gcf,strcat(currentFolder,'\MeanAutoZ_painter.fig'));

    figure
    hold on
%     plot(Powerfreq2(1:65),mean(pwelch_fft_z(:,1:65),1),'k','linewidth',2.5);
    errorbar(Powerfreq2(1:65),mean(pwelch_fft_z(:,1:65),1),std(pwelch_fft_z(:,1:65),1)./sqrt(N),'color',lg,'linewidth',1.5);
    plot(Powerfreq2(1:65),mean(pwelch_fft_z(:,1:65),1),'k','linewidth',2.5);
    ylabel('Z-scored PSD (mV^2/Hz)');
    % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
    % title(MVL(1))
    % ylabel('PSD (mV^2/Hz)');pwelch_fft_z
    % ylabel('PSD ');
    xlabel('Frequency (Hz)');
    hold on
    errorbar(Powerfreq2(1:65),mean(Powerspec2_totshuffle_z(:,1:65),1),std(Powerspec2_totshuffle_z(:,1:65)./sqrt(N),1),'color',lr,'linewidth',1.5);
    plot(Powerfreq2(1:65),mean(Powerspec2_totshuffle_z(:,1:65),1),'r','linewidth',2.5);
    set(gca,'fontsize',16,'YColor','k','XColor','k');
    axis square
    box off
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\preliminary\session',num2str(w),'\MeanPopZ_session',num2str(w),'.png'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\MeanPSDZ_session',num2str(w),'.png'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\MeanPSDZ_session',num2str(w),'.svg'));
%     saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\MeanPSDZ_session',num2str(w),'.fig'));
%     saveas(gcf,strcat(currentFolder,'\MeanPSDZ.png'));
%     saveas(gcf,strcat(currentFolder,'\MeanPSDZ.svg'));
%     saveas(gcf,strcat(currentFolder,'\MeanPSDZ.fig'));

    for i=1:10

        currentFolder_2 = sprintf('C:\\Users\\xscogno\\Dropbox\\Paper Waves\\Figures Waves Paper\\Revision - 1st round\\Extended 3\\session%d\\cell%d',w,i);
        mkdir(currentFolder_2)


        %Autocorrelation - Exp. data
        figure
        plot((alfa(sorting_ps2(i),:)),'k','linewidth',2.5)
        xticks([1 maxlag maxlag*2])
        xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
        xlabel('Time lag (s)');
        ylabel('Autocorrelation');
        axis([-inf inf 0 0.3])
        set(gca,'fontsize',20,'YColor','k','XColor','k');
        title(['Cell ',num2str(sorting_ps2(i))]);
        % h = gca; h.YAxis.Visible = 'off';
        axis square
        box off
%         saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\preliminary\session',num2str(w),'\DataAuto_session',num2str(w),'_cell',num2str((i)),'.png'));

        if plot_all==1
%             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\DataAuto_session',num2str(w),'_cell',num2str((i)),'.png'));
%             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\DataAuto_session',num2str(w),'_cell',num2str((i)),'.fig'));
%             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\DataAuto_session',num2str(w),'_cell',num2str((i)),'.svg'));
%             saveas(gcf,strcat(currentFolder_2,'\DataAuto.png'));
%             saveas(gcf,strcat(currentFolder_2,'\DataAuto.svg'));
%             saveas(gcf,strcat(currentFolder_2,'\DataAuto.fig'));

            %PSD - Exp. Data
            figure
            plot(Powerfreq2(1:65),(pwelch_fft(sorting_ps2(i),1:65)),'k','linewidth',2.5);
            hold on
            l=xline(peak_freq(sorting_ps2(i)));
            l.LineStyle='--';
            l.Color=[170 170 170]/255;
            l.LineWidth=1.5;
            % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
            % title(MVL(1))
            ylabel('PSD (mV^2/Hz)');
            % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
            % title(MVL(1))
            % ylabel('PSD (mV^2/Hz)');
            % ylabel('PSD ');
            xlabel('Frequency (Hz)');
            set(gca,'fontsize',20,'YColor','k','XColor','k');
            box off
            title(['Cell ',num2str(sorting_ps2(i))]);
            % h=xline(freq_wave,'--r','linewidth',2.5);
            % h=xline(freq_wave/2,'--b','linewidth',2.5);
            axis square
%             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\DataPSD_session',num2str(w),'_cell',num2str((i)),'.png'));
%             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\DataPSD_session',num2str(w),'_cell',num2str((i)),'.fig'));
%             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\DataPSD_session',num2str(w),'_cell',num2str((i)),'.svg'));
%             saveas(gcf,strcat(currentFolder_2,'\DataPSD.png'));
%             saveas(gcf,strcat(currentFolder_2,'\DataPSD.svg'));
%             saveas(gcf,strcat(currentFolder_2,'\DataPSD.fig'));

            %Circular shuffle
            rand_num_circ=floor((T-120/(1/fs_120))*rand) + floor(120/(1/fs_120)); % I shifht by at least 2 minutes
            spk=FR(sorting_ps2(i),:);
            spk_circ_shuffle=circshift(spk,rand_num_circ);
            alfa_circshuffle=xcorr(spk_circ_shuffle,spk_circ_shuffle,maxlag,'coeff');
            [Powerspec2_circshuffle,~] = doPwelch(alfa_circshuffle,fs_120,2*4096);
            %             figure
            %             plot(alfa_circshuffle,'k','linewidth',2)
            %             xticks([1 maxlag maxlag*2])
            %             xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
            %             xlabel('Time lag (s)');
            %             ylabel('Autocorrelation');
            %             axis([-inf inf 0 0.3])
            %             set(gca,'fontsize',20);
            %             title(['Cell ',num2str(sorting_ps2(i))]);
            %             axis square
            %             box off
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CircShAuto_session',num2str(w),'_cell',num2str((i)),'.png'));
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CircShAuto_session',num2str(w),'_cell',num2str((i)),'.fig'));
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CircShAuto_session',num2str(w),'_cell',num2str((i)),'.svg'));
            %             figure
            %             plot(Powerfreq2(1:65),Powerspec2_circshuffle(1:65),'k','linewidth',2.5);
            %             hold on
            %             l=xline(peak_freq(sorting_ps2(i)));
            %             l.LineStyle='--';
            %             l.Color=[170 170 170]/255;
            %             l.LineWidth=1.5;
            %             ylabel('PSD (mV^2/Hz)');
            %             xlabel('Frequency (Hz)');
            %             set(gca,'fontsize',20,'YColor','k','XColor','k');
            %             box off
            %             title(['Cell ',num2str(sorting_ps2(i))]);
            %             axis square
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CircShPSD_session',num2str(w),'_cell',num2str((i)),'.png'));
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CircShPSD_session',num2str(w),'_cell',num2str((i)),'.fig'));
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CircShPSD_session',num2str(w),'_cell',num2str((i)),'.svg'));

            %Total shuffle
            spk_full_shuffle=shuffle(spk);
            alfa_totshuffle=xcorr(spk_full_shuffle,spk_full_shuffle,maxlag,'coeff');
            [Powerspec2_totshuffle,~] = doPwelch(alfa_totshuffle,fs_120,2*4096);
            %             figure
            %             plot(alfa_totshuffle,'k','linewidth',2)
            %             xticks([1 maxlag maxlag*2])
            %             xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
            %             xlabel('Time lag (s)');
            %             ylabel('Autocorrelation');
            %             axis([-inf inf 0 0.3])
            %             set(gca,'fontsize',20);
            %             title(['Cell ',num2str(sorting_ps2(i))]);
            %             axis square
            %             box off
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\TotShAuto_session',num2str(w),'_cell',num2str((i)),'.png'));
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\TotShAuto_session',num2str(w),'_cell',num2str((i)),'.fig'));
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\TotShAuto_session',num2str(w),'_cell',num2str((i)),'.svg'));
            %             figure
            %             plot(Powerfreq2(1:65),Powerspec2_totshuffle(1:65),'k','linewidth',2.5);
            %             hold on
            %             l=xline(peak_freq(sorting_ps2(i)));
            %             l.LineStyle='--';
            %             l.Color=[170 170 170]/255;
            %             l.LineWidth=1.5;
            %             ylabel('PSD (mV^2/Hz)');
            %             xlabel('Frequency (Hz)');
            %             set(gca,'fontsize',20,'YColor','k','XColor','k');
            %             box off
            %             title(['Cell ',num2str(sorting_ps2(i))]);
            %             axis square
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\TotShPSD_session',num2str(w),'_cell',num2str((i)),'.png'));
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\TotShPSD_session',num2str(w),'_cell',num2str((i)),'.fig'));
            %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\TotShPSD_session',num2str(w),'_cell',num2str((i)),'.svg'));

            %Combined plots
            figure
            plot(alfa_totshuffle,'r','linewidth',2.5)
            xticks([1 maxlag maxlag*2])
            xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
            xlabel('Time lag (s)');
            ylabel('Autocorrelation');
            axis([-inf inf 0 0.3])
            set(gca,'fontsize',20);
            title(['Cell ',num2str(sorting_ps2(i))]);
            axis square
            box off
            hold on
            plot(alfa_circshuffle,'b','linewidth',2.5)
            xticks([1 maxlag maxlag*2])
            xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
            xlabel('Time lag (s)');
            ylabel('Autocorrelation');
            axis([-inf inf 0 0.3])
            set(gca,'fontsize',20,'YColor','k','XColor','k');
            title(['Cell ',num2str(sorting_ps2(i))]);
            axis square
            box off
%             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CombinedShAuto_session',num2str(w),'_cell',num2str((i)),'.png'));
%             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CombinedShAuto_session',num2str(w),'_cell',num2str((i)),'.fig'));
%             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CombinedShAuto_session',num2str(w),'_cell',num2str((i)),'.svg'));
%             saveas(gcf,strcat(currentFolder_2,'\CombinedShAuto.png'));
%             saveas(gcf,strcat(currentFolder_2,'\CombinedShAuto.svg'));
%             saveas(gcf,strcat(currentFolder_2,'\CombinedShAuto.fig'));

            figure
            plot(Powerfreq2(1:65),Powerspec2_totshuffle(1:65),'r','linewidth',2.5)
            %             xticks([1 maxlag maxlag*2])
            %             xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
            xlabel('Frequency (Hz)');
            ylabel('PSD (mV^2/Hz)');
            %             axis([-inf inf 0 0.3])
            set(gca,'fontsize',20);
            title(['Cell ',num2str(sorting_ps2(i))]);
            axis square
            box off
            hold on
            plot(Powerfreq2(1:65),Powerspec2_circshuffle(1:65),'b','linewidth',2.5)
            %             xticks([1 maxlag maxlag*2])
            %             xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
            xlabel('Frequency (Hz)');
            ylabel('PSD (mV^2/Hz)');
            %             axis([-inf inf 0 0.3])
            set(gca,'fontsize',20,'YColor','k','XColor','k');
            title(['Cell ',num2str(sorting_ps2(i))]);
            axis square
            box off
%             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CombinedShPSD_session',num2str(w),'_cell',num2str((i)),'.png'));
%             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CombinedShPSD_session',num2str(w),'_cell',num2str((i)),'.fig'));
% %             saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Waves\Revision\New figures\USOs in single cells and for all sessions\session',num2str(w),'\CombinedShPSD_session',num2str(w),'_cell',num2str((i)),'.svg'));
%             saveas(gcf,strcat(currentFolder_2,'\CombinedShPSD.png'));
%             saveas(gcf,strcat(currentFolder_2,'\CombinedShPSD.svg'));
%             saveas(gcf,strcat(currentFolder_2,'\CombinedShPSD.fig'));
        end
        close all
    end
    clear alfa_totshuffle alfa_circshuffle alfa cells_osc cells_d FR FRp lags peak_freq peak_psd period_ps Powerfreq2 Powerspec2_circshuffle ...
        Powerspec2_totshuffle pwelch_fft pwelch_fft_section quality signal sorting_0 sortind_descend sorting_ascend sorting_freq2 sorting_ps2 ...
        spikes spikes_d_s spk spk_full_shuffle spk_circ_shuffle val val_sorting_freq2 val_sorting_ps2 dates pwelch_fft_section pwelch_fft_z ...
        Powerspec2_circshuffle_z Powerfreq2z Powerspec2z Powerspec2_totshuffle_z alfa_circshuffle_z alfa_totshuffle_z
end

%% Figures autocorrelograms for all cells


%PCA sorting
figure
imagesc(val(sortind_descend,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 400])
axis square
colorbar
caxis([0 1.5])

%Frequency sorting
figure
imagesc(val(sorting_freq2,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 400])
axis square
colorbar


%Power sorting
figure
imagesc(val(sorting_ps2,:))
caxis([0 0.5])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 400])
axis square
colorbar

% title('PS sorting'); 

%PSD for oscillatory cells

aux=find(quality==0);
sorting_ps2_c=sorting_ps2;
sorting_ps2_c(aux)=[];

figure
imagesc(val(sorting_ps2_c,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([100 300])
axis square
hcb = colorbar;
hcb.Title.String = 'zscore(autocorrelation)';
title('PS sorting'); 

%PSD for non-oscillatory cells

aux2=find(quality==1);
sorting_ps2_c2=sorting_ps2;
sorting_ps2_c2(aux2)=[];

figure
imagesc(val(sorting_ps2_c2,:))
caxis([0 1])
colormap cividis;
set(gca,'fontsize',18)
xticks([1 maxlag maxlag*2])
xticklabels({floor(-maxlag/fs_120) 0 ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Neurons #');
yticks([10 30])
axis square
hcb = colorbar;
% hcb.Title.String = 'zscore(autocorrelation)';
% title('PS sorting'); 

% [~,sorting_w,~]=get_sorting(spikes);
% 
% 
% figure
% imagesc(val(sorting_w,:))
% caxis([0 1])
% colormap cividis;
% axis square
% set(gca,'fontsize',18)
% xticks([1 maxlag maxlag*2])
% xticklabels({-maxlag/2 0 maxlag/2});
% xlabel('Time (s)');
% ylabel('Neurons #');
% yticks([100 400])
% axis square
% colorbar
% title('PCA sorting'); 



%% Mean autocorrelation

%Mean autocorrelation
mean_auto=mean(alfa(find(quality>0),:));
std_auto=std(alfa(find(quality>0),:));
% std_auto=std(alfa(find(quality>0),:))/sqrt(count);

figure
x=1:size(alfa,2);
x2=[x,fliplr(x)];
inBetween = [mean_auto+std_auto, fliplr(mean_auto-std_auto)];
fill(x2, inBetween, [252 215 215]./255);
hold on
plot(mean_auto+std_auto,'w','linewidth',1.5);
hold on
plot(mean_auto-std_auto,'w','linewidth',1.5);
plot((mean(alfa(find(quality>0),:))),'k','linewidth',1);
xticks([1,length(mean_auto)/2,length(mean_auto)]);
xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
axis([-inf inf -inf 0.2])
xlabel('Time (s)')
ylabel('Mean autocorrelation');
% set(gca, 'YScale', 'log')
set(gca,'fontsize',16);


%Fraction of oscillatory cells
osc=count/N;
non_osc=(N-count)/N;

figure
bar([1,2],[osc,non_osc]*100);
axis([0.5 2.5 0 100]);
xticks([1,2]);
xticklabels({'Osc','Non-Osc'});
ylabel('Cells percentage %')
box off
set(gca,'fontsize',16)

%Period distribution for cells that oscillate
periods=1./peak_freq(cells_osc);
figure
histogram(periods,0:10:500);
ylabel('Counts');
xlabel('Cell period (s)');
set(gca,'fontsize',16);


% figure
% subplot(1,2,1)
% plot(peak_freq)
% title('Freq')
% subplot(1,2,2)
% plot(1./peak_freq)
% title('Period')
% 




%% Example of individual cells using PSD


% [table_u,~,~]=identify_waves_latestversion_6_f(mouse,day,10,dt,0,spikes);      
% 
% num_waves=size(table_u,1);
% freq_wave=num_waves/(T/fs_120);

[val_peak_psd,peak_psd_ind]=sort(peak_psd, 'descend');

N_non_osc=50;
N_osc=200;
for i=1:N_non_osc
    cell_non_osc(i)=peak_psd_ind(i); 
end

non_osc=N-count;
peak_psd_ind_c=peak_psd_ind;
peak_psd_ind_c(1:non_osc)=[];
val_peak_psd_c=val_peak_psd;
val_peak_psd_c(1:non_osc)=[];

for i=1:N_osc
    cell_osc(i)=peak_psd_ind_c(i); 
end

countnonosc=0;
for i=(find(quality==0))
    countnonosc=countnonosc+1;
    if countnonosc<30
        figure
        plot((alfa(i,:)),'k','linewidth',2)
        % set(gca, 'YScale', 'log')
        xticks([1 maxlag maxlag*2])
        xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
        xlabel('Time (s)');
        ylabel('Autocorrelation');
        axis([-inf inf 0 0.2])
        set(gca,'fontsize',20);
        % h = gca; h.YAxis.Visible = 'off';
        axis square
        box off
        figure
        plot(Powerfreq2(1:65),(pwelch_fft(i,1:65)),'k','linewidth',2.5);
        hold on
        l=xline(peak_freq(i));
        l.LineStyle='--';
        l.Color=[170 170 170]/255;
        l.LineWidth=1.5;
        % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
        % title(MVL(1))
        ylabel('PSD (mV^2/Hz)');
%         ylabel('PSD ');
        xlabel('Frequency (Hz)');
        set(gca,'fontsize',20,'YColor','k','XColor','k');
        box off
        % h=xline(freq_wave,'--r','linewidth',2.5);
        % h=xline(freq_wave/2,'--b','linewidth',2.5);
        axis square
    end
end

%Cells with high peak in PSD
for i=1:90
% figure
% subplot(1,2,1)
% plot((alfa(cell_osc(i),:)),'k','linewidth',2)
% xticks([1 maxlag maxlag*2])
% xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
% xlabel('Time (s)');
% ylabel('Autocorrelation');
% axis([-inf inf 0 0.3])
% set(gca,'fontsize',20);
% % h = gca; h.YAxis.Visible = 'off';
% axis square
% box off
figure
plot(Powerfreq2(1:65),(pwelch_fft(cell_osc(i),1:65)),'k','linewidth',2.5);
hold on
 l=xline(peak_freq(cell_osc(i)));
        l.LineStyle='--';
        l.Color=[170 170 170]/255;
        l.LineWidth=1.5;
        % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
        % title(MVL(1))
        ylabel('PSD (mV^2/Hz)');
% plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
% title(MVL(1))
% ylabel('PSD (mV^2/Hz)');
% ylabel('PSD ');
xlabel('Frequency (Hz)');
set(gca,'fontsize',20,'YColor','k','XColor','k');
box off
% h=xline(freq_wave,'--r','linewidth',2.5);
% h=xline(freq_wave/2,'--b','linewidth',2.5);
axis square
% title(round(mvl(cell_osc(i)),2))
end

%Isolated examples

high_freq=find(peak_freq>0.02);

for i=high_freq
figure
% subplot(1,2,1)
plot((alfa(i,:)),'k','linewidth',2)
xticks([1 maxlag maxlag*2])
xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Autocorrelation');
axis([-inf inf 0 0.3])
set(gca,'fontsize',20);
% h = gca; h.YAxis.Visible = 'off';
axis square
box off
figure
plot(Powerfreq2(1:65),(pwelch_fft(i,1:65)),'k','linewidth',2.5);
hold on
% plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
% title(MVL(1))
hold on
 l=xline(peak_freq(i));
        l.LineStyle='--';
        l.Color=[170 170 170]/280;
        l.LineWidth=1.5;
        % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
        % title(MVL(1))
        ylabel('PSD (mV^2/Hz)');% ylabel('PSD ');
xlabel('Frequency (Hz)');
set(gca,'fontsize',20,'YColor','k','XColor','k');
box off
% h=xline(freq_wave,'--r','linewidth',2.5);
% h=xline(freq_wave/2,'--b','linewidth',2.5);
axis square
% title(round(mvl(cell_osc(i)),2))
end


high_freq=find(peak_freq<0.0067);

for i=129%101:150%high_freq
figure
% subplot(1,2,1)
plot((alfa(high_freq(i),:)),'k','linewidth',2)
xticks([1 maxlag maxlag*2])
xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
xlabel('Time (s)');
ylabel('Autocorrelation');
axis([-inf inf 0 0.3])
set(gca,'fontsize',20);
% h = gca; h.YAxis.Visible = 'off';
axis square
box off
figure
plot(Powerfreq2(1:65),(pwelch_fft(high_freq(i),1:65)),'k','linewidth',2.5);
hold on
% plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
% title(MVL(1))
hold on
 l=xline(peak_freq(high_freq(i)));
        l.LineStyle='--';
        l.Color=[170 170 170]/280;
        l.LineWidth=1.5;
        % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
        % title(MVL(1))
        ylabel('PSD (mV^2/Hz)');% ylabel('PSD ');
xlabel('Frequency (Hz)');
set(gca,'fontsize',20,'YColor','k','XColor','k');
box off
% h=xline(freq_wave,'--r','linewidth',2.5);
% h=xline(freq_wave/2,'--b','linewidth',2.5);
axis square
% title(round(mvl(cell_osc(i)),2))
end


for i=144%0:150%high_freq
% figure
% % subplot(1,2,1)
% plot((alfa(high_freq(i),:)),'k','linewidth',2)
% xticks([1 maxlag maxlag*2])
% xticklabels({-ceil(maxlag/fs_120),0,ceil(maxlag/fs_120)});
% xlabel('Time (s)');
% ylabel('Autocorrelation');
% axis([-inf inf 0 0.3])
% set(gca,'fontsize',20);
% % h = gca; h.YAxis.Visible = 'off';
% axis square
% box off
figure
plot(Powerfreq2(1:65),(pwelch_fft(high_freq(i),1:65)),'k','linewidth',2.5);
hold on
% plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
% title(MVL(1))
hold on
 l=xline(peak_freq(high_freq(i)));
        l.LineStyle='--';
        l.Color=[170 170 170]/280;
        l.LineWidth=1.5;
        % plot(f(1:70),fft_(cell_high_locking,1:70),'k','linewidth',2.5);
        % title(MVL(1))
        ylabel('PSD (mV^2/Hz)');% ylabel('PSD ');
        ylim([0 0.8])
xlabel('Frequency (Hz)');
set(gca,'fontsize',20,'YColor','k','XColor','k');
box off
% h=xline(freq_wave,'--r','linewidth',2.5);
% h=xline(freq_wave/2,'--b','linewidth',2.5);
axis square
% title(round(mvl(cell_osc(i)),2))
end

