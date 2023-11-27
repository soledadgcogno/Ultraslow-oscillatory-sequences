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


%% Wave identification

downsampling_factor=4;
fs_120 = 7.73;
fs_30 = 30.95;
fs = fs_120/downsampling_factor;
maxlag=floor(560*fs_120);%240;%250 - 520;
N_sh=200;
count_c=0;

epoch_length=1;%1/8; %in seconds
win_size=ceil(fs_120)*epoch_length;

%Identifying waves
for w=1:length(waves)
    disp(w);
    row_w=waves(w);
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);
    munit=big_table(waves(w),5);
    dt=floor(big_table(waves(w),8));
    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_spk,'-mat');
    spikes=full(spikes_d_s);

    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);
    number_of_cells(w)=N;
    number_of_sequences(w)=size(table_u,1);

    %Prepares new spike matrix by keeping time bins with sequences only
    spikes_w=[];
    for i=1:size(table_u,1)
        spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
    end
    full_length(w)=size(spikes,2);
    wave_length(w)=size(spikes_w,2);

    %Calculates PSDs  
    FR=spikes;
    count=0;

    for i=1:N
        count_c=count_c+1;
        [alfa(i,:),lags]=xcorr(spikes_w(i,:),spikes_w(i,:),maxlag,'coeff');
        [val(i,:),~]=zscore(alfa(i,:)); %Check whether I need to zscore
        signal=alfa(i,:);

        %Calculate spectogram using pwelch method
        clear Powerspec2 Powerfreq2; [Powerspec2,Powerfreq2] = doPwelch((signal),fs_120,2*4096);
        [peak_freq,peak_psd,quality]=check_peak_quality_4_f(Powerspec2,Powerfreq2,0.04,find(Powerfreq2>0.003,1),find(Powerfreq2>0.1,1));
%         clear aux; aux=find(lags>20,1);
%         max(signal(aux:end))
       
% % %         Shuffling - Method 1
% %         clear spk; spk=spikes_w(i,:);
% %         clear Powerspec2_sh Powerfreq2_sh;
% %         for sh=1:N_sh
% %             random_ordering=randperm(length(spk));
% %             clear signal_sh; [signal_sh]=xcorr(spk(random_ordering),spk(random_ordering),maxlag,'coeff');
% %             [Powerspec2_sh(sh,:),Powerfreq2_sh(sh,:)] = doPwelch(signal_sh,fs_120,2*4096);
% % %             [peak_freq_sh(i,sh),peak_psd_sh(i,sh),quality_sh(i,sh)]=check_peak_quality_3c_f(Powerspec2_sh(sh,:)',Powerfreq2_sh(sh,:)',0.04);       
% %         end
% %         ub=max(find(Powerfreq2<0.1));
% %         pct=100-(5/ub);
% %         thr=prctile(Powerspec2_sh,pct,1);
% %         cell_ind(count_c)=0; 
% %         for l=1:ub 
% %             if (Powerspec2(l)>thr(l)) cell_ind(count_c)=1;    end
% %         end

        %Shuffling - Method 2
        clear spk; spk=spikes_w(i,:);
        n_chunks=floor(size(spikes_w,2)/win_size);
        spk_d=discretize([1:size(spikes_w,2)],n_chunks);
       
        for sh=1:N_sh
            clear spk_sh sh_o;
            spk_sh=[];
            sh_o=randperm(n_chunks);
            for c=1:n_chunks
                clear aux; aux=find(spk_d==sh_o(c));
%                 spk_sh((c-1)*length(aux) +1 : c*length(aux))=spk(aux);
                spk_sh=[spk_sh,spk(aux)];
            end

            clear signal_sh; [signal_sh]=xcorr(spk_sh,spk_sh,maxlag,'coeff');
            [Powerspec2_sh(sh,:),Powerfreq2_sh] = doPwelch((signal_sh),fs_120,2*4096);
        end


        clear aux; aux=find(Powerfreq2==peak_freq);

        thr95=prctile(Powerspec2_sh,95,1);
        thr99=prctile(Powerspec2_sh,99,1);

        cell_ind_95(count_c)=0;
        cell_ind_99(count_c)=0;
        cell_ind_win1(count_c)=0;
        if (Powerspec2(aux)>thr95(aux)) cell_ind_95(count_c)=1;    end
        if (Powerspec2(aux)>thr99(aux)) cell_ind_99(count_c)=1;    end
        if (peak_freq>0.1) cell_ind_win1(count_c)=1;    end
      


        %
        %         for win=1:3
%             ub=aux-2:aux+2;
%             pct=100-(5/length(ub));
%             thr=prctile(Powerspec2_sh,pct,1);
%             cell_ind_win1(count_c)=0;
%             for l=ub
%                 %                 disp(l)
%                 if l>0  if (Powerspec2(l)>thr(l)) cell_ind_win1(count_c)=1;    end ;end
%             end
% 
%             ub=aux-5:aux+5;
%             pct=100-(5/length(ub));
%             thr=prctile(Powerspec2_sh,pct,1);
%             cell_ind_win2(count_c)=0;
%             for l=ub
%                 %                 disp(l)
%                 if l>0 if (Powerspec2(l)>thr(l)) cell_ind_win2(count_c)=1;    end; end
%             end
% 
%             ub=aux-10:aux+10;
%             pct=100-(5/length(ub));
%             thr=prctile(Powerspec2_sh,pct,1);
%             cell_ind_win3(count_c)=0;
%             for l=ub
%                 if l>0
%                     if (Powerspec2(l)>thr(l)) cell_ind_win3(count_c)=1;    end
%                 end
%             end
% 
%         end


        clear  Powerspec2 P2 P1  Y Powerfreq2 signal laf alfa val thr spk Powerspec2_sh Powerfreq2_sh spk_d
    end

    clear table_u spikes dt spikes_d_s row_w mat_b spikes_W Ens_b sorting_before FR alfa val lags ...
        spikes_w table_u count N Powerspec2 Powerfreq2 spk cells_d
end