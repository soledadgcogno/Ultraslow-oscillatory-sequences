function [fraction,osc,tot] = get_oscillatory_cells_npx(spikes_w,bin_size,epoch_length,prtile)

fs_120 = 1/bin_size;
maxlag=floor(250*fs_120);%240;%250 - 520;
N_sh=200;
count_c=0;
N=size(spikes_w,1);
% epoch_length=20; %in seconds
win_size=floor(fs_120)*epoch_length;

for i=1:N
    count_c=count_c+1;
    clear signal; [signal,~]=xcorr(spikes_w(i,:),spikes_w(i,:),floor(maxlag),'coeff');

    %Calculate spectogram using pwelch method
    clear Powerspec2 Powerfreq2; [Powerspec2,Powerfreq2] = doPwelch(signal,fs_120,1*4096);
    [peak_freq,~,~]=check_peak_quality_4_f(Powerspec2,Powerfreq2,0.04,find(Powerfreq2>0.003,1),find(Powerfreq2>0.1,1));

    %Shuffling
    clear spk; spk=spikes_w(i,:);
    n_chunks=floor(size(spikes_w,2)/win_size);
    spk_d=discretize([1:size(spikes_w,2)],n_chunks);
    for sh=1:N_sh
        clear spk_sh sh_o;
        spk_sh=[];
        sh_o=randperm(n_chunks);
        for c=1:n_chunks
            clear aux; aux=find(spk_d==sh_o(c));
            spk_sh=[spk_sh,spk(aux)];
        end

        clear signal_sh; [signal_sh]=xcorr(spk_sh,spk_sh,maxlag,'coeff');
        [Powerspec2_sh(sh,:),Powerfreq2_sh] = doPwelch((signal_sh),fs_120,1*4096);
    end


    clear aux; aux=find(Powerfreq2==peak_freq);

    thr=prctile(Powerspec2_sh,prtile,1);
%     thr99=prctile(Powerspec2_sh,99,1);

    cell_ind(count_c)=0;
%     cell_ind_99(count_c)=0;
    cell_ind_win1(count_c)=0;
    if (Powerspec2(aux)>thr(aux)) cell_ind(count_c)=1;    end
%     if (Powerspec2(aux)>thr99(aux)) cell_ind_99(count_c)=1;    end
    if (peak_freq>0.1) cell_ind_win1(count_c)=1;    end

    fraction=sum(cell_ind)/length(cell_ind);
    osc=sum(cell_ind);
    tot=length(cell_ind);

    clear  Powerspec2 P2 P1  Y Powerfreq2 signal laf alfa val thr spk Powerspec2_sh Powerfreq2_sh spk_d
end
end