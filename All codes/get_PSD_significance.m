function [Powerspec2, Powerfreq2,peak_freq,peak_psd,Powerspec2_sh,alfa] = get_PSD_significance(spikes_w, maxlag,fs_120,win_size,N_sh)

[alfa,lags]=xcorr(spikes_w,spikes_w,maxlag,'coeff');
[val,~]=zscore(alfa); %Check whether I need to zscore
signal=alfa;

%Calculate spectogram using pwelch method
clear Powerspec2 Powerfreq2; [Powerspec2,Powerfreq2] = doPwelch((signal),fs_120,2*4096);
[peak_freq,peak_psd,quality]=check_peak_quality_4_f(Powerspec2,Powerfreq2,0.04,find(Powerfreq2>0.003,1),find(Powerfreq2>0.1,1));

%Shuffling - Method 2
clear spk; spk=spikes_w;
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
    [Powerspec2_sh(sh,:),Powerfreq2_sh] = doPwelch((signal_sh),fs_120,2*4096);
end

end