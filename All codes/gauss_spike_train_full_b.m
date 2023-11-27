function [sd,bin] = gauss_spike_train_full_b(spk_times,start,stop,dt,thr_SD,downsample_new_bin)

%kernel_size has to be in seconds
%interval is in seconds

N=length(spk_times);
startT=start*1000; % start in ms
stopT=stop*1000; % stop in ms
bin_size=120;
dt_bin=dt/(bin_size*0.001);
trialTime=stopT-startT;
nBins=floor(trialTime/bin_size); 
interval=floor(downsample_new_bin*1000/bin_size);


for iCell=1:length(spk_times)
    normSpk=spk_times{iCell}.*1000; % sec to ms
    normSpk=normSpk-startT; % normalise to start of trial
    normSpk=normSpk(normSpk<=stopT);
    edgesT=linspace(0, stopT-startT, nBins+1);
    binnedSpks=histcounts(normSpk,edgesT);
    binnedSpikes(iCell,1)={binnedSpks};
end

% new_bin=interval/bin_size;
for i=1:N
    s(i,:) = full(fire_rate(binnedSpikes{i,:},dt_bin,'g'));

    sd(i,:) = downsample(s(i,:),interval);

    thr=mean(s(i,:) + thr_SD*std(s(i,:)));
    bin(i,find(sd(i,:)>thr))=1;
    bin(i,find(sd(i,:)<=thr))=0;

end

