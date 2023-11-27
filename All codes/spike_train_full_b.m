function [full,bin] = spike_train_full_b(spk_times,start,stop,bin_size_s,thr_SD)

%bin_size_s,startT,stopT have to be in seconds

N=length(spk_times);
startT=start*1000; % start in ms
stopT=stop*1000; % stop in ms
bin_size=bin_size_s*1000; % bin size in ms
trialTime=stopT-startT;
nBins=floor(trialTime/bin_size); 


for iCell=1:length(spk_times)
    normSpk=spk_times{iCell}.*1000; % sec to ms
    normSpk=normSpk-startT; % normalise to start of trial
    normSpk=normSpk(normSpk<=stopT);
    edgesT=linspace(0, stopT-startT, nBins+1);
    binnedSpks=histcounts(normSpk,edgesT);
    binnedSpikes(iCell,1)={binnedSpks};
end

for i=1:N
    full(i,:)=binnedSpikes{i,:};
    thr=mean(full(i,:) + thr_SD*std(full(i,:)));
    bin(i,:)=full(i,:);
    bin(i,find(full(i,:)>thr))=1;
    bin(i,find(full(i,:)<=thr))=0;
end

