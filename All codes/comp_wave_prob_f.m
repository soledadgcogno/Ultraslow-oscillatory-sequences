function [probs_wave,probs_wave_amp] = comp_wave_prob_f(spikes_w,dt,clus)

if dt==0
    [~,sorting,~]=get_sorting(spikes_w);
else
    [~,sorting,~]=get_sorting_smoothed(spikes_w,dt);
end

mat=spikes_w(sorting,:);

if dt==0
    FRp = spikes_downsample(mat,clus,1);
else
    FRp = spikes_downsample(mat,clus,dt);
end

[~,signal]=max(FRp); %Calculate signal by keeping the cluster with largest activity
signal_s(1)=signal(1);
signal_s(2:length(find(diff(signal)~=0))+1) = signal(find(diff(signal)~=0) + 1);

%Checks that the wave has an upwards direction
aux=diff(signal_s);
pos=length(find(aux>6));
neg=length(find(aux<-6));
if pos<neg
else
    cells_per_ens=floor(size(spikes_w,1)/clus);
    sorting_s=sorting(1:clus*cells_per_ens);
    sorting_sf=flip(sorting_s);
    %     sorting_sf=flip(sorting);    
    if dt==0
        FRpi = spikes_downsample(spikes_w(sorting_sf,:),clus,1);
    else
        FRpi = spikes_downsample(spikes_w(sorting_sf,:),clus,dt);
    end        
    clear signal signal_s
    [~,signal]=max(FRpi);
    signal_s(1)=signal(1);
    signal_s(2:length(find(diff(signal)~=0))+1) = signal(find(diff(signal)~=0) + 1);
end

%Calculates all possible wave lengths
[A]=nchoosek(1:clus,2); %All possible transitions found in "clus" clusters
long_s_ten_clus=zeros(1,clus);
long_s_ten_clus_amp=zeros(1,clus);
for i=1:size(A,1)
    
    ini = find(ismember(signal_s,A(i,1))==1);
    en = find(ismember(signal_s,A(i,2))==1);
    in = sort([ini,en],'ascend');
    
    for j=1:size(in,2)-1
        temp=diff(signal_s(in(j):in(j+1)));
        d_temp=floor(sum(temp>0)./length(temp)); %% So all are positive
        if d_temp==1
            long_s_ten_clus(1,length(signal_s(in(j):in(j+1))))=long_s_ten_clus(1,length(signal_s(in(j):in(j+1)))) + 1;
            long_s_ten_clus_amp(1,1+signal_s(in(j+1))-signal_s(in(j)))=long_s_ten_clus_amp(1,1+signal_s(in(j+1))-signal_s(in(j))) + 1;

        end
        clear temp d_temp
    end
end

probs_wave=long_s_ten_clus(1,:)./sum(long_s_ten_clus(1,:));
probs_wave_amp=long_s_ten_clus_amp(1,:)./sum(long_s_ten_clus_amp(1,:));

end

                                    

    