function [oscil,peak_f,max_timelag] = comp_timelag_prob_3D_npx_b(spk_chunk,spikes)

% Lagged correlation
maxlag=250;%250; % 480 for calcium data
downsampling_factor=4;
N=size(spikes,1);
FRp = spikes_downsample(spikes,N,downsampling_factor);
FR=FRp;
% FR=spikes;
sf=8/downsampling_factor;

count=0;
for i=1:N
    for j=i+1:N
        [val,time]=(xcorr(FR(i,:),FR(j,:),maxlag)); %Check whether I need to zscore
        
        [v,in]=max(val);
        if time(in)>=0
            Corr_mat(i,j)=v;
            %             Corr_mat(j,i)=-v;
            Corr_mat(j,i)=v;
            lag_mat(i,j)=time(in);
            lag_mat(j,i)=-time(in);
        else
            %             Corr_mat(i,j)=-v;
            Corr_mat(i,j)=v;
            Corr_mat(j,i)=v;
            lag_mat(i,j)=time(in);
            lag_mat(j,i)=-time(in);
        end
        clear val time
    end
end

% Position and distance on the ring
[coeff]=pca(spikes');
% [coeff]=pca(spk_chunk');
% [coeff]=pca(spikes_w');

angle_pca=atan2(coeff(:,2),coeff(:,1));

for i=1:N
    for j=1:N
        delta_ring_pca(i,j)=(angdiff(angle_pca(i),angle_pca(j)));
    end
end
delta_ring_pca_vec=delta_ring_pca(:);
Corr_vec=Corr_mat(:);
lag_vec=lag_mat(:);

aux=eye(N,N);
aux2=aux(:);
template=find(aux2==1);

delta_ring_pca_vec(template)=[];
% delta_ring_umap_vec(template)=[];
Corr_vec(template)=[];
lag_vec(template)=[];
corrm=nan(30000,11);
dist_ring_edges=-3.14:2*3.14/11:3.14;
for i=1:length(dist_ring_edges)-1
    a=find(delta_ring_pca_vec>dist_ring_edges(i));
    b=find(delta_ring_pca_vec<=dist_ring_edges(i+1));
    c=intersect(a,b);
    corrm(1:length(lag_vec(c)),i)=lag_vec(c);
    clear a b c
end

% edges=[-maxlag:10:maxlag];
vals=[-maxlag:4:maxlag]./sf;
clear y aux
figure
for i=1:length(dist_ring_edges)-1
    figure
    aux=histogram((corrm(:,i)),-maxlag:4:maxlag);%,'Normalization','Probability');
    %     aux=histogram(abs(corrm(:,i)),10:90);%,'Normalization','Probability');
    %     y(i,:)=aux.Values./sum(aux.Values);
    y(i,:)=aux.Values;
    [a,b]=max(y(i,2:end-1));
    max_prob(i)=b;
    prob_wave(i)=nansum(-y(i,:).*log( y(i,:)));
end
% close all

y=y./sum(sum(y));
%
figure
imagesc(flip(y)')
xticks([1 6 11])
xticklabels({'-3.1 - -2.6','-0.3 - 0.3','2.6 - 3.14'})
% yticks([1 13 24])
% yticklabels({60,0,-60})
xlabel('Distance on ring (rad)');
ylabel('Time lag (s)');
co=colorbar;
% colormap puor
colormap gmt_ocean;
% colormap cividis
% colormap gmt_nighttime;
caxis([0 0.001]);
set(gca,'fontsize',16)
xtickangle(45)
axis square


for j=1:size(y,1)
    signal=y(j,2:end-1);
    signal_c=signal;
    signal=[signal,signal,signal];

%     aux=find(zscore(signal_c)<0);
%     sigp=zscore(signal_c);
%     sigp(aux)=0;
%     signal=[sigp,sigp,sigp];

    clear Powerspec  Powerfreq quality
    [Powerspec,Powerfreq] = doPwelch(signal,2,512/2);
    min_freq=3/(length(signal)/2);
    quality=check_peak_quality_3_npx_b(Powerfreq,Powerspec,min_freq);
%     figure; subplot(1,2,1); plot((1:length(signal))./2,signal); subplot(1,2,2); plot(Powerfreq,Powerspec);

    if (quality==1)
        oscil(j)=1;
        [a,b]=findpeaks(Powerspec);
        c=find(max(a));
        peak_f(j)=Powerfreq(b(c));
        [~,a2]=max(signal_c);
        max_timelag(j)=vals(a2);
    else
        oscil(j)=0;
        peak_f(j)=0;
        max_timelag(j)=0;
    end
end

% close all
end



    