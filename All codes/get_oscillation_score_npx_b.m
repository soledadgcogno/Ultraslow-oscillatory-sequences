function [fraction,dt] = get_oscillation_score_npx_b(spk_chunk,spk_all,bin_size)

fs_120=1/bin_size;
quality=NaN;
fft_peak=NaN;
dt=NaN;
fraction=NaN;

[coeff,score,~]=pca(spk_chunk');
[~,b]=get_sorting_npx(spk_chunk);
spy(spk_all(b,:))
pbaspect([25 2.5 1]);

score1_full=spk_all'*coeff(:,1);
score2_full=spk_all'*coeff(:,2);

phase=atan2(score(:,2),score(:,1));
phase_f=atan2(score2_full,score1_full);

phase=phase_f;
signal=sin(phase);

[Powerspec,Powerfreq] = doPwelch(signal,fs_120,2048);

flag=1;
upper_bound=find(Powerfreq>0.1,1);
fourier=Powerspec(2:upper_bound); %3:15
[m,v]=find_peaks_smooth(fourier',1,0);
sel=[find(m==1) find(m==numel(fourier))];
m(sel)=[]; v(sel)=[];
[~,i]=max(v);
[m1,v1]=find_peaks_smooth(-fourier',1,0);
if numel(m)==0
    flag=0;
    return;
end

imin=find(m1<m(i));
[~,imin2]=min(-v1(imin));

window=find(Powerfreq>0.02,1);
imax=find(m1>m(i),1,'first');
tail=mean(Powerspec(m1(imax)+1:m1(imax)+1+window*2));

if  v(i)<4*tail ||v(i)<-4*v1(imin2)% || %Used to be 9
    flag=0;
end

% disp(flag)
if flag==0
    quality=0;
    fft_peak=0;
    dt=0;
    fraction=0;
elseif flag==1
    [oscil,~,~] = comp_timelag_prob_3D_npx_b(spk_chunk,spk_all);
    WS = sum(oscil)./(length(oscil)-3);
    fraction=sum(oscil)/length(oscil);
    if WS>=1
        quality=1;
        [Powerspec,Powerfreq] = doPwelch(signal,fs_120,1024*2);
        P1_max=max(Powerspec(3:end));
        ind=find(Powerspec==P1_max);
        period=1/Powerfreq(ind);
        dt_aux=period*fs_120/10; %10 is such that one sequence evolved in 10 time bins
        fft_peak=Powerfreq(ind);
        dt=dt_aux;
    else
        quality=0;
        fft_peak=0;
        dt=0;
    end
end
