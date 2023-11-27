function [peak_freq,peak_psd,goodness]=check_peak_quality_3b_f(fourier,freq)
goodness=1;

fourier_c=fourier;
fourier=fourier(5:50); %3:15
[m,v]=find_peaks_smooth(fourier',1,0);
sel=[find(m==1) find(m==numel(fourier))];
m(sel)=[]; v(sel)=[];
[~,i]=max(v);
[m1,v1]=find_peaks_smooth(-fourier',1,0);
if numel(m)==0
    goodness=0;
    peak_freq=inf;
    peak_psd=inf;
    return;
end

% figure
% plot(fourier(1:40))

imin=find(m1<m(i));

if numel(imin)==0
    goodness=0;
    peak_freq=inf;
    peak_psd=inf;
    return;
end

[~,imin2]=min(-v1(imin));

imax=find(m1>m(i),1,'first');
tail=mean(fourier_c(m1(imax):end));

%if  v(i)<5*tail || v(i)<-1.5*v1(imin2)% ||
if  v(i)<5*tail || v(i)<0.04 || v(i)<-1.5*v1(imin2)% ||    
    goodness=0;
%    peak_freq=inf;
%    peak_psd=inf;
 peak_psd=v(i);%fourier_c(m(i)+2);
    in=fourier_c==peak_psd;
    peak_freq=freq(in);
else    
    %peak_freq=freq(m(i)+2);
    peak_psd=v(i);%fourier_c(m(i)+2);
    in=find(fourier_c==peak_psd);
    peak_freq=freq(in);
end

return;
end