%% Regular wave score
clear all
close all

rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

clusters=10;
mice_number=12;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];

index_c=0;
index_t=0;
window=50;
tic
fs_120=7.73;
for mo=1:mice_number
    
    if mice(mo,1)~= 'L'
        mouse=mice(mo,:);
    else
        if mice(mo,2)=='0'
            mouse=[mice(mo,1),mice(mo,3:5)];
            mouse_l=mouse(2);
            mouse_a=mouse(4);
        else
            mouse=mice(mo,:);
            mouse_l=mouse(2:3);
            mouse_a=mouse(5);
        end
    end
    
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    
    quality=NaN(dates.daysnum,4);
    max_timelag=NaN(dates.daysnum,4);
    rmse=NaN(dates.daysnum,4);
    
    for day=1:dates.daysnum
        disp(day)
        for s=1:dates.sesnum(day)
            munit=dates.ses{day}(s);
            file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
            
            if exist(file_name)==2
                load(file_name,'-mat');
                spikes=full(spikes_d_s);
                [N,T]=size(spikes);
                
                if T>7000 && N>150
                    [coeff,score,latent]=pca(spikes');
                    %                     phase=atan2(latent(2)*score(:,2),latent(1)*score(:,1));
                    phase=atan2(score(:,2),score(:,1));
                    signal=sin(phase);
                    
                    if T>16384*1.1
                        [Powerspec,Powerfreq] = doPwelch(signal,fs_120,16384);
                        
                        flag=1;
                        %                         fourier_c=fourier;
                        fourier=Powerspec(2:140); %3:15
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
                        
                        imax=find(m1>m(i),1,'first');
                        tail=mean(Powerspec(m1(imax)+1:m1(imax)+1+window*2));
                        
                        if  v(i)<9*tail ||v(i)<-9*v1(imin2)% ||
                            flag=0;
                        end
                        
                    else
                        [Powerspec,Powerfreq] = doPwelch(signal,fs_120,8192);
                        flag=1;
                        %                         fourier_c=fourier;
                        fourier=Powerspec(2:70); %3:15
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
                        
                        imax=find(m1>m(i),1,'first');
                        tail=mean(Powerspec(m1(imax)+1:m1(imax)+1+window));
                        
                        if  v(i)<9*tail ||v(i)<-9*v1(imin2)% ||
                            flag=0;
                        end
                        
                    end
                       
                    disp(flag)
                    if flag==0
                        quality(day,s)=0;
                        fft_peak(day,s)=0;
                        dt(day,s)=0;
                        fraction(day,s)=0;
                    elseif flag==1
                        
                        [oscil,~,~] = comp_timelag_prob_3D(spikes);
                        WS = sum(oscil)./(length(oscil)-3);
                        fraction(day,s)=sum(oscil)/length(oscil);
                        if WS>=1
                            quality(day,s)=1;               
                            
                            [Powerspec,Powerfreq] = doPwelch(signal,fs_120,8192);
                            P1_max=max(Powerspec(3:end));
                            ind=find(Powerspec==P1_max);
                            period=1/Powerfreq(ind);
                            
                            dt_aux=period*fs_120/10; %10 is such that one sequence evolved in 10 time bins
                            fft_peak(day,s)=Powerfreq(ind);
                            dt(day,s)=dt_aux;
                        else
                            quality(day,s)=0;
                            fft_peak(day,s)=0;
                            dt(day,s)=0;
                        end
                    end
          
                else
                    quality(day,s)=NaN;
                    fft_peak(day,s)=NaN;
                    dt(day,s)=NaN;
                    fraction(day,s)=NaN;

                end
            end
            clear signal Powerspec Powerfreq spikes N T oscil score phase coeff spikes_d_s f gof ratio peaks cells_d fourierY L P2 P1 f P1_max ind dt_aux period
        end
    end    
    WS_stat.fraction=fraction;
    WS_stat.WS=quality;
    WS_stat.FFT=fft_peak;
    WS_stat.dt=dt;       
%     save([save_data_path,'WS_Osc_15_sf7p73II_',mice(mo,:)],'WS_stat');    
    clear optimal_dt wave_score_ent dates mouse thr Entr_up WS_stat WS
    clear peak_f osc H max_timelag rmse quality fft_peak dt fraction 
    
end
toc

