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
                    %                     figure
                    %                     plot(Powerspec)
                    %                     hold on
                    %                     refline(0,9*tail)
                    %                     c=refline(0,-9*v1(imin2));
                    %                     c.Color = 'r';
                    %                     title([day flag])
                    
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
                            
                            %                             Y = fft(signal);
                            %                             L=length(signal);
                            %                             P2 = abs(Y/L);
                            %                             P1 = P2(1:L/2+1);
                            %                             P1(2:end-1) = 2*P1(2:end-1);
                            %                             f = 8*(0:(L/2))/L;
                            %                             P1_max=max(P1(3:end));
                            %                             ind=find(P1==P1_max);
                            %                             period=1/f(ind);
                            
                            
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
                    
                    %
                    %                     disp(flag)
                    %                     if flag==0
                    %                         quality(day,s)=0;
                    %
                    %                     elseif flag==1
                    %                         [oscil,~,~] = comp_timelag_prob_3D(spikes);
                    %                         WS = sum(oscil)./(length(oscil)-3);
                    %
                    %                         if WS>=1
                    %                             quality(day,s)=1;
                    %                             max_timelag(day,s)=NaN; %Improve
                    %                         else
                    %                             quality(day,s)=0;
                    %                             max_timelag(day,s)=NaN; %Improve
                    %                         end
                    %                     end
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
    %     WS_stat.max_timelag=max_timelag;
    %     WS_stat.rmse=rmse;
    WS_stat.FFT=fft_peak;
    WS_stat.dt=dt;
    
    
    save([save_data_path,'WS_Osc_15_sf7p73II_',mice(mo,:)],'WS_stat');
    
    clear optimal_dt wave_score_ent dates mouse thr Entr_up WS_stat WS
    clear peak_f osc H max_timelag rmse quality fft_peak dt fraction 
    
end
toc

%% Wavescore after shuffling spike trains


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
N_sh=200;

tic
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
    
    quality=NaN(dates.daysnum,N_sh,4);
    mean_quality=NaN(dates.daysnum,4);
    std_quality=NaN(dates.daysnum,4);
    %     max_timelag=NaN(dates.daysnum,4);
    %     rmse=NaN(dates.daysnum,4);
    
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
                    
                    for sh=1:N_sh
                        disp(sh)
                        mat_sh=shuffle(spikes')';
                        [coeff,score,latent]=pca(mat_sh');
                        phase=atan2(score(:,2),score(:,1));
                        signal=sin(phase);
                        
                        
                        if T>16384*1.1
                            [Powerspec,Powerfreq] = doPwelch(signal,8,16384);
                            
                            flag=1;
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
                            [Powerspec,Powerfreq] = doPwelch(signal,8,8192);
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
                        
                        if flag==0
                            quality(day,sh,s)=0;
                        else
                            [oscil,~,~] = comp_timelag_prob_3D(spikes);
                            WS = sum(oscil)./(length(oscil)-3);
                            
                            if WS>=1
                                quality(day,sh,s)=1;
                            else
                                quality(day,sh,s)=0;
                            end
                        end
                        clear signal Powerspec Powerfreq N oscil score phase coeff  f gof ratio peaks cells_d mat_sh latent
                    end
                    clear spikes spikes_d_s T
                end
            end
            mean_quality(day,s)=mean(quality(day,:,s));
            std_quality(day,s)=std(quality(day,:,s));
            
        end
    end
    
    WS_stat_shuffle.WS=quality;
    WS_stat_shuffle.WS_mean=mean_quality;
    WS_stat_shuffle.WS_std=std_quality;
    
    %     WS_stat.max_timelag=max_timelag;
    %     WS_stat.rmse=rmse;
    
    save([save_data_path,'WS_Osc_14_shuffle_',mice(mo,:)],'WS_stat_shuffle');
    
    clear optimal_dt wave_score_ent dates mouse thr Entr_up WS_stat WS
    clear peak_f osc H max_timelag rmse quality std_quality mean_quality
    
end
toc

%% Figure of wave score in real session VS shuffle
% clear all
% close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';


mice_number=12;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];

ml=[2.64,3.7,2.7,2.64,3.2,3.3,2.52,3.3,3.5,-10,-10,-10]; %ML coordinates

count=0;
for m=1:mice_number
    
    if mice(m,1)~= 'L'
        mouse=mice(m,:);
    else
        if mice(m,2)=='0'
            mouse=[mice(m,1),mice(m,3:5)];
            mouse_l=mouse(2);
            mouse_a=mouse(4);
        else
            mouse=mice(m,:);
            mouse_l=mouse(2:3);
            mouse_a=mouse(5);
        end
    end
    
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    load([save_data_path ['WS_Osc_14_',mice(m,:),'.mat']]);
    load([save_data_path ['WS_Osc_14_shuffle_',mice(m,:),'.mat']]);
    
    
    for day=1:dates.daysnum
        for s=1:dates.sesnum(day)
            disp(s)
            munit=dates.ses{day}(s);
            
            file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
            
            if (exist(file_name_spk) == 2)
                disp(day)
                
                if s==1 %|| s==2
                    count=count+1;
                    
                    if isfield(dates,'actual_day')  == 1
                        big_table(count,3)=dates.actual_day(day);   %Day # on the wheel
                    else
                        big_table(count,3)=day; %Day # on the wheel
                    end
                    
                    if mouse(1) == '9' || mouse(1) == '6'
                        big_table(count,1)=1;%str2num(mice(m,2:3));
                        big_table(count,10)=-10; %ML position
                        if mouse=='92227'
                            big_table(count,2)=1;
                        elseif mouse == '92229'
                            big_table(count,2)=2;
                        elseif mouse == '60961'
                            big_table(count,2)=3;
                        end
                    else
                        big_table(count,1)=str2num(mice(m,2:3)); %Litter
                        big_table(count,2)=str2num(mice(m,end)); %Mouse number
                        big_table(count,10)=ml(m); %ML position
                    end
                    
                    big_table(count,4)=s; %Session
                    big_table(count,5)=munit; %munit
                    big_table(count,6)=WS_stat.WS(day,s); %WS
                    big_table(count,11)=WS_stat_shuffle.WS_mean(day,s); %Mean shuffle
                    big_table(count,12)=WS_stat_shuffle.WS_std(day,s); %Mean shuffle
                    big_table(count,13)=length(find(WS_stat_shuffle.WS(day,:,s)>0))/N_sh; %Mean shuffle
                    
                    
                    threshold_kl=1;
                    if ( big_table(count,6)>=threshold_kl )
                        big_table(count,7)=1;
                    else
                        big_table(count,7)=0;
                    end
                    
                    if (isnan(big_table(count,6)))
                        big_table(count,7)=NaN;
                        big_table(count,8)=NaN; %dt
                        big_table(count,9)=NaN; %fft
                    else
                        big_table(count,8)=WS_stat.dt(day,s); %dt
                        big_table(count,9)=WS_stat.FFT(day,s); %fft
                    end
                    
                end
                
            end
        end
        
    end
    clear WS_stat WS_stat_shuffle
end

adults=find(big_table(:,3)>15);%Sessions inthe table with waves
waves_ses=find(big_table(:,6)==1); %Sessions in the table with waves
nowaves_ses=find(big_table(:,6)==0); %Sessions in the table with NO waves
waves=intersect(waves_ses,adults);
no_waves=intersect(nowaves_ses,adults);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures
figure
hold on
for i=1:length(waves)
    errorbar([0,1],[big_table(waves(i),6),big_table(waves(i),11)],[0,big_table(waves(i),12)],'-*','linewidth',2);
end
axis([-0.5 1.5 -0.2 1.2])
% legend boxoff
ylabel('Wave score');
xticks([0 1]);
xticklabels({'Session','Shuffle'});
set(gca,'fontsize',16)
yticks([0 1])

figure
hold on
for i=1:length(no_waves)
    errorbar([0,1],[big_table(no_waves(i),6),big_table(no_waves(i),11)],[0,big_table(no_waves(i),12)],'-*','linewidth',2);
end
axis([-0.5 1.5 -0.2 1.2])
% legend boxoff
ylabel('Wave score');
xticks([0 1]);
xticklabels({'Session','Shuffle'});
set(gca,'fontsize',16)
yticks([0 1])


figure
hold on
for i=1:length(waves)
    plot(i*[1,1],[big_table(waves(i),6),big_table(waves(i),11)],'--','color',[115,115,115]./255,'linewidth',2);
    % e.Color= [115,115,115]./255;
end
for i=1:length(waves)
    e=errorbar(i*[1,1],[big_table(waves(i),6),big_table(waves(i),11)],[0,big_table(waves(i),12)],'ko','linewidth',2);
    e.Color= 'k';
end
hold on
for i=1:length(waves)
    scatter(i,[big_table(waves(i),6)],80,'o','Markerfacecolor','k','Markeredgecolor','k');
end
hold on
for i=1:length(waves)
    scatter(i,[big_table(waves(i),11)],80,'o','Markerfacecolor','w','Markeredgecolor','k');
end
axis([0 14 -0.3 1.3]);
ylabel('Wave score');
yticks([0 1]);
xlabel('Session #');
set(gca,'fontsize',18);
xticks([1 7 13])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Histograms
vec_wave=[];
vec_nowave=[];

for i=1:length(waves)
    vec_wave(i)=big_table(waves(i),13);
end

for i=1:length(no_waves)
    vec_nowave(i)=big_table(no_waves(i),13);
end


figure
H=histogram(vec_wave,[0:0.005:0.15]);
figure
bar(H.BinEdges(1:end-1),H.Values,'k');
axis([-0.005 0.04 0 12])
ylabel('Counts_W_a_v_e_ _s_e_s_s_i_o_n_s');
xlabel('Fraction of shuffles for which Wave score = 1');
set(gca,'fontsize',18)
xticks([0 0.02 0.04])
yticks([0 11]);

figure
H=histogram(vec_nowave,[0:0.005:0.15]);
figure
bar(H.BinEdges(1:end-1),H.Values,'k');
axis([-0.005 0.04 0 52])
ylabel('Counts_N_o_ _W_a_v_e_ _s_e_s_s_i_o_n_s');
xlabel('Fraction of shuffles for which Wave score = 1');
set(gca,'fontsize',18)
xticks([0 0.02 0.04])
yticks([0 50]);

figure
H=histogram(vec_wave,[0:0.005:0.15],'Normalization','Probability');
figure
bar(H.BinEdges(1:end-1),H.Values,'k');
axis([-0.005 0.04 0 1])
ylabel('Fraction of Wave sessions');
xlabel('Fraction of shuffles for which Wave score = 1');
set(gca,'fontsize',18)
xticks([0 0.02 0.04])
yticks([0 1]);

figure
H=histogram(vec_nowave,[0:0.005:0.15],'Normalization','Probability');
figure
bar(H.BinEdges(1:end-1),H.Values,'k');
axis([-0.005 0.04 0 1])
ylabel('Fraction of No Wave sessions');
xlabel('Fraction of shuffles for which Wave score = 1');
set(gca,'fontsize',18)
xticks([0 0.02 0.04])
yticks([0 1]);

%% Wavescore after subsampling cells

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
N_sh=200;

tic
for mo=2:mice_number
    
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
    
    quality0p9=NaN(dates.daysnum,N_sh,4);
    quality0p7=NaN(dates.daysnum,N_sh,4);
    quality0p5=NaN(dates.daysnum,N_sh,4);
    fraction_quality0p5=NaN(dates.daysnum,4);
    fraction_quality0p7=NaN(dates.daysnum,4);
    fraction_quality0p9=NaN(dates.daysnum,4);
    
    for day=15:dates.daysnum
        disp(day)
        for s=1:dates.sesnum(day)
            munit=dates.ses{day}(s);
            file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
            
            if exist(file_name)==2
                load(file_name,'-mat');
                spikes=full(spikes_d_s);
                [N,T]=size(spikes);
                
                if T>7000 && N>150
                    
                    count=0;
                    for prctage = [0.5,0.7,0.9]
                        N_new= floor(N*prctage);
                        for sh=1:N_sh
                            disp(sh)
                            cells=randperm(N,N_new);
                            mat_sh=spikes(cells,:);
                            [coeff,score,latent]=pca(mat_sh');
                            phase=atan2(score(:,2),score(:,1));
                            signal=sin(phase);
                            
                            
                            if T>16384*1.1
                                [Powerspec,Powerfreq] = doPwelch(signal,8,16384);
                                
                                flag=1;
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
                                [Powerspec,Powerfreq] = doPwelch(signal,8,8192);
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
                                                   
                            if (prctage == 0.5)
                                if flag==0
                                    quality0p5(day,sh,s)=0;
                                else
                                    [oscil,~,~] = comp_timelag_prob_3D(spikes);
                                    WS = sum(oscil)./(length(oscil)-3);
                                    
                                    if WS>=1
                                        quality0p5(day,sh,s)=1;
                                    else
                                        quality0p5(day,sh,s)=0;
                                    end
                                end
                            elseif (prctage == 0.7)
                                if flag==0
                                    quality0p7(day,sh,s)=0;
                                else
                                    [oscil,~,~] = comp_timelag_prob_3D(spikes);
                                    WS = sum(oscil)./(length(oscil)-3);
                                    
                                    if WS>=1
                                        quality0p7(day,sh,s)=1;
                                    else
                                        quality0p7(day,sh,s)=0;
                                    end
                                end
                                
                            elseif (prctage == 0.9)
                                if flag==0
                                    quality0p9(day,sh,s)=0;
                                else
                                    [oscil,~,~] = comp_timelag_prob_3D(spikes);
                                    WS = sum(oscil)./(length(oscil)-3);
                                    
                                    if WS>=1
                                        quality0p9(day,sh,s)=1;
                                    else
                                        quality0p9(day,sh,s)=0;
                                    end
                                end
                            end
                            
                            
                            clear signal Powerspec Powerfreq oscil score phase coeff  f gof ratio peaks cells_d mat_sh latent cells v1 v m1 m fourier
                        end                        
                    end
                    
                    clear spikes spikes_d_s T N
                end
            end
            fraction_quality0p5(day,s)=length(find(quality0p5(day,:,s)>0))/N_sh;
            fraction_quality0p7(day,s)=length(find(quality0p7(day,:,s)>0))/N_sh;
            fraction_quality0p9(day,s)=length(find(quality0p9(day,:,s)>0))/N_sh;

        end
    end
    
    WS_stat_dropcells.score0p5=quality0p5;
    WS_stat_dropcells.score0p7=quality0p7;
    WS_stat_dropcells.score0p9=quality0p9;
    WS_stat_dropcells.fraction0p5=fraction_quality0p5;
    WS_stat_dropcells.fraction0p7=fraction_quality0p7;
    WS_stat_dropcells.fraction0p9=fraction_quality0p9;

    
    save([save_data_path,'WS_Osc_14_dropcells_',mice(mo,:)],'WS_stat_dropcells');
    
    clear optimal_dt wave_score_ent dates mouse thr Entr_up WS_stat WS
    clear peak_f osc H max_timelag rmse quality std_quality mean_quality quality0p9 quality0p7 quality0p5 fraction_quality0p5 fraction_quality0p7 ...
        fraction_quality0p9
    
end
toc
