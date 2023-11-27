clear all
close all

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
    load([save_data_path ['WS_Osc_15_sf7p73II_',mice(m,:),'.mat']]);
%     ws_ent=load([save_data_path ['WS_Entropy_',mice(m,:),'.mat']]);
    ws_ent=load([save_data_path ['WS_Entropy_dt66_',mice(m,:),'.mat']]);
    ws_prob=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_sf7p73_',mice(m,:),'.mat']]);
    ws_prob_sig=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_sf7p73_',mice(m,:),'with_significance.mat']]);

    
    for day=1:dates.daysnum
        for s=1:dates.sesnum(day)
            disp(s)
            munit=dates.ses{day}(s);
            
            file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
            
            if (exist(file_name_spk) == 2)
                disp(day)
                spk=load(file_name_spk);
                
                if s<100
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
                    big_table(count,6)=WS_stat.WS(day,s); %oscillation score binary
                    big_table(count,11)=ws_ent.WS_stat.wave_score_ent(day,s); % Entropy - super old
                    big_table(count,12)=size(spk.spikes_d_s,1); %N
                    big_table(count,13)=WS_stat.fraction(day,s); %fraction - Osc score
                    big_table(count,14)=ws_prob.WS_stat.seq_score_prob(day,s); %Sequence score
                    big_table(count,15)=ws_prob.WS_stat.seq_score_prob(day,s); %Sequence score
                    big_table(count,16)=ws_prob_sig.WS_stat.seq_score_prob_sig(day,s); %Sequence score with significance



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
               clear spk 
            end
        end
        
    end
    clear WS_stat ws_ent ws_prob
end


%% Wave sessions


figpath='C:\Users\xscogno\MATLAB\Flavio2\Waves\Semi final scripts\Raster Plots all MEC sessions\';
MEC_ses=find(big_table(:,10)>3); %Sessions inthe in MEC
adults=find(big_table(:,3)>15); %Sessions inthe of adults
sessions=intersect(MEC_ses,adults);
count=0;

ws=big_table(sessions,6);
index_nan=isnan(ws);
sessions(index_nan)=[];
session_number=big_table(sessions,4);
days=big_table(sessions,3);
repeated_days=find(diff(days)==0);
clear index_to_delete

count=0;
for w=1:length(repeated_days)
    ws1=big_table(sessions(repeated_days(w)),6);
    ws2=big_table(sessions(repeated_days(w)+1),6);
    
    if (ws1==1 && ws2==0)
        count=count+1;
        index_to_delete(count)=repeated_days(w)+1;
    elseif  (ws1==0 && ws2==1)
        count=count+1;
        index_to_delete(count)=repeated_days(w);
    elseif (ws1==0 && ws2==0)
        count=count+1;
        index_to_delete(count)=repeated_days(w);
    else
        print('Other case');
    end   
end

sessions(index_to_delete)=[];
big_table_2=big_table(sessions,:);
mec_sessions=sessions;
waves=mec_sessions(find(big_table(mec_sessions,6)==1));

%% Variance


% Phase and ring
figure
countf=0;
for w=1:length(waves)

    row_w=waves(w);
    disp(w)

    count=count+1;
    mouse=['L',num2str(big_table(row_w,1)),'M',num2str(big_table(row_w,2))];
    day=big_table(row_w,3);
    s=big_table(row_w,4);
    munit=big_table(row_w,5);

    clus=10;
    disc_phase=10;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

    load(file_name_spk,'-mat'); %Spike times
    spikes_d=full(spikes_d_s);

    [coeff,score,latent,tsquared,explained]=pca(spikes_d');

    variance_explained(w)=sum(explained(1:2));
    ration_eigen12(w)=explained(2)/explained(1);
    ration_eigen13(w)=explained(3)/explained(1);
    var(w,:)=explained(1:3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Condition on having waves
    sf=7.73;
    N=size(spikes_d,1);
    num_clus_discr=10;
    make_fig=0;
    dt=floor(big_table(row_w,8));

    %     FRp = spikes_downsample(spikes_d,N,1);
    FRp = spikes_d;

    if w==7
        for i=1:N
            FR(i,:)=(full(fire_rate(FRp(i,:),2*floor(dt/1),'g')));
        end
    else
        for i=1:N
            FR(i,:)=(full(fire_rate(FRp(i,:),4*floor(dt/1),'g')));
        end
    end

    clear coeff score latent tsquared explained
    [coeff,score,latent,tsquared,explained]=pca(FR');

    variance_explainedFR(w)=sum(explained(1:2));
    ration_eigen12FR(w)=explained(2)/explained(1);
    ration_eigen13FR(w)=explained(3)/explained(1);
    varFR(w,:)=explained(1:3);


    FRp = spikes_downsample(spikes_d,N,4);
    clear coeff score latent tsquared explained
    [coeff,score,latent,tsquared,explained]=pca(FRp');

    variance_explainedDown(w)=sum(explained(1:2));
    ration_eigen12Down(w)=explained(2)/explained(1);
    ration_eigen13Down(w)=explained(3)/explained(1);
    varDOWN(w,:)=explained(1:3);


    clear X Y spikes spikes_d_s score mouse i FRp FR coeff cells_d T window score FR FRp epoch ICI_ses table_u spikes spikes_d ICI_ses_t
    clear time_points  spikes_r FR FRp

end

%%
figure
plot(variance_explained,'-*')
xlabel('Session #');
ylabel('Variance explained by first 2 PCs %');
title('Matrix of calcium activity')

figure
plot(variance_explainedFR,'-*')
xlabel('Session #');
ylabel('Variance explained by first 2 PCs %');
title('Smoothed matrix of calcium activity')

figure
plot(variance_explainedDown,'-*')
xlabel('Session #');
ylabel('Variance explained by first 2 PCs %');
title('Downsampled matrix of calcium activity')

figure
boxplot(var)
ylabel('Variance explained by PC')
xlabel('PC')
title('Matrix of calcium activity')

[p,h,stat]=ranksum(var(:,1),var(:,2));
[p,h,stat]=ranksum(var(:,1),var(:,3));
[p,h,stat]=ranksum(var(:,2),var(:,3));


figure
boxplot(varFR)
ylabel('Variance explained by PC')
xlabel('PC')
title('Smoothed matrix of calcium activity')
[p,h,stat]=ranksum(varFR(:,1),varFR(:,2));
[p,h,stat]=ranksum(varFR(:,1),varFR(:,3));
[p,h,stat]=ranksum(varFR(:,2),varFR(:,3));

figure
boxplot(varDOWN)
[p,h,stat]=ranksum(varDOWN(:,1),varDOWN(:,2));
[p,h,stat]=ranksum(varDOWN(:,1),varDOWN(:,3));
[p,h,stat]=ranksum(varDOWN(:,2),varDOWN(:,3));

