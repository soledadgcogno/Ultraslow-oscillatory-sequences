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
    ws_prob=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_sf7p73_',mice(m,:),'.mat']]);
    ws_prob_sig=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_sf7p73_',mice(m,:),'with_significance.mat']]);
    ws_ent=load([save_data_path ['WS_Entropy_dt66_',mice(m,:),'.mat']]);
    ws_prob_old=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_',mice(m,:),'.mat']]);
    ws_prob_sig_old=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_',mice(m,:),'with_significance.mat']]);
    
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
                    big_table(count,6)=WS_stat.WS(day,s); %WS - Osc
                    big_table(count,11)=ws_ent.WS_stat.wave_score_ent(day,s); %WS - Ent
                    big_table(count,12)=size(spk.spikes_d_s,1); %N
                    big_table(count,13)=WS_stat.fraction(day,s); %fraction - Osc
                    big_table(count,14)=ws_prob.WS_stat.seq_score_prob(day,s); %WS - Prob
                    big_table(count,15)=ws_prob_sig.WS_stat.seq_score_prob_sig(day,s); %WS - Prob
                    big_table(count,16)=ws_prob_old.WS_stat.wave_score_prob(day,s); %WS - Prob
                    big_table(count,17)=ws_prob_sig_old.WS_stat.wave_score_prob_sig(day,s); %WS - Prob


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
    clear WS_stat ws_ent ws_prob ws_prob_sig ws_prob_old ws_prob_sig_old
end


%% Waveness, wave score and fraction - Data

%------------------------- V1
V1_ses=find(big_table(:,10)<0); %Sessions inthe in MEC
adults=find(big_table(:,3)>15); %Sessions inthe of adults
sessions=intersect(V1_ses,adults);
count=0;

ws=big_table(sessions,6);
index_nan=isnan(ws);
sessions(index_nan)=[];
session_number=big_table(sessions,4);
days=big_table(sessions,3);
repeated_days=find(diff(days)==0);
clear index_to_delete

index_to_delete=[];
count=0;
for w=1:length(repeated_days)
    N1=big_table(sessions(repeated_days(w)),12);
    N2=big_table(sessions(repeated_days(w)+1),12);
    
    if (N1>N2)
        count=count+1;
        index_to_delete(count)=repeated_days(w)+1;
    elseif  (N1<N2)
        count=count+1;
        index_to_delete(count)=repeated_days(w);    
    else
        print('Other case');
    end   
end

sessions(index_to_delete)=[];
big_table_2=big_table(sessions,:);
V1_sessions=sessions;
for w=1:length(V1_sessions)
    
    row_w=V1_sessions(w);
    disp(w)
    
    count=count+1;
    mouse=['L',num2str(big_table(row_w,1)),'M',num2str(big_table(row_w,2))];
    day=big_table(row_w,3);
    s=big_table(row_w,4);
    munit=big_table(row_w,5);
    ws_binary_V1(w)=big_table(row_w,6);
    ws_entropy_V1(w)=big_table(row_w,11);
    ws_fraction_V1(w)=big_table(row_w,13);
    ws_probability_V1(w)=big_table(row_w,14);

%     ML_pos(w)=big_table(row_w,10);
    
end


%------------------------- PaS

PaS_ses1=find(big_table(:,10)<=3); %Sessions inthe in MEC
PaS_ses2=find(big_table(:,10)>0); %Sessions inthe in MEC
PaS_ses=intersect(PaS_ses1,PaS_ses2);
adults=find(big_table(:,3)>15); %Sessions inthe of adults
sessions=intersect(PaS_ses,adults);
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
    N1=big_table(sessions(repeated_days(w)),12);
    N2=big_table(sessions(repeated_days(w)+1),12);
    
    if (N1>N2)
        count=count+1;
        index_to_delete(count)=repeated_days(w)+1;
    elseif  (N1<N2)
        count=count+1;
        index_to_delete(count)=repeated_days(w);    
    else
        print('Other case');
    end   
end

sessions(index_to_delete)=[];

big_table_2=big_table(sessions,:);

PaS_sessions=sessions;

for w=1:length(PaS_sessions)
    
    row_w=PaS_sessions(w);
    disp(w)    
    count=count+1;
    mouse=['L',num2str(big_table(row_w,1)),'M',num2str(big_table(row_w,2))];
    day=big_table(row_w,3);
    s=big_table(row_w,4);
    munit=big_table(row_w,5);
    ws_binary_PaS(w)=big_table(row_w,6);
    ws_entropy_PaS(w)=big_table(row_w,11);
    ws_fraction_PaS(w)=big_table(row_w,13);
    ws_probability_PaS(w)=big_table(row_w,14);
 
%     ML_pos(w)=big_table(row_w,10);    
end


%------------------------- MEC


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

for w=1:length(mec_sessions)
    
    row_w=mec_sessions(w);
    disp(w)
    
    count=count+1;
    mouse=['L',num2str(big_table(row_w,1)),'M',num2str(big_table(row_w,2))];
    day=big_table(row_w,3);
    s=big_table(row_w,4);
    munit=big_table(row_w,5);
    ws_binary_mec(w)=big_table(row_w,6);
    ws_entropy_mec(w)=big_table(row_w,11);
    ws_fraction_mec(w)=big_table(row_w,13);
    ws_probability_mec(w)=big_table(row_w,14);


%     ML_pos(w)=big_table(row_w,10);
    
end

waveness_total=[ws_entropy_mec,ws_entropy_PaS,ws_entropy_V1]';
wave_score_total=[ws_binary_mec,ws_binary_PaS,ws_binary_V1]';
fraction_total=[ws_fraction_mec,ws_fraction_PaS,ws_fraction_V1]';
waveness_prob_total=[ws_probability_mec,ws_probability_PaS,ws_probability_V1]';

%% Concatenated sessions
sessions_to_keep=[mec_sessions;PaS_sessions;V1_sessions];

table=big_table(sessions_to_keep,[1,2,3,4,5,13,6,14,15]);