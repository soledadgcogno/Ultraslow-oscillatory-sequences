% Fig6g
% Extended data Fig 10q
% Extended data Fig. 10r
% Extended data Fig. 12f
%% 
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
                    big_table(count,6)=WS_stat.WS(day,s); %Oscillation score - binary
                    big_table(count,11)=ws_ent.WS_stat.wave_score_ent(day,s); %OLD : entropy
                    big_table(count,12)=size(spk.spikes_d_s,1); %N
                    big_table(count,13)=WS_stat.fraction(day,s); %Oscillation score 
                    big_table(count,14)=ws_prob.WS_stat.seq_score_prob(day,s); %Sequence score
                    big_table(count,15)=ws_prob_sig.WS_stat.seq_score_prob_sig(day,s); %Sequence score with significance

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

%% sequence score and oscillatio score - Data

%------------------------- VIS
V1_ses=find(big_table(:,10)<0); 
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
    seq_score_sig_v1(w)=big_table(row_w,15);
%     ML_pos(w)=big_table(row_w,10);    
end


%------------------------- PaS
PaS_ses1=find(big_table(:,10)<=3); 
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
    seq_score_sig_pas(w)=big_table(row_w,15);
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
    seq_score_sig_mec(w)=big_table(row_w,15);    
end

fraction_total=[ws_fraction_mec,ws_fraction_PaS,ws_fraction_V1]'; %oscillation score
sig_seq_score=[seq_score_sig_mec,seq_score_sig_pas,seq_score_sig_v1]';


%% Figure sequence score VS oscillation score for MEC sessions - Extended data Fig 10q

figure
scatter(ws_fraction_mec,ws_probability_mec,50,'ko','filled');
alpha 0.5
axis([-0.05 1.05 0.2 0.7])
xticks([0 1])
ylabel('Sequence score');
xlabel('Oscillation score');
set(gca,'fontsize',16)
yticks([0.1 0.4 0.7])

%% Figures of histograms total

% Histogram of fractions in adults for all brain areas
edges=[0-0.025:0.05:1+0.025];
figure
h=histcounts(fraction_total,edges);
bar(edges(2:end)-((edges(2)-edges(1))/2),h,1,'FaceColor','k');
ylabel('Counts');
xlabel('Oscillation score')
box off
axis([-0.1 1.1 0 50])
y1=get(gca,'ylim');
hold on
% plot([8/11 8/11],y1,':','linewidth',3,'Color',[139 197 63]./255);
plot([8/11 8/11],y1,':','linewidth',3,'Color','b');
l.linestyle=':';
xticks([0 0.2 0.4 0.6 0.8 1])
set(gca,'fontsize',18,'YColor','k','XColor','k') 


% Histogram of oscillation score for all brain areas - Extended data Fig. 12f
figure
edges=[0-0.025:0.05:1+0.025];
h_mec=histcounts(ws_fraction_mec,edges);
bar(edges(2:end)-((edges(2)-edges(1))/2),h_mec,1,'FaceColor',[34 170 226]/255);
hold on
h_pas=histcounts(ws_fraction_PaS,edges);
bar(edges(2:end)-((edges(2)-edges(1))/2),h_pas,1,'FaceColor',[252 176 60]/255);
h_v1=histcounts(ws_fraction_V1,edges);
bar(edges(2:end)-((edges(2)-edges(1))/2),h_v1,1,'FaceColor',[140 198 54]/255);
alpha 0.5
ylabel('Counts');
xlabel('Oscillation score')
box off
axis([-0.1 1.1 0 25])
y1=get(gca,'ylim');
hold on
% plot([8/11 8/11],y1,':','linewidth',3,'Color',[139 197 63]./255);
plot([8/11 8/11],y1,':','linewidth',3,'Color','b');
l.linestyle=':';
xticks([0 0.2 0.4 0.6 0.8 1])
set(gca,'fontsize',18,'YColor','k','XColor','k') 
legend('MEC','PaS','V1');

% Histogram of oscillation score - Extended data Fig 5d
figure
edges=[0-0.025:0.05:1+0.025];
h=histcounts(ws_fraction_mec,edges);
bar(edges(2:end)-((edges(2)-edges(1))/2),h,1,'FaceColor','k');
% bar(-0.1:0.05:1.1,h,1,'FaceColor','k')
ylabel('Counts');
xlabel('Oscillation score')
box off
axis([-0.1 1.1 0 10])
y1=get(gca,'ylim');
hold on
% plot([8/11 8/11],y1,':','linewidth',3,'Color',[139 197 63]./255);
plot([8/11 8/11],y1,':','linewidth',3,'Color','b');
l.linestyle=':';
xticks([0 0.2 0.4 0.6 0.8 1])
set(gca,'fontsize',18,'YColor','k','XColor','k') 


%% Percentage of sessions with significant sequence score in MEC - with and without oscillations
%Extended data Fig. 10r

mec_waves=find(ws_binary_mec==1);
mec_no_waves=find(ws_binary_mec==0);

seq_score_sig_mec_osc=seq_score_sig_mec(mec_waves);
seq_score_sig_mec_no_osc=seq_score_sig_mec(mec_no_waves);

figure
bar([100*sum(seq_score_sig_mec_osc)/length(seq_score_sig_mec_osc) , 100*sum(seq_score_sig_mec_no_osc)/length(seq_score_sig_mec_no_osc)],'facecolor','k');
set(gca,'fontsize',16,'ycolor','k','xcolor','k');
xticklabels({'MEC/Osc','MEC/No Osc'});
ylabel({'Sessions with significant';'Sequence Score (%)'})
box off


%% Bar figures for oscillatory and non-oscillatory sessions (Figure 6g)

%PaS
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
PaS_sessions_waves=length(find(big_table(PaS_sessions,6)>0));
PaS_sessions_nowaves=length(find(big_table(PaS_sessions,6)==0));
% figure
% bar([1,2],[PaS_sessions_waves,PaS_sessions_nowaves]);
% xticklabels({'Wave','No wave'})
% ylabel('Number of sessions')
% set(gca,'fontsize',20)
% box off

%MEC
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
mec_sessions_waves=length(find(big_table(mec_sessions,6)>0));
mec_sessions_nowaves=length(find(big_table(mec_sessions,6)==0));
% figure
% bar([1,2],[mec_sessions_waves,mec_sessions_nowaves]);
% xticklabels({'Wave','No wave'})
% ylabel('Counts')
% set(gca,'fontsize',16)
% box off

%VIS
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
v1_sessions=sessions;
v1_sessions_waves=length(find(big_table(v1_sessions,6)>0));
v1_sessions_nowaves=length(find(big_table(v1_sessions,6)==0));
% figure
% bar([1,2],[v1_sessions_waves,v1_sessions_nowaves]);
% xticklabels({'Wave','No wave'})
% ylabel('Counts')
% set(gca,'fontsize',16)
% box off


%Full histogram - Figure 6g
figure
bar([1,2],[PaS_sessions_waves,v1_sessions_waves,mec_sessions_waves;PaS_sessions_nowaves,v1_sessions_nowaves,mec_sessions_nowaves]);
xticklabels({'Oscillation','No oscillation'})
ylabel('Number of sessions')
set(gca,'fontsize',20)
box off
legend({'PaS' ; 'V1' ; 'MEC'})
legend boxoff

