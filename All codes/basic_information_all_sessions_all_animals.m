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
                    big_table(count,6)=WS_stat.WS(day,s); %WS - Osc
                    big_table(count,11)=ws_ent.WS_stat.wave_score_ent(day,s); %WS - Ent
                    big_table(count,12)=size(spk.spikes_d_s,1); %N
                    big_table(count,13)=WS_stat.fraction(day,s); %fraction - Osc
                    big_table(count,14)=ws_prob.WS_stat.seq_score_prob(day,s); %WS - Prob


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


%% Figure sequence score VS oscillation score MEC

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


% % Histogram of fractions in adults for all brain areas - overlap figure
figure
edges=[0-0.025:0.05:1+0.025];
h_v1=histcounts(ws_fraction_V1,edges);
bar(edges(2:end)-((edges(2)-edges(1))/2),h_mec,1,'FaceColor',[34 170 226]/255);
hold on
bar(edges(2:end)-((edges(2)-edges(1))/2),h_pas,1,'FaceColor',[252 176 60]/255);
bar(edges(2:end)-((edges(2)-edges(1))/2),h_v1,1,'FaceColor',[140 198 54]/255);
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


% Histogram of fractions in adults for MEC
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


% Histogram of fractions in adults for MEC - overlap figure
figure
edges=[0-0.025:0.05:1+0.025];
h_v1=histcounts(ws_fraction_V1,edges);
bar(edges(2:end)-((edges(2)-edges(1))/2),h_mec,1,'FaceColor',[34 170 226]/255);
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


% Histogram of fractions in adults for PaS - overlap figure
figure
edges=[0-0.025:0.05:1+0.025];
h_v1=histcounts(ws_fraction_V1,edges);
bar(edges(2:end)-((edges(2)-edges(1))/2),h_pas,1,'FaceColor',[252 176 60]/255);
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


% Histogram of fractions in adults for V1 - overlap figure
figure
edges=[0-0.025:0.05:1+0.025];
h_v1=histcounts(ws_fraction_V1,edges);
bar(edges(2:end)-((edges(2)-edges(1))/2),h_v1,1,'FaceColor',[140 198 54]/255);
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




% Histogram of wave score in adults
figure
h=hist(wave_score_total,[-0.1:0.05:1.1]);
bar(-0.1:0.05:1.1,h,1,'FaceColor',[112,128,144]./255)
ylabel('Counts');
xlabel('Binary score')
box off
y1=get(gca,'ylim');
hold on
% plot([8/11 8/11],y1,'k--','linewidth',3);
% l.linestyle='--';
set(gca,'fontsize',18) 
axis([-0.1 1.1 0 80])

% Histogram of waveness in adults
figure
h=hist(waveness_prob_total,0:0.1:3);
bar(0:0.1:3,h,1,'FaceColor',[112,128,144]./255)
ylabel('Counts');
xlabel('Waveness')
box off
y1=get(gca,'ylim');
hold on
% plot([threshold_kl threshold_kl],y1,'k--','linewidth',3);
% l.linestyle='--';
set(gca,'fontsize',18) 
axis([0 3 0 40])



%% Figures of histograms per brain area

%Histogram of fraction MEC
figure
axis([-0.1 1.1 0 15])
h=histcounts(ws_fraction_mec,[0:0.1:1.1]);
bar(-0:0.1:1,h,1,'FaceColor',[112,128,144]./255)
hold on
xlabel('Wave score');
ylabel('Counts (MEC sessions)');
y1=get(gca,'ylim');
hold on
plot([8/11 8/11],[0 15],'k--','linewidth',3);
l.linestyle='--';
set(gca,'fontsize',16)
axis([-0.1 1.1 0 15])
box off


% figure
% scatter(ws_fraction_mec,ws_probability_mec,50,'ko','filled');
% alpha 0.5
% axis([-0.5 1.5 0.1 0.7])
% xticks([0 1])
% xticklabels({'No waves','Waves'});
% ylabel('Waveness');
% set(gca,'fontsize',16)
% yticks([0.1 0.4 0.7])


[p,h,stat]=ranksum(waveness_prob_total(),wave_score_total);

% ------------------------------------------------------------------------

% figure
% histogram(ws_fraction_mec,[-0.1:0.1:1.1])
% hold on 
% histogram(ws_fraction_V1,[-0.1:0.1:1.1])
% histogram(ws_fraction_PaS,[-0.1:0.1:1.1])
% xlabel('Wave score');
% ylabel('Counts');
% legend({'MEC','V1','PaS'});
% plot([8/11 8/11],y1,'k--','linewidth',3);
% l.linestyle='--';
% box off
% set(gca,'fontsize',16)

% 
% figure
% histogram(ws_binary_mec,[-0.1:0.05:1.1])
% hold on 
% histogram(ws_binary_V1,[-0.1:0.05:1.1])
% histogram(ws_binary_PaS,[-0.1:0.05:1.1])
% xlabel('Wave score');
% ylabel('Counts');
% legend({'MEC','V1','PaS'});
% 
% figure
% histogram(ws_entropy_mec,[0:0.1:3]);
% hold on 
% histogram(ws_entropy_V1,[0:0.1:3]);
% histogram(ws_entropy_PaS,[0:0.1:3]);
% xlabel('Waveness');
% ylabel('Counts');
% legend({'MEC','V1','PaS'});

%% Percentage of sessions with significant sequence score in MEC - with and without oscillations

mean_waveness_mec=mean(ws_probability_mec);


mec_waves=find(ws_binary_mec==1);
mec_no_waves=find(ws_binary_mec==0);

%% Comparison sequence score across cortical areas

mean_waveness_mec=mean(ws_probability_mec);
mean_waveness_v1=mean(ws_probability_V1);
mean_waveness_pas=mean(ws_probability_PaS);
sem_waveness_mec=std(ws_probability_mec)/sqrt(length(ws_probability_mec));
sem_waveness_v1=std(ws_probability_V1)/sqrt(length(ws_probability_V1));
sem_waveness_pas=std(ws_probability_PaS)/sqrt(length(ws_probability_PaS));

mec_waves=find(ws_binary_mec==1);
mec_no_waves=find(ws_binary_mec==0);

mean_waveness_mec_waves=mean(ws_probability_mec(mec_waves));
mean_waveness_mec_nowaves=mean(ws_probability_mec(mec_no_waves));
sem_waveness_mec_waves=std(ws_probability_mec(mec_waves))/sqrt(length(mec_waves));
sem_waveness_mec_nowaves=std(ws_probability_mec(mec_no_waves))/sqrt(length(mec_no_waves));



figure
bar([1,2,3],[mean_waveness_pas,mean_waveness_v1,mean_waveness_mec]);
hold on 
errorbar([1,2,3],[mean_waveness_pas,mean_waveness_v1,mean_waveness_mec],[sem_waveness_pas,sem_waveness_v1,sem_waveness_mec],'.k');
xticklabels({'PaS','V1','MEC'});
ylabel('Waveness');
set(gca,'fontsize',18);
box off

figure
bar([1,2],[mean_waveness_mec_nowaves,mean_waveness_mec_waves]);
hold on 
errorbar([1,2],[mean_waveness_mec_nowaves,mean_waveness_mec_waves],[sem_waveness_mec_nowaves,sem_waveness_mec_waves],'.k');
xticklabels({'No waves','Waves'});
ylabel('Waveness MEC (bits)');
set(gca,'fontsize',18);
box off

figure
bar([1,2,3,4],[mean_waveness_pas,mean_waveness_v1,mean_waveness_mec_nowaves,mean_waveness_mec_waves]);
hold on 
errorbar([1,2,3,4],[mean_waveness_pas,mean_waveness_v1,mean_waveness_mec_nowaves,mean_waveness_mec_waves]...
    ,[sem_waveness_pas,sem_waveness_v1,sem_waveness_mec_nowaves,sem_waveness_mec_waves],'.k');
xticklabels({'PaS','V1','MEC NW','MEC W'});
ylabel('Waveness');
set(gca,'fontsize',18);
box off

[p,h]=ranksum(ws_probability_PaS,ws_probability_mec(mec_no_waves))
[p,h]=ranksum(ws_probability_PaS,ws_probability_V1)

%% Comparison waveness-entropy across cortical areas

mean_waveness_mec=mean(ws_entropy_mec);
mean_waveness_v1=mean(ws_entropy_V1);
mean_waveness_pas=mean(ws_entropy_PaS);
sem_waveness_mec=std(ws_entropy_mec)/sqrt(length(ws_entropy_mec));
sem_waveness_v1=std(ws_entropy_V1)/sqrt(length(ws_entropy_V1));
sem_waveness_pas=std(ws_entropy_PaS)/sqrt(length(ws_entropy_PaS));

mec_waves=find(ws_binary_mec==1);
mec_no_waves=find(ws_binary_mec==0);

mean_waveness_mec_waves=mean(ws_entropy_mec(mec_waves));
mean_waveness_mec_nowaves=mean(ws_entropy_mec(mec_no_waves));
sem_waveness_mec_waves=std(ws_entropy_mec(mec_waves))/sqrt(length(mec_waves));
sem_waveness_mec_nowaves=std(ws_entropy_mec(mec_no_waves))/sqrt(length(mec_no_waves));



figure
bar([1,2,3],[mean_waveness_pas,mean_waveness_v1,mean_waveness_mec]);
hold on 
errorbar([1,2,3],[mean_waveness_pas,mean_waveness_v1,mean_waveness_mec],[sem_waveness_pas,sem_waveness_v1,sem_waveness_mec],'.k');
xticklabels({'PaS','V1','MEC'});
ylabel('Waveness');
set(gca,'fontsize',18);
box off

figure
bar([1,2],[mean_waveness_mec_nowaves,mean_waveness_mec_waves]);
hold on 
errorbar([1,2],[mean_waveness_mec_nowaves,mean_waveness_mec_waves],[sem_waveness_mec_nowaves,sem_waveness_mec_waves],'.k');
xticklabels({'No waves','Waves'});
ylabel('Waveness MEC (bits)');
set(gca,'fontsize',18);
box off

figure
bar([1,2,3,4],[mean_waveness_pas,mean_waveness_v1,mean_waveness_mec_nowaves,mean_waveness_mec_waves]);
hold on 
errorbar([1,2,3,4],[mean_waveness_pas,mean_waveness_v1,mean_waveness_mec_nowaves,mean_waveness_mec_waves]...
    ,[sem_waveness_pas,sem_waveness_v1,sem_waveness_mec_nowaves,sem_waveness_mec_waves],'.k');
xticklabels({'PaS','V1','MEC NW','MEC W'});
ylabel('Waveness');
set(gca,'fontsize',18);
box off

[p,h]=ranksum(ws_entropy_PaS,ws_entropy_mec(mec_no_waves))
[p,h]=ranksum(ws_entropy_PaS,ws_entropy_V1)

%% Bar figures for waves no waves

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


figure
bar([1,2],[PaS_sessions_waves,PaS_sessions_nowaves]);
xticklabels({'Wave','No wave'})
ylabel('Number of sessions')
title('Lateral');
set(gca,'fontsize',20)
box off

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

figure
bar([1,2],[mec_sessions_waves,mec_sessions_nowaves]);
xticklabels({'Wave','No wave'})
ylabel('Counts')
title('Medial')
set(gca,'fontsize',16)
box off

%V1

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


figure
bar([1,2],[v1_sessions_waves,v1_sessions_nowaves]);
xticklabels({'Wave','No wave'})
ylabel('Counts')
title('V1')
set(gca,'fontsize',16)
box off


%Full histogram
figure
bar([1,2],[PaS_sessions_waves,v1_sessions_waves,mec_sessions_waves;PaS_sessions_nowaves,v1_sessions_nowaves,mec_sessions_nowaves]);
xticklabels({'Wave','No wave'})
ylabel('Number of sessions')
set(gca,'fontsize',20)
box off
legend({'PaS' ; 'V1' ; 'MEC'})
legend boxoff
% 
% %% Figures for all sessions
% 
% % Histogram of wave scores in adults
% aux2=find(big_table(:,3)>0);
% table_ad=big_table(aux2,:);
% figure
% h=hist(table_ad(:,6),0:0.05:1);
% bar(0:0.05:1,h,'FaceColor',[112,128,144]./255)
% ylabel('Counts');
% xlabel('Wave score')
% box off
% axis([-0.1 1.1 0 60])
% y1=get(gca,'ylim');
% hold on
% plot([threshold_kl threshold_kl],y1,'k--','linewidth',3);
% l.linestyle='--';
% set(gca,'fontsize',18) 
% 
% 
% %Wave sesions
% clus=10;
% waves_ses=find(big_table(:,6)==1); %Sessions inthe table with waves
% adults=find(big_table(:,3)>0); %Sessions inthe table with waves
% waves=intersect(waves_ses,adults);
% count=0;
% 
% % Figures of waves score and ML axis
% % table_aux=big_table;
% % table_aux(find(big_table(:,7)==-10),:)=[];
% aux2=find(big_table(:,3)>15);
% table_ad=big_table(aux2,:);
% figure
% scatter(table_ad(:,8),table_ad(:,6),80,'ko','filled')
% alpha 0.5
% ylabel('Wave score');
% xlabel('ML axis (mm)');
% box off
% set(gca,'fontsize',18)
% axis([-11 5 -0.1 1.1])



%% Figures for adult sessions, without distinguishing between brain areas and FoVs

% Histogram of wave scores in adults for all brain areas
aux2=find(big_table(:,3)>15);
table_ad=big_table(aux2,:);
figure
h=hist(table_ad(:,6),0:0.05:1.4);
bar(0:0.05:1.4,h,5,'FaceColor',[112,128,144]./255)
ylabel('Counts');
% xlabel('Wave score')
xlabel('');
box off
y1=get(gca,'ylim');
hold on
% plot([threshold_kl threshold_kl],y1,'k--','linewidth',3);
% l.linestyle='--';
set(gca,'fontsize',18) 
xticks([0 1]);
xticklabels({'No waves','Waves'})
axis([-0.5 1.5 0 80])

% Histogram of waveness in adults
aux2=find(big_table(:,3)>15);
table_ad=big_table(aux2,:);
figure
h=hist(table_ad(:,14),0:0.1:3);
bar(0:0.1:3,h,1,'FaceColor',[112,128,144]./255)
ylabel('Counts');
xlabel('Waveness')
box off
y1=get(gca,'ylim');
hold on
% plot([threshold_kl threshold_kl],y1,'k--','linewidth',3);
% l.linestyle='--';
set(gca,'fontsize',18) 
axis([0 3 0 40])

% Histogram of fractions in adults
aux2=find(big_table(:,3)>15);
table_ad=big_table(aux2,:);
figure
h=hist(table_ad(:,13),0:0.1:3);
bar(0:0.1:3,h,1,'FaceColor',[112,128,144]./255)
ylabel('Counts');
xlabel('Fraction')
box off
y1=get(gca,'ylim');
hold on
plot([8/11 8/11],y1,'k--','linewidth',3);
l.linestyle='--';
set(gca,'fontsize',18) 
axis([-0.1 2 0 80])

% Histogram of wave score in adults
aux2=find(big_table(:,3)>15);
table_ad=big_table(aux2,:);
figure
h=hist(table_ad(:,6),0:0.1:3);
bar(0:0.1:3,h,1,'FaceColor',[112,128,144]./255)
ylabel('Counts');
xlabel('Wave score')
box off
y1=get(gca,'ylim');
hold on
% plot([8/11 8/11],y1,'k--','linewidth',3);
% l.linestyle='--';
set(gca,'fontsize',18) 
axis([-0.1 2 0 80])



% 
% %Wave sesions
% clus=10;
% waves_ses=find(big_table(:,7)==1); %Sessions inthe table with waves
% adults=find(big_table(:,3)>15); %Sessions inthe table with waves
% waves=intersect(waves_ses,adults);
% count=0;
% 
% % Figures of waves score and ML axis
% % table_aux=big_table;
% % table_aux(find(big_table(:,7)==-10),:)=[];
% aux2=find(big_table(:,3)>15);
% table_ad=big_table(aux2,:);
% figure
% scatter(table_ad(:,10),table_ad(:,6),80,'ko','filled')
% alpha 0.5
% ylabel('Wave score');
% xlabel('ML axis (mm)');
% box off
% set(gca,'fontsize',18)
% axis([-11 5 -0.1 1.14])
% 
% 
% % Figures of waves score and ML axis
% % table_aux=big_table;
% % table_aux(find(big_table(:,7)==-10),:)=[];
% aux2=find(big_table(:,3)>15);
% table_ad=big_table(aux2,:);
% figure
% scatter(table_ad(:,10),table_ad(:,6),80,'ko','filled')
% alpha 0.5
% ylabel('Wave score');
% xlabel('ML axis (mm)');
% box off
% set(gca,'fontsize',18)
% axis([2 4 -0.1 1.4])