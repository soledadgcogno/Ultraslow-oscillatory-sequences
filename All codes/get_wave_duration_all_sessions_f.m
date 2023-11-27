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
    ws_ent=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_',mice(m,:),'.mat']]);
 
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
                    big_table(count,11)=ws_ent.WS_stat.wave_score_prob(day,s); %WS - Ent
                    big_table(count,12)=size(spk.spikes_d_s,1);
                    
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
    clear WS_stat 
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


%% Wave identification

%Identifying waves
sf=7.73; 
duration_waves=nan(100,length(waves));
for w=1:length(waves)
    disp(w);   
    row_w=waves(w);
    count=count+1;
    mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
    day=big_table(waves(w),3);
    s=big_table(waves(w),4);
    munit=big_table(waves(w),5);   
    dt=floor(big_table(waves(w),8));
        
    num_clus_discr=10;
    downsample_factor=1;
    make_fig=0;
    
    file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    load(file_name_spk,'-mat');
    spikes=full(spikes_d_s);

    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);

    spikes_w=[];
    for i=1:size(table_u,1)
        spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
    end
    
    number_waves(w)=size(table_u,1);
    duration_waves(1:number_waves(w),w)=(table_u(:,2)-table_u(:,1))*downsample_factor/sf;
    session_duration(w)=T/sf;
    fraction_time_waves(w)=nansum(duration_waves(:,w))./(T/sf);
    number_cells(w)=size(spikes,1);
    duration_waves_ses{w}=(table_u(:,2)-table_u(:,1))*downsample_factor/sf;
    mean_duration(w)=mean((table_u(:,2)-table_u(:,1))*downsample_factor/sf);
    std_duration(w)=std((table_u(:,2)-table_u(:,1))*downsample_factor/sf);
    sem_duration(w)=std((table_u(:,2)-table_u(:,1))*downsample_factor/sf)/sqrt(number_waves(w));
    mouse_id(w,:)=mouse;
    relevant_time_scale(w)=dt;
    wave_epochs_duration(w)=size(spikes_w,2)/sf;
    fraction_with_waves(w)=size(spikes_w,2)/T;


    clear table_u spikes dt spikes_d_s row_w spikes_w
 
end

% wave_rate_wave_epoch=number_waves./wave_epochs_duration; This can also be
% found in sine fitting
mouse_idx=[1,1,1,2,2,2,3,3,3,3,4,4,4,4,4];

for i=1:15
min_length_s(i)=min(duration_waves_ses{i});
max_length_s(i)=max(duration_waves_ses{i});
end

mean_osc_bin_size=mean(relevant_time_scale);
%% Figures

min_duration=min(min(duration_waves));
max_duration=max(max(duration_waves));

% %Session L9M4 Day 17
% 
% figure
% plot(duration_waves(:,7),'k-*','linewidth',1.5);
% ylabel('Duration (s)');
% xlabel('Wave #');
% yticks([60,120,180]);
% box off
% set(gca,'fontsize',16)
% axis([1 25 70 210])
% 
% %Session L8M2 Day 19
% 
% figure
% plot(duration_waves(:,2),'k-*','linewidth',1.5);
% ylabel('Duration (s)');
% xlabel('Wave #');
% yticks([5,25,50]);
% box off
% set(gca,'fontsize',16)
% axis([0 47 10 60])



%Histogram all sessions
figure
y=histogram(duration_waves(:),15);
bar(y.BinEdges(2:end),y.Values);
ylabel('Counts');
xlabel('Duration (s)');
box off
set(gca,'fontsize',16)
yticks([0 60 120])
xticks([0 50 100 150 200])

% Group by animal
animal1=[duration_waves(:,1);duration_waves(:,2);duration_waves(:,3)];
animal2=[duration_waves(:,4);duration_waves(:,5);duration_waves(:,6)];
animal3=[duration_waves(:,7);duration_waves(:,8);duration_waves(:,9);duration_waves(:,10)];
animal4=[duration_waves(:,11);duration_waves(:,12);duration_waves(:,13);duration_waves(:,14);duration_waves(:,15)];

figure
hold on
scatter(1*ones(1,length(animal1)),animal1,50,'filled','k');
scatter(2*ones(1,length(animal2)),animal2,50,'filled','k');
scatter(3*ones(1,length(animal3)),animal3,50,'filled','k');
scatter(4*ones(1,length(animal4)),animal4,50,'filled','k');
alpha 0.3
axis([-0 5 0 220]);
xticks([1 2 3 4]);
yticks([0 100 200]);
ylabel('Duration (s)');
xlabel('Animal #');
set(gca,'fontsize',16);
xticklabels({'L8m2','L9m1','L9m4','L5m5'})

% Variability of duration
figure
hold on
for i=1:length(waves)
errorbar(1:length(waves),mean_duration,std_duration,'k','linewidth',2);
end
ylabel('Wave duration');
xlabel('Session #');
set(gca,'fontsize',16);
axis([0 16 20 160])

%colormap
cc=colorcube(30);
 for i=1:length(waves)
     if mouse_idx(i)==1
         color_mouse(i,:)=cc(9,:);
     elseif mouse_idx(i)==2
         color_mouse(i,:)=cc(12,:);
     elseif mouse_idx(i)==3
         color_mouse(i,:)=cc(14,:);
     elseif mouse_idx(i)==4
           color_mouse(i,:)=cc(21,:);
     end
 end

% Number of cells and duration
a=[];
b=[];
figure
hold on
for i=1:length(waves)
    h=scatter(number_cells(i)*ones(1,length(duration_waves_ses{i})),duration_waves_ses{i},[],color_mouse(i,:),'filled');
    alpha 0.8
    a=[a,number_cells(i)*ones(1,length(duration_waves_ses{i}))];
    b=[b;duration_waves_ses{i}];
    
    if (i==1 || i==4 || i==7 || i==11);
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
ylabel('Wave duration (s)');
xlabel('Number of imaged cells');
set(gca,'fontsize',16);
legend('Animal 1','Animal 2','Animal 3','Animal 4');
legend boxoff

[rho,pval]=corr(a',b,'type','spearman');

% Wave rate
figure
hold on
h=plot(number_waves./(session_duration./60),'-k','markersize',1,'linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
for w=1:15
    h=plot(w,number_waves(w)./(session_duration(w)./60),'o','markersize',10,'color',color_mouse(w,:),'markerfacecolor',color_mouse(w,:));
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Wave rate (waves/minute)');
xlabel('Session #')
axis([0 16 0 1.8])
xticks([1 7 15])
box off
legend('Animal 1','Animal 2','Animal 3','Animal 4');
legend boxoff


% Number of waves VS session duration

figure
hold on
for w=1:15
h=scatter(session_duration(w)/60,number_waves(w),60,color_mouse(w,:),'filled');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Number of waves per session');
xlabel('Session duration (minutes)')
axis([25 65 0 80])
xticks([30 45 60])
box off
legend('Animal 1','Animal 2','Animal 3','Animal 4');
legend boxoff
[rho,pval]=corr(session_duration',number_waves','type','spearman');

short=find(session_duration/60<40);
long=find(session_duration/60>40);

[p,h,stst]=ranksum(number_waves(short)',number_waves(long)');


mat=nan(8,2);
mat(:,1)=number_waves(short)';
mat(1:7,2)=number_waves(long)';

figure
boxplot(mat)
ylabel('Number of cycles per session');
xticks([1 2])
xticklabels({'30','60'});
xlabel('Session length (minutes)');
set(gca,'fontsize',16,'YColor','k','XColor','k');
box off

% Session length and mean cycle wave duration
figure
hold on
for w=1:15
h=errorbar(session_duration(w)/60,mean_duration(w),std_duration(w),'o','markersize',6,'linewidth',1.5,'MarkerEdgeColor',color_mouse(w,:),'MarkerFaceColor',color_mouse(w,:));
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    h.Color=color_mouse(w,:);
end
ylabel('Cycle length (s)');
xlabel('Session length (minutes)');
set(gca,'fontsize',16);
axis([25 65 0 200])
xticks([30 45 60])
box off
legend('Animal 1','Animal 2','Animal 3','Animal 4');
legend boxoff

[rho,pval]=corr((session_duration/60)',(mean_duration)','type','spearman');

short=find(session_duration/60<40);
long=find(session_duration/60>40);

mat=nan(8,2);
mat(:,1)=mean_duration(short)';
mat(1:7,2)=mean_duration(long)';

figure
boxplot(mat)
ylabel('Mean vycle length (s)');
xticks([1 2])
xticklabels({'30','60'});
xlabel('Session length (minutes)');
set(gca,'fontsize',16,'YColor','k','XColor','k');
box off


[p,h,stst]=ranksum(mean_duration(short)',mean_duration(long)','tail','left');

% Wave rate vs 1/duration
figure
set(gca,'fontsize',16);
hold on
for w=1:15
h=scatter(1./(mean_duration(w)/60),number_waves(w)./(session_duration(w)./60),60,color_mouse(w,:),'filled');

    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
l=refline(1,0);
l.LineStyle='--';
l.Color='k';
l.LineWidth=1.5;
ylabel('Wave rate (waves/minute)');
xlabel('1/(mean wave duration) (1/minute)');
box off
legend('Animal 1','Animal 2','Animal 3','Animal 4');
legend boxoff

% Relevant time scale
bin_size=1/sf;
figure
hold on
h=plot(relevant_time_scale/sf,'-k','markersize',1,'linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
for w=1:15
    h=plot(w,relevant_time_scale(w)/sf,'o','markersize',10,'color',color_mouse(w,:),'markerfacecolor',color_mouse(w,:));
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Relevant time scale (s)');
xlabel('Session #')
axis([0 16 2 18])
yticks([2 10 18])
box off
legend('Animal 1','Animal 2','Animal 3','Animal 4');
legend boxoff


% RNumber of sequences

figure
scatter(1:3,[length(duration_waves_ses{1}),length(duration_waves_ses{2}),length(duration_waves_ses{3})],60,color_mouse(1,:),'filled');
hold on
scatter(1:3,[length(duration_waves_ses{4}),length(duration_waves_ses{5}),length(duration_waves_ses{6})],60,color_mouse(4,:),'filled');
scatter(1:4,[length(duration_waves_ses{7}),length(duration_waves_ses{8}),length(duration_waves_ses{9}),length(duration_waves_ses{10})],60,color_mouse(7,:),'filled');
scatter(1:5,[length(duration_waves_ses{11}),length(duration_waves_ses{12}),length(duration_waves_ses{13}),length(duration_waves_ses{14}),length(duration_waves_ses{15})],60,color_mouse(11,:),'filled');
xlim([0 6]);
ylim([0 80]);
ylabel('Total number of cycles per session');
xlabel('Session #');
legend()
legend('60355','60584','60585','59914');
legend boxoff
set(gca,'fontsize',16,'YColor','k','XColor','k');


%% Figure of wave progression for one example session

w=8;

dur_waves=duration_waves_ses{w};
mean_dur=mean(dur_waves);
std_dur=std(dur_waves);

% dur_waves_113=(table_113(:,2)-table_113(:,1))*downsample_factor/sf;
% dur_waves_117=(table_117(:,2)-table_117(:,1))*downsample_factor/sf;

figure
scatter(1:length(dur_waves),dur_waves,60,'black','filled');
axis([0 length(dur_waves)+1 0 250 ]);
ylabel('Duration of individual waves (s)');
xlabel('Wave index #');
set(gca,'fontsize',16);

figure
boxplot(dur_waves);
ylabel('Cycle length (s)');
% xlabel('Ensemble #')
set(gca,'fontsize',16,'YColor','k','XColor','k');
xticks([1])
xticklabels({'Cycles in example session'});
ylim([0 250])
% yticks([])
box off

dur_waves_mediansubstracted=dur_waves-median(dur_waves);
[p,h,stats] = signrank(dur_waves_mediansubstracted);


pct_change=(max(dur_waves)-min(dur_waves))/min(dur_waves);

% mean(dur_waves);
% std(dur_waves)/sqrt(length(dur_waves));
% 
% min_duration_1session=min(dur_waves);
% max_duration_1session=max(dur_waves);


%% Figures that I separate per animal

% Wave rate

figure
hold on
h=plot(1:3,number_waves(1:3)./(session_duration(1:3)./60),'k','linewidth',2);
for w=1:3
    hold on
   	plot(w,number_waves(w)./(session_duration(w)./60),'-o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Wave rate (waves/minute)');
xlabel('Session #')
axis([0.5 3.5 0 1.8])
xticks([1 2 3])
box off

figure
hold on
h=plot(1:3,number_waves(4:6)./(session_duration(4:6)./60),'k','linewidth',2);
for w=4:6
    hold on
    plot(w-3,number_waves(w)./(session_duration(w)./60),'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Wave rate (waves/minute)');
xlabel('Session #')
axis([0.5 3.5 0 1.8])
xticks([1 2 3])
box off

figure
hold on
h=plot(1:4,number_waves(7:10)./(session_duration(7:10)./60),'k','linewidth',2);
for w=7:10
%     h=plot(w-6,number_waves(w)./(session_duration(w)./60),'o','markersize',10,'color',color_mouse(w,:),'markerfacecolor',color_mouse(w,:));
    h=plot(w-6,number_waves(w)./(session_duration(w)./60),'o','markersize',10,'color','k','markerfacecolor','k');

    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Wave rate (waves/minute)');
xlabel('Session #')
axis([0.5 4.5 0 1.8])
xticks([1 2 3 4])
box off

figure
hold on
h=plot(1:5,number_waves(11:15)./(session_duration(11:15)./60),'k','linewidth',2);
for w=11:15
%     h=plot(w-10,number_waves(w)./(session_duration(w)./60),'o','markersize',10,'color',color_mouse(w,:),'markerfacecolor',color_mouse(w,:));
    h=plot(w-10,number_waves(w)./(session_duration(w)./60),'o','markersize',10,'color','k','markerfacecolor','k');

    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Wave rate (waves/minute)');
xlabel('Session #')
axis([0.5 5.5 0 1.8])
xticks([1 2 3 4 5])
box off


% Relevant time scale
bin_size=1/sf;

figure
hold on
h=plot(relevant_time_scale(1:3)/sf,'-k','markersize',1,'linewidth',2);
for w=1:3
    h=plot(1:3,relevant_time_scale(1:3)/sf,'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16,'YColor','k','XColor','k');
ylabel('Oscillation bin size (s)');
xlabel('Session #')
axis([0.5 3.5 0 20])
yticks([0 10 20])
xticks([1 2 3 ])
box off

figure
hold on
h=plot(relevant_time_scale(4:6)/sf,'-k','markersize',1,'linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
for w=4:6
    h=plot(1:3,relevant_time_scale(4:6)/sf,'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16,'YColor','k','XColor','k');
ylabel('Oscillation bin size (s)');
xlabel('Session #')
axis([0.5 3.5 0 20])
yticks([0 10 20])
xticks([1 2 3 ])
box off

figure
hold on
h=plot(relevant_time_scale(7:10)/sf,'-k','markersize',1,'linewidth',2);
for w=7:10
    h=plot(1:4,relevant_time_scale(7:10)/sf,'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16,'YColor','k','XColor','k');
ylabel('Oscillation bin size (s)');
xlabel('Session #')
axis([0.5 4.5 0 20])
yticks([0 10 20])
xticks([1 2 3 4 ])
box off

figure
hold on
h=plot(relevant_time_scale(11:15)/sf,'-k','markersize',1,'linewidth',2);
for w=11:15
    h=plot(1:5,relevant_time_scale(11:15)/sf,'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16,'YColor','k','XColor','k');
ylabel('Oscillation bin size (s)');
xlabel('Session #')
axis([0.5 5.5 0 20])
yticks([0 10 20])
xticks([1 2 3 4 5])
box off

% Fraction of session with waves


figure
hold on
h=plot(fraction_with_waves(1:3),'-k','markersize',1,'linewidth',2);
for w=1:3
    h=plot(1:3,fraction_with_waves(1:3),'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Fraction of session with oscillation');
xlabel('Session #')
axis([0.5 3.5 0 1])
yticks([0 0.5 1])
xticks([1 2 3 ])
box off

figure
hold on
h=plot(fraction_with_waves(4:6),'-k','markersize',1,'linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
for w=4:6
    h=plot(1:3,fraction_with_waves(4:6),'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Fraction of session with oscillation');
xlabel('Session #')
axis([0.5 3.5 0 1])
yticks([0 0.5 1])
xticks([1 2 3 ])
box off

figure
hold on
h=plot(fraction_with_waves(7:10),'-k','markersize',1,'linewidth',2);
for w=7:10
    h=plot(1:4,fraction_with_waves(7:10),'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Fraction of session with oscillation');
xlabel('Session #')
axis([0.5 4.5 0 1])
yticks([0 0.5 1])
xticks([1 2 3 4 ])
box off

figure
hold on
h=plot(fraction_with_waves(11:15),'-k','markersize',1,'linewidth',2);
for w=11:15
    h=plot(1:5,fraction_with_waves(11:15),'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Fraction of session with oscillation');
xlabel('Session #')
axis([0.5 5.5 0 1])
yticks([0 0.5 1])
xticks([1 2 3 4 5])
box off


% Duration of combined wave epochs

figure
hold on
h=plot(wave_epochs_duration(1:3)/60,'-k','markersize',1,'linewidth',2);
for w=1:3
    h=plot(1:3,wave_epochs_duration(1:3)/60,'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Duration of oscillatory epoch (min)');
xlabel('Session #')
axis([0.5 3.5 20 30])
yticks([20 25 30])
xticks([1 2 3 ])
box off

figure
hold on
h=plot(wave_epochs_duration(4:6)/60,'-k','markersize',1,'linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
for w=4:6
    h=plot(1:3,wave_epochs_duration(4:6)/60,'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Duration of oscillatory epoch (min)');
xlabel('Session #')
axis([0.5 3.5 15 60])
yticks([20 40 60])
xticks([1 2 3 ])
box off

figure
hold on
h=plot(wave_epochs_duration(7:10)/60,'-k','markersize',1,'linewidth',2);
for w=7:10
    h=plot(1:4,wave_epochs_duration(7:10)/60,'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Duration of oscillatory epoch (min)');
xlabel('Session #')
axis([0.5 4.5 15 60])
yticks([20 40 60])
xticks([1 2 3 4 ])
box off

figure
hold on
h=plot(wave_epochs_duration(11:15)/60,'-k','markersize',1,'linewidth',2);
for w=11:15
    h=plot(1:5,wave_epochs_duration(11:15)/60,'o','markersize',10,'color','k','markerfacecolor','k');
    if (w==1 || w==4 || w==7 || w==11)
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
set(gca,'fontsize',16);
ylabel('Duration of oscillatory epoch (min)');
xlabel('Session #')
axis([0.5 5.5 5 25])
yticks([5 15 25])
xticks([1 2 3 4 5])
box off

% Cycle length 
%mouse 1

dur_waves_1=duration_waves_ses{1};
dur_waves_2=duration_waves_ses{2};
dur_waves_3=duration_waves_ses{3};

mat_dur=nan(max([length(dur_waves_1),length(dur_waves_2),length(dur_waves_3)]),3);
mat_dur(1:length(dur_waves_1),1)=dur_waves_1';
mat_dur(1:length(dur_waves_2),2)=dur_waves_2';
mat_dur(1:length(dur_waves_3),3)=dur_waves_3';

figure
boxplot(mat_dur);
ylabel('Cycle length (s)');
xlabel('Session #');
set(gca,'fontsize',16,'YColor','k','XColor','k');
xticks([1,2,3])
% xticklabels({'Session #'});
ylim([0 120])
% yticks([])
box off

clear mat_dur dur_waves_1 dur_waves_2 dur_waves_3

%mouse 2
dur_waves_1=duration_waves_ses{4};
dur_waves_2=duration_waves_ses{5};
dur_waves_3=duration_waves_ses{6};

mat_dur=nan(max([length(dur_waves_1),length(dur_waves_2),length(dur_waves_3)]),3);
mat_dur(1:length(dur_waves_1),1)=dur_waves_1';
mat_dur(1:length(dur_waves_2),2)=dur_waves_2';
mat_dur(1:length(dur_waves_3),3)=dur_waves_3';

figure
boxplot(mat_dur);
ylabel('Cycle length (s)');
xlabel('Session #');
set(gca,'fontsize',16,'YColor','k','XColor','k');
xticks([1,2,3])
% xticklabels({'Session #'});
ylim([0 200])
% yticks([])
box off

clear mat_dur dur_waves_1 dur_waves_2 dur_waves_3

%mouse 3
dur_waves_1=duration_waves_ses{7};
dur_waves_2=duration_waves_ses{8};
dur_waves_3=duration_waves_ses{9};
dur_waves_4=duration_waves_ses{10};

mat_dur=nan(max([length(dur_waves_1),length(dur_waves_2),length(dur_waves_3),length(dur_waves_4)]),4);
mat_dur(1:length(dur_waves_1),1)=dur_waves_1';
mat_dur(1:length(dur_waves_2),2)=dur_waves_2';
mat_dur(1:length(dur_waves_3),3)=dur_waves_3';
mat_dur(1:length(dur_waves_4),4)=dur_waves_4';

figure
boxplot(mat_dur);
ylabel('Cycle length (s)');
xlabel('Session #');
set(gca,'fontsize',16,'YColor','k','XColor','k');
xticks([1,2,3,4])
% xticklabels({'Session #'});
ylim([0 250])
% yticks([])
box off

clear mat_dur dur_waves_1 dur_waves_2 dur_waves_3 dur_waves_4

%mouse 4
dur_waves_1=duration_waves_ses{11};
dur_waves_2=duration_waves_ses{12};
dur_waves_3=duration_waves_ses{13};
dur_waves_4=duration_waves_ses{14};
dur_waves_5=duration_waves_ses{15};

mat_dur=nan(max([length(dur_waves_1),length(dur_waves_2),length(dur_waves_3),length(dur_waves_4),length(dur_waves_5)]),5);
mat_dur(1:length(dur_waves_1),1)=dur_waves_1';
mat_dur(1:length(dur_waves_2),2)=dur_waves_2';
mat_dur(1:length(dur_waves_3),3)=dur_waves_3';
mat_dur(1:length(dur_waves_4),4)=dur_waves_4';
mat_dur(1:length(dur_waves_5),5)=dur_waves_5';

figure
boxplot(mat_dur);
ylabel('Cycle length (s)');
xlabel('Session #');
set(gca,'fontsize',16,'YColor','k','XColor','k');
xticks([1,2,3,4,5])
% xticklabels({'Session #'});
ylim([0 250])
% yticks([])
box off

clear mat_dur dur_waves_1 dur_waves_2 dur_waves_3 dur_waves_4 dur_waves_5

%% Quantification of cycle length within VS between
duration_pooled=[];
for w=1:15
duration_pooled=[duration_pooled;duration_waves_ses{w}];
end

for w=1:15
    fraction_s(w)=min(duration_waves_ses{w})/max(duration_waves_ses{w});
    std_s(w)=std(duration_waves_ses{w});
end

for sh=1:500
    duration_pooled_c=duration_pooled;
    for w=1:15
        %disp(length(duration_pooled_c))
        in=randperm(length(duration_pooled_c),number_waves(w));
        cycles_s=duration_pooled_c(in);
        fraction_s_sh(w,sh)=min(cycles_s)/max(cycles_s);
        std_s_sh(w,sh)=std(cycles_s);

        duration_pooled_c(in)=[];

        clear cycles_s in 
    end
end

mat_cyclelength=nan(length(fraction_s_sh(:)),2);
mat_cyclelength(1:15,1)=fraction_s(:);
mat_cyclelength(1:end,2)=fraction_s_sh(:);

% figure
% boxplot(mat_cyclelength)
% xticklabels({'Data','Shuffle'});
% title('fraction ')
% ylim([0 1])
% [p,h,stats]=ranksum(fraction_s(:),fraction_s_sh(:));

clear mat_cyclelength
mat_cyclelength=nan(length(std_s_sh(:)),2);
mat_cyclelength(1:15,1)=std_s(:);
mat_cyclelength(1:end,2)=std_s_sh(:);

figure
boxplot(mat_cyclelength)
xticklabels({'Data','Shuffle'});
ylabel('S.D. of cycle length within a session (s) ');
set(gca,'fontsize',16,'Ycolor','k','Xcolor','k');
ylim([0 90])
box off
[p,h,stats]=ranksum(std_s(:),std_s_sh(:),'tail','left');
[p,h,stats]=ranksum(std_s(:),std_s_sh(:));

%--------------------- Another quantification
count_w=0;
for w=1:15
    counts=0;
    clear delta_t_w_s
    for c=1:number_waves(w)
        for k=c+1:(number_waves(w))
            count_w=count_w+1;
            counts=counts+1;

            delta_t_w(count_w)=min(duration_waves_ses{w}(c),duration_waves_ses{w}(k))...
                /max(duration_waves_ses{w}(c),duration_waves_ses{w}(k));
            delta_t_w_s(counts)=delta_t_w(count_w);
        end
    end
    delta_w_mean(w)=mean(delta_t_w_s);
end

count_b=0;
for w=1:14
    counts=0;
    clear delta_t_b_s
    for c=1:number_waves(w)
        for w2=w+1:15
            if w2~=w
                for k=1:(number_waves(w2))
                    count_b=count_b+1;
                    counts=counts+1;

                    delta_t_b(count_b)=min(duration_waves_ses{w}(c),duration_waves_ses{w2}(k))...
                        /max(duration_waves_ses{w}(c),duration_waves_ses{w2}(k));
                    delta_t_b_s(counts)=delta_t_b(count_b);

                end
            end
        end
    end
    delta_b_mean(w)=mean(delta_t_b_s);
end

mat_deltat=nan(length(delta_t_b),2);
mat_deltat(1:length(delta_t_w),1)=delta_t_w(:);
mat_deltat(1:end,2)=delta_t_b(:);

figure
boxplot(mat_deltat);
xticklabels({'Within session','Between sessions'});
ylabel({'Fraction between';'shortest and longest cycle length'});
ylim([0 1]);
set(gca,'fontsize',16,'Ycolor','k','Xcolor','k');
box off

[p,h,stats]=ranksum(delta_t_w(:),delta_t_b(:),'tail','right');
[p,h,stats]=ranksum(delta_t_w(:),delta_t_b(:));
[h,p,stats]=ttest2(delta_t_w(:),delta_t_b(:));

%--------------------- Here I change the normalization for the second quantification

count_w=0;
for w=1:15
    counts=0;
    clear delta_t_w_s
    for c=1:number_waves(w)
        for k=c+1:(number_waves(w))
            count_w=count_w+1;
            counts=counts+1;

            delta_t_w(count_w)=min(duration_waves_ses{w}(c),duration_waves_ses{w}(k))...
                /max(duration_waves_ses{w}(c),duration_waves_ses{w}(k));
            delta_t_w_s(counts)=delta_t_w(count_w);
        end
    end
    delta_w_mean(w)=mean(delta_t_w_s);
end

count_b=0;
for w=1:15
    counts=0;
    clear delta_t_b_s
    for c=1:number_waves(w)
        for w2=1:15
            if w2~=w
                for k=1:(number_waves(w2))
                    count_b=count_b+1;
                    counts=counts+1;

                    delta_t_b(count_b)=min(duration_waves_ses{w}(c),duration_waves_ses{w2}(k))...
                        /max(duration_waves_ses{w}(c),duration_waves_ses{w2}(k));
                    delta_t_b_s(counts)=delta_t_b(count_b);

                end
            end
        end
    end
    delta_b_mean(w)=mean(delta_t_b_s);

end

mat_deltat=nan(length(delta_w_mean),2);
mat_deltat(1:length(delta_w_mean),1)=delta_w_mean(:);
mat_deltat(1:end,2)=delta_b_mean(:);

figure
boxplot(mat_deltat);
xticklabels({'Within session','Between sessions'});
ylabel({'Fraction between';'shortest and longest cycle length'});
ylim([0 1]);
set(gca,'fontsize',16,'Ycolor','k','Xcolor','k');
box off

[p,h,stats]=ranksum(delta_w_mean(:),delta_b_mean(:),'tail','right');
[p,h,stats]=ranksum(delta_w_mean(:),delta_b_mean(:));
[h,p,stats]=ttest2(delta_w_mean(:),delta_b_mean(:));

%% Firing rate per epoch of the wave - Normalization 1 - Equal intervals


%Identifying waves
sf=7.73; 
% N_p=3;
% figure
p_1=nan(length(waves),7);
p_2=nan(length(waves),7);
for N_p=10

    sum_spikes=nan(100,N_p,length(waves));
    sum_spikes_n2=nan(100,N_p,length(waves));
    sum_spikes_n3=nan(100,N_p,length(waves));
    sum_spikes_n4=nan(100,N_p,length(waves));
    sum_spikes_n5=nan(100,N_p,length(waves));

    for w=1:length(waves)
        disp(w);
        row_w=waves(w);
        count=count+1;
        mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
        day=big_table(waves(w),3);
        s=big_table(waves(w),4);
        munit=big_table(waves(w),5);
        dt=floor(big_table(waves(w),8));

        num_clus_discr=10;
        downsample_factor=1;
        make_fig=0;

        file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
        load(file_name_spk,'-mat');
        spikes=full(spikes_d_s);
        [~,sorting_w,~]=get_sorting(spikes);

        [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);

        spikes_r=[];
        for wa=1:size(table_u,1)
            spikes_r=[spikes_r,spikes(:,table_u(wa,1):table_u(wa,2))];
            segment=floor((table_u(wa,2)-table_u(wa,1))/N_p);

            for p=1:N_p
                if p == N_p
                    aux=[table_u(wa,1)+(p-1)*segment : table_u(wa,2)];
                else
                    aux=[table_u(wa,1)+(p-1)*segment : table_u(wa,1)+p*segment-1];
                end

%                 sum_spikes(wa,p,w)=sum(sum(spikes(:,aux)));
%                 sum_spikes_n2(wa,p,w)=sum(sum(spikes(:,aux)))./sum(sum(spikes(:,table_u(wa,1):table_u(wa,2))));
%                 sum_spikes_n3(wa,p,w)=sum(sum(spikes(:,aux)))./(length(aux)*(1/sf));
                sum_spikes_n4(wa,p,w)=sum(sum(spikes(:,aux)))./(N*(length(aux))*(1/sf));
%                 sum_spikes_n5(wa,p,w)=sum(sum(spikes(:,aux)))./N;

                clear aux
            end
            clear segment
        end

        sum_spikes_n4_average(w,:)=nanmean(sum_spikes_n4(:,:,w));

        FR(w)=sum(sum(spikes_r))/(size(spikes_r,2)*(1/sf)*N);
        sum_spikes_n6(:,:,w)=sum_spikes_n4(:,:,w)/FR(w);

        sum_spikes_n6_average(w,:)=nanmean(sum_spikes_n6(:,:,w));

        % p_1(w) = anova1(sum_spikes_n6(:,:,w));


        clear table_u spikes dt spikes_d_s row_ws spikes_r 

    end
    clear table_u spikes dt spikes_d_s row_ws spikes_r sum_spikes_n sum_spikes
end


% Simple normalization with multcompare n=15
figure
boxplot(sum_spikes_n4_average)
ylabel({'Mean firing rate (Hz)'});
xlabel('Sequence interval');
set(gca,'fontsize',16);
axis([0.5 10.5 0 0.28])
box off

title('n=15')
[p_value,~,stats] = friedman(sum_spikes_n4_average,1);
% [p_value,~,stats] = kruskalwallis(sum_spikes_n4_average);

[results,means] = multcompare(stats, 'Alpha',0.05,'CType','bonferroni');
pvalues=results(:,6);
alpha_adj=0.05;
nonsig=length(find(pvalues>alpha_adj));
sig=length(find(pvalues<=alpha_adj));


[results,~,stats]=multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');

clear p p_value stats results means pvalues nonsig sig

%percentage change
means=mean(sum_spikes_n4_average);
pc_change=(max(means)-min(means))/min(means);
% Simple normalization with wilcoxon n=15

count=0;
for i=1:10
    for j=i+1:10
        count=count+1;
        [p(count)]=ranksum(sum_spikes_n4_average(:,i),sum_spikes_n4_average(:,j));
    end
end
nonsig=length(find(p>(0.05/45)));

% Simple normalization with multcompare n=421
sum_spikes_n_r=[];
for i=1:15
    sum_spikes_n_r=[sum_spikes_n_r;sum_spikes_n4(:,:,i)];
end

aux=find(isnan(sum_spikes_n_r(:,1))==1);
sum_spikes_n_r(aux,:)=[];

figure
boxplot(sum_spikes_n_r)
ylabel('Normalized activity during wave interval');
xlabel('Interval');

[p_value_10,~,stats_10] = friedman(sum_spikes_n_r,1);

%[p_value_10,~,stats_10] = kruskalwallis(sum_spikes_n_r);
[results,~] = multcompare(stats_10, 'Alpha',0.05,'CType','bonferroni');
pvalues=results(:,6);
alpha_adj=0.05;
nonsig=length(find(pvalues>alpha_adj));
sig=length(find(pvalues<=alpha_adj));

[results]=multcompare(stats_10,'Alpha',0.05,'CType','tukey-kramer')

% Simple normalization with wilcoxon n=421
clear p p_value stats results means pvalues nonsig sig

count=0;
for i=1:10
    for j=i+1:10
        count=count+1;
        [p(count)]=ranksum(sum_spikes_n_r(:,i),sum_spikes_n_r(:,j));
    end
end
nonsig=length(find(p>alpha_adj));


% Normalization by firing rate
% 
% figure
% boxplot(sum_spikes_n6_average)
% ylabel({'Normalized activity';'during wave interval (Hz)'});
% xlabel('Wave interval');
% set(gca,'fontsize',16);
% box off
% 
% [p_value,~,stats] = kruskalwallis(sum_spikes_n6_average);
% [results,means] = multcompare(stats, 'Alpha',0.05,'CType','bonferroni');
% pvalues=results(:,6);
% alpha_adj=0.05/45;
% 
% nonsig=length(find(pvalues>alpha_adj));
% sig=length(find(pvalues<=alpha_adj));


% 
% % 
% sum_spikes_n_r=[];
% for i=1:15
%     sum_spikes_n_r=[sum_spikes_n_r;sum_spikes_n6(:,:,i)];
% end
% 
% aux=find(isnan(sum_spikes_n_r(:,1))==1);
% sum_spikes_n_r(aux,:)=[];
% figure
% boxplot(sum_spikes_n_r)
% ylabel('Normalized activity during wave interval');
% xlabel('Interval');
% 
% 
% [p_value_10,~,stats_10] = anova1(sum_spikes_n_r);
% [results,means] = multcompare(stats_10, 'Alpha',0.05,'CType','bonferroni');
% pvalues=results(:,6);
% nonsig=length(find(pvalues>0.05));
% sig=length(find(pvalues<=0.05));

% %% Firing rate per epoch of the wave - Normalization 2 - Using the phase of the wave
% 
% %Identifying waves
% sf=7.73;
% phase_bins=[-pi:2*pi/10:pi];
% 
% N_p=10;
% sum_spikes=nan(100,N_p,length(waves));
% sum_spikes_n=nan(100,N_p,length(waves));
% sum_spikes_n2=nan(100,N_p,length(waves));
% sum_spikes_n4=nan(100,N_p,length(waves));
% sum_spikes_n5=nan(100,N_p,length(waves));
% sum_spikes_n6=nan(100,N_p,length(waves));
% sum_spikes_n7=nan(100,N_p,length(waves));
% 
% % length_interval=nan(100,N_p,length(waves));
% % sum_spikes_n3=nan(200000,10);
% % table_phase_spikes=[];
% % countp=zeros(1,10);
% 
% for w=1:length(waves)
%     disp(w);
%     row_w=waves(w);
% 
%     count=count+1;
%     mouse=['L',num2str(big_table(waves(w),1)),'M',num2str(big_table(waves(w),2))];
%     day=big_table(waves(w),3);
%     s=big_table(waves(w),4);
%     munit=big_table(waves(w),5);
%     dt=floor(big_table(waves(w),8));
% 
%     num_clus_discr=10;
%     downsample_factor=1;
%     make_fig=0;
% 
%     file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
%     load(file_name_spk,'-mat');
%     spikes=full(spikes_d_s);
%     %     [~,sorting_w,~]=get_sorting(spikes);
%     N=size(spikes,1);
% 
%     clear FRp
%     for i=1:N
%         FRp(i,:)=full(fire_rate(spikes(i,:),1*dt,'g')); %smooth using as kernel the dt chosen for each session
%     end
% 
%     [~,scoret,~] = pca(FRp');
%     phase_f=(atan2(smooth(scoret(:,2),floor(1*dt)),smooth(scoret(:,1),floor(1*dt))));
% 
%     [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);
% 
%     spikes_r=[];
% %     phase_r=[]; %Reduced phase; only contains wave epochs
% 
%     for wa=1:size(table_u,1)
%         spikes_r=[spikes_r,spikes(:,table_u(wa,1):table_u(wa,2))];
% 
%         phase_d=discretize(phase_f(table_u(wa,1):table_u(wa,2)),phase_bins);
%         segment=spikes(:,table_u(wa,1):table_u(wa,2));
% 
% %         spikes_segment=sum(segment)/sum(sum(segment));
% 
% for p=1:N_p
% 
%     aux=find(phase_d==p);
%     if aux>2
%         %             sum_spikes(wa,p,w)=sum(sum(segment(:,aux)));
%         %             sum_spikes_n2(wa,p,w)=sum(sum(segment(:,aux)))/sum(sum(segment));
%         %             sum_spikes_n4(wa,p,w)=sum(sum(segment(:,aux)))/(length(aux)*(1/sf));
%         %             sum_spikes_n5(wa,p,w)=sum(sum(segment(:,aux)))/sum(sum(segment))/(length(aux)*(1/sf));
%         sum_spikes_n6(wa,p,w)=sum(sum(segment(:,aux)))/(N*(length(aux)*(1/sf)));
% 
%         if(sum_spikes_n6(wa,p,w)>2)
%             disp([w,wa]);
%         end
% 
%         length_interval(wa,p,w)=(length(aux)*(1/sf));
%         spikes_interval(wa,p,w)=sum(sum(segment));
%     end
%     clear aux
% end
% 
%         clear segment phase_d
% 
%     end
% 
%     sum_spikes_n6_average(w,:)=nanmean(sum_spikes_n6(:,:,w),1);
% 
%     FR(w)=sum(sum(spikes_r))/(N*size(spikes_r,2)*(1/sf));
%     sum_spikes_n7(:,:,w)=sum_spikes_n6(:,:,w)/FR(w);
% 
%     sum_spikes_n7_average(w,:)=nanmean(sum_spikes_n7(:,:,w));
% 
% 
%     clear table_u spikes dt spikes_d_s row_ws spikes_r sum_spikes_n sum_spikes FRp phase_f sorting_w ...
%         phase_r coefft scoret
% end
% 
% % Simple normalization n=15 - Wilcoxon
% 
% figure
% boxplot(sum_spikes_n6_average)
% ylabel({'Normalized activity';'during wave interval (Hz)'});
% xlabel('Phase interval');
% set(gca,'fontsize',16);
% box off
% 
% p = friedman(sum_spikes_n6_average,1);
% 
% p = kruskalwallis(sum_spikes_n6_average);
% count=0;
% for i=1:10
%     for j=i+1:10
%         count=count+1;
%         [p(count)]=ranksum(sum_spikes_n6_average(:,i),sum_spikes_n6_average(:,j));
%     end
% end
% alpha_adj=0.05/45;
% nonsig=length(find(p>alpha_adj));
% 
% 
% % Simple normalization n=15 - multcompare
% [p_value,~,stats]  = kruskalwallis(sum_spikes_n6_average);
% [results,means] = multcompare(stats, 'Alpha',0.05,'CType','bonferroni');
% pvalues=results(:,6);
% alpha_adj=0.05/45;
% nonsig=length(find(pvalues>alpha_adj));
% sig=length(find(pvalues<=alpha_adj));
% 
% clear p_value results means pvalues nonsig sig stats
% 
% 
% % Simple normalization n=421 - wilcoxon
% sum_spikes_n_r=[];
% for i=1:15
%     sum_spikes_n_r=[sum_spikes_n_r;sum_spikes_n6(:,:,i)];
% end
% 
% clear aux
% aux=find(isnan(sum_spikes_n_r(:,1))==1);
% sum_spikes_n_r(aux,:)=[];
% 
% figure
% boxplot(sum_spikes_n_r)
% ylabel('Normalized activity during wave interval');
% xlabel('Phase Interval');
% 
% [p_value,~,stats] = kruskalwallis(sum_spikes_n_r);
% p = friedman(sum_spikes_n_r,1);
% 
% 
% count=0;
% for i=1:10
%     for j=i+1:10
%         count=count+1;
% 
%         [p(count)]=ranksum(sum_spikes_n_r(:,i),sum_spikes_n_r(:,j));
%     end
% end
% nonsig=length(find(p>(0.05/45)));
% 
% % Simple normalization n=421 - multcompare
% [results,means] = multcompare(stats, 'Alpha',0.05,'CType','bonferroni');
% pvalues=results(:,6);
% alpha_adj=0.05/45;
% nonsig=length(find(pvalues>alpha_adj));
% sig=length(find(pvalues<=alpha_adj));

% 