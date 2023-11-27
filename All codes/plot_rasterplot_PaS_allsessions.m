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
    load([save_data_path ['WS_Osc_14_',mice(m,:),'.mat']]);
    ws_ent=load([save_data_path ['WS_Entropy_',mice(m,:),'.mat']]);
 
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

%%

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
for w=5:length(repeated_days)
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
    ws=big_table(row_w,6);
    ws_ent=big_table(row_w,11);
 
    clus=10;
    disc_phase=10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load files
    
    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    file_name_snr=[dpath ['SNR_DFF_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    
    if (isfield(dates,'folder_name')==1)
        day_index=find (dates.actual_day==day);
        file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day_index),'_MUnit',num2str(munit),'.mat']];
    else
        file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
    end
    
     if(exist(file_name_spk)==2)
        
        %     load(file_name_dff,'-mat'); %DFF
        load(file_name_spk,'-mat'); %Spike times
        spikes_d=full(spikes_d_s);
        [N,T]=size(spikes_d);        
        [sorting_ascend,sortind_descend,~]=get_sorting(spikes_d);
        
%         fig=figure;
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
%         hold on
%         spy(flip(spikes_d(sorting_descend,:)),'k')
%         pbaspect([25 2.5 1]);
%         xlabel('Time');
%         ylabel('Cells');
%         title([mouse,'    -    Day=',num2str(day),'    -    Session=',num2str(s),'    -    MUnit=',num2str(munit),'    -    Wave score = ',num2str(ws)...
%             ,'    -    Wave score Ent= ',num2str(ws_ent)]);
%         set(gca,'fontsize',18);
        
        
        fig=figure;
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
        hold on
        for i=1:size(spikes_d,1)
            %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_descend(i),:),5,'k','filled')
            scatter((1:size(spikes_d,2))./8,i*spikes_d(sortind_descend(i),:),5,'k','filled')            
            alpha 0.3
        end
        axis([-inf inf 1 inf]);
        %title([mouse,' Day',num2str(day)]);
        xlabel('Time (s)');
        ylabel('Neurons #');
        title([mouse,'    -    Day=',num2str(day),'    -    Session=',num2str(s),'    -    MUnit=',num2str(munit),'    -    Wave score = ',num2str(ws)...
            ,'    -    Waveness= ',num2str(ws_ent)]);
        set(gca,'fontsize',18);
        yticks([100 300])             

%         C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July 2021\Raster plots
        
        saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July 2021\Raster plots\PaS\',mouse,'_Day',num2str(day),'_Session',num2str(s),'.svg'));
        saveas(gcf,strcat('C:\Users\xscogno\Dropbox\Paper Waves\Figures Waves Paper\June-July 2021\Raster plots\PaS\',mouse,'_Day',num2str(day),'_Session',num2str(s),'.fig'));

        close all
     end
     
     clear spikes_d_s spikes_d sorting_descend sorting_ascend sorting_0 N mouse cells_d
end

%% Compare wave scores

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
    ws_binary(w)=big_table(row_w,6);
    ws_entropy(w)=big_table(row_w,11);
    ML_pos(w)=big_table(row_w,10);
    
end
% 
figure
scatter(ws_binary,ws_entropy,50,'ko','filled');
alpha 0.5
axis([-0.5 1.5 0.6 2.4])
xlabel('Wave score');
ylabel('Waveness');
set(gca,'fontsize',16)
xticks([0 1])
% 
% % figure
% % scatter(ML_pos,ws_binary,50,'ko','filled');
% % alpha 0.5
% % % axis([-0.5 1.5 0.6 1.6])
% % xlabel('Binary wave score');
% % ylabel('Continuous wave score');
% % set(gca,'fontsize',16)
% % 
% % figure
% % scatter(ML_pos,ws_entropy,50,'ko','filled');
% % alpha 0.5
% % % axis([-0.5 1.5 0.6 1.6])
% % xlabel('Binary wave score');
% % ylabel('Continuous wave score');
% % set(gca,'fontsize',16)
