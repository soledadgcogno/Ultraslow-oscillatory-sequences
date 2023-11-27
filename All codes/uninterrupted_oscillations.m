%Uninterrupted oscillations

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

%% Oscillations


countf=0;
for w=1:length(waves)
    row_w=waves(w);
    disp(w)

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
    [N,T]=size(spikes_d);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Condition on having waves
    sf=7.73;
    sf_d=7.73/4;
    num_clus_discr=10;
    make_fig=1;
    dt=floor(big_table(row_w,8));
    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes_d);

    for i=2:size(table_u,1)
        ICI_ses(i-1)=(table_u(i,1)-table_u(i-1,2));
    end

    epoch=[];
    count=0;
    for i=1:size(ICI_ses,2)
        %disp(i);
        if ICI_ses(i)==1
            count=count+1;
            epoch(count,1)=i;
            for j=i+1:size(ICI_ses,2)
                if ICI_ses(j)==1
                else
                    epoch(count,2)=j;%includes the last cycle of the succession of uninterrupted cycles
                    break
                    %                     i=j+1;
                end
            end

        end
    end

    aux=find(epoch(:,2)==0);
    epoch(aux,2)=size(table_u,1);
    [a,r]=max(abs(epoch(:,2)-epoch(:,1)))

    max_num_uninterrupted_cycles(w)=epoch(r,2)+1-(epoch(r,1));
    max_time_uninterrupted_cycles(w)=((table_u(epoch(r,2),2) - table_u(epoch(r,1),1))/sf)/60; %in min

    cycle_length=(table_u(:,2)-table_u(:,1))/sf/60;
    ini=epoch(r,1);
    endi=epoch(r,2);
    max_time_uninterrupted_cycles_2(w)=sum(cycle_length(ini:endi))+((max_num_uninterrupted_cycles(w)-1)*0.129/60);


    clear X Y spikes spikes_d_s score mouse i FRp FR coeff cells_d T window score FR FRp epoch ICI_ses table_u spikes spikes_d ICI_ses_t ...
        aux

end

%% Figures

figure
hold on
h=plot(1:3,max_time_uninterrupted_cycles(1:3),'-o','markersize',10,'color','k','markerfacecolor','k','linewidth',2);
set(gca,'fontsize',16,'XColor','k','YColor','k');
ylabel({'Max length of uninterrupted';'oscillations (min)'});
xlabel('Session #')
axis([0.5 3.5 0 16])
xticks([1 2 3])
box off

figure
hold on
h=plot(1:3,max_time_uninterrupted_cycles(4:6),'-o','markersize',10,'color','k','markerfacecolor','k','linewidth',2);
set(gca,'fontsize',16,'XColor','k','YColor','k');
ylabel({'Max length of uninterrupted';'oscillations (min)'});
xlabel('Session #')
axis([0.5 3.5 0 15])
xticks([1 2 3])
box off

figure
hold on
h=plot(1:4,max_time_uninterrupted_cycles(7:10),'-o','markersize',10,'color','k','markerfacecolor','k','linewidth',2);
set(gca,'fontsize',16,'XColor','k','YColor','k');
ylabel({'Max length of uninterrupted';'oscillations (min)'});
xlabel('Session #')
axis([0.5 4.5 0 27])
xticks([1 2 3 4])
box off

figure
hold on
h=plot(1:5,max_time_uninterrupted_cycles(11:15),'-o','markersize',10,'color','k','markerfacecolor','k','linewidth',2);
set(gca,'fontsize',16,'XColor','k','YColor','k');
ylabel({'Max length of uninterrupted';'oscillations (min)'});
xlabel('Session #')
axis([0.5 5.5 0 6])
xticks([1 2 3 4 5])
box off
