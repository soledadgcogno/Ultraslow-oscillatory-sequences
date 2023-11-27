% First cell with table and wave sessions
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
    ws_prob=load([save_data_path ['WS_Prob_for more than 3 ensembles_dt66_',mice(m,:),'.mat']]);

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
                    big_table(count,11)=ws_prob.WS_stat.wave_score_prob(day,s); %WS - Ent
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

% Wave sessions


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


%% PI and frequency

sf=7.73;
downsampling_factor=4;
fs_120 = 7.73;
fs_30 = 30.95;
fs = fs_120/downsampling_factor;
maxlag=floor(560*fs_120);%240;%250 - 520;

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
    [N,T]=size(spikes_d);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Condition on having waves

    dt=floor(big_table(row_w,8));  
   
    num_clus_discr=10;
    make_fig=0;
    [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes_d);

    spikes_r=[]; %Reduced spike matrix; only contains wave epochs
    phase_r=[]; %Reduced phase; only contains wave epochs


    for i=1:N
        FRp(i,:)=full(fire_rate(spikes_d(i,:),1*dt,'g')); %smooth using as kernel the dt chosen for each session
    end

    [coefft,scoret,~] = pca(FRp');
    phase_f=(atan2(smooth(scoret(:,2),floor(1*dt)),smooth(scoret(:,1),floor(1*dt))));
    radius_f=sqrt(coefft(:,2).*coefft(:,2)+coefft(:,1).*coefft(:,1));
    for wa=1:size(table_u,1)
        spikes_r=[spikes_r,spikes_d(:,table_u(wa,1):table_u(wa,2))];
        phase_r=[phase_r;phase_f(table_u(wa,1):table_u(wa,2))];
    end

    spikes=spikes_d;
    phase=phase_r;
    spikes_d=[];
    spikes_d=spikes_r;
    T_d=size(spikes_r,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MVL
    
    for i=1:N
        p=phase(find(spikes_d(i,:)));
        mean_p(i)=circ_mean(p);
        MVL{w}(i) = circ_r(p);
        clear p H
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PI


    for i=1:size(spikes,1) %Loop on elements of mat
        calc=(spikes(i,:)); %spike train of row i

        for wa=1:size(table_u,1) %Spikes of each cell per wave
            spikes_per_wave(i,wa)=sum(calc(table_u(wa,1):table_u(wa,2)));
        end

        aux=spikes_per_wave(i,:);
        aux2=find((cumsum(sort(aux,'descend'))./sum(aux))>0.90,1);%Number of waves needed to account for 95% of the spikes
        number_events{w}(i)=sum(aux)/(T_d/sf);

        if isempty(aux2)==1
            aux2=0;
        end

        PI{w}(i)=aux2/size(table_u,1);
    end


    %%%%%%%%%%%%%%%%%%%%%%%% Frequency of single cells

    for i=1:size(spikes,1) %Loop on elements of mat

        [alfa(i,:),lags]=xcorr(spikes_r(i,:),spikes_r(i,:),maxlag,'coeff');
        [val(i,:),~]=zscore(alfa(i,:)); %Check whether I need to zscore
        signal=alfa(i,:);

        clear Powerspec2 Powerfreq2

        %Calculate spectogram using pwelch method
        [Powerspec2,Powerfreq2] = doPwelch(signal,fs_120,2*4096);
%         pwelch_fft(i,:)=Powerspec2;
%         pwelch_fft_section(i,:)=Powerspec2(5:100);
        [peak_freq(i),peak_psd(i),quality(i)]=check_peak_quality_3c_f(Powerspec2,Powerfreq2,0.04);
%         period_ps(i)=1/peak_freq(i);

        clear  Powerspec2 P2 P1  Y

    end

    % Frequency of the global oscillation
    full_length_oscillation=sum(table_u(:,2)-table_u(:,1))/sf;
    full_length_oscillation2=size(spikes_d,2)/sf;

    freq1(w)=size(table_u,1)/full_length_oscillation;
    freq2(w)=size(table_u,1)/full_length_oscillation2;

    rel_freq{w}=peak_freq/freq1(w);
    freq_sc{w}(1:length(peak_freq))=peak_freq;
    freq_s=abs(rel_freq{w}-1);

    [a,b]=sort(freq_s, 'ascend');

    % Analysis of fraction
    countf=0;
    for fraction=[0.05,0.1,0.2,0.3,0.4,0.5]
        countf=countf+1;
        n_cell=floor(fraction*length(rel_freq{w}));
        cells_sim=b(1:n_cell);
        cells_diff=b(end-n_cell+1:end);

        mvlsim(w,countf)=mean(MVL{w}(cells_sim));
        mvldiff(w,countf)=mean(MVL{w}(cells_diff));
        pisim(w,countf)=mean(PI{w}(cells_sim));
        pidiff(w,countf)=mean(PI{w}(cells_diff));

        %mdl = fitlm(number_events{w}(cells_sim),PI{w}(cells_sim));
        table{w,countf}(:,1)=number_events{w}(cells_sim);
        table{w,countf}(:,2)=PI{w}(cells_sim);


        clear cells_sim  cells_diff
    end

    clear cells cells_di cells_d coeff dff exclude1 exclude2 latent_pca MI_TB phase phase_d phase_di score signal_dff snr SNR spikes spikes_d spikes_d_s p 
    clear not_locked locked locking locking_sh locking_sh_mean locking_sh_99 locking_sh_1 MVL_sh p_sh sp_do phase phase_d mean_p mean_p_sh prob_phase_firing
    clear std_p std_p__sh var_p var_p__sh calc_di coefft FRp MI MI_b MI_withbias phase_down phase_f phase_r radius_f scoret spikes_r spk_do ...
        table_u Anat d_within d_across delta_tissue delta_tissue_vec ens Ens ense_n  sorting_w r_i r_j PR new_mat new_mat2 spikes_sorted ...
        PR_locked PR_not_locked index_sorting_locked sorting_w_not_locked sorting_w_locked delta_phase_mean index_sorting_not_locked...
        Sens_sh AUC_sh TC2_sh TC3_sh Arg_sens2_sh Sens AUC TC2 TC3 Arg_sens2 p_1 p_0 mean_fr calc cell_wave_part prob_k sp_do ind...
        pwelch_fft_section pwelch_fft quality val alfa lags peak_freq peak_psd quality signal spikes_per_wave ...
        Powerspec2 Powerfreq2 period_ps

end


%% Figure for example session


figure
scatter(rel_freq{8},PI{8},'filled','k');
ylabel('Participation index');
xlabel('Frequency relative to oscillation');
set(gca,'fontsize',16,'ycolor','k','xcolor','k');
axis([0 8 0 1])


cc=jet(4);
ev_d=discretize(number_events{8},[0,prctile(number_events{8},25),prctile(number_events{8},50),prctile(number_events{8},75),max(number_events{8})]);
figure
scatter(rel_freq{8},PI{8},80,cc(ev_d,:),'filled');
ylabel('Participation index');
xlabel('Frequency relative to oscillation');
set(gca,'fontsize',16,'ycolor','k','xcolor','k');
axis([0 8 0 1])
colormap jet(4)
co=colorbar('XTick',0:1,'XTickLabel',{'min','max'});
co.Label.String = 'Event rate';

%     
% figure
% scatter(peak_freq,MVL,'filled','k');
% ylabel('Participation index');
% xlabel('Frequency of single cells (Hz)');
% set(gca,'fontsize',16,'ycolor','k','xcolor','k');
% axis([0 0.06 0 1])

figure
scatter(rel_freq{8},MVL{8},'filled','k');
ylabel('Locking degree');
xlabel('Frequency relative to oscillation');
set(gca,'fontsize',16,'ycolor','k','xcolor','k');
axis([0 7 0 1])

freq_s=abs(rel_freq{8}-1);
[a,b]=sort(freq_s, 'ascend');

fraction=0.1;
n_cell=floor(fraction*length(rel_freq{8}));
cells_sim=b(1:n_cell);
cells_diff=b(end-n_cell+1:end);

figure
boxplot([MVL{8}(cells_sim)',MVL{8}(cells_diff)']);
ylabel('Locking degree');
xticklabels({'Rel. frequency ~ 1','Rel. frequency > 1'});
set(gca,'fontsize',16,'XColor','k','YColor','k');
axis([0.5 2.5 0 1]);
box off

[p,h,stat]=ranksum(MVL{8}(cells_sim)',MVL{8}(cells_diff)');

figure
boxplot([PI{8}(cells_sim)',PI{8}(cells_diff)']);
ylabel('Participation index');
xticklabels({'Rel. frequency ~ 1','Rel. frequency > 1'});
set(gca,'fontsize',16,'XColor','k','YColor','k');
axis([0.5 2.5 0 1]);
box off
[p,h,stat]=ranksum(PI{8}(cells_sim)',PI{8}(cells_diff)');


figure
scatter(number_events{8}(cells_sim),PI{8}(cells_sim),80,'k','filled');
ylabel('Participation index');
xlabel('Event rate (Hz)');
set(gca,'fontsize',16,'XColor','k','YColor','k');
axis([0 0.35 0 1])
box off
[p_lf,S_lf] = polyfit(number_events{8}(cells_sim),PI{8}(cells_sim),1);
[f,delta] = polyval(p_lf,number_events{8}(cells_sim),S_lf);
hold on
plot(number_events{8}(cells_sim),f,'-m','linewidth',2);
clear p_lf    S_lf f delta

mdl = fitlm(number_events{8}(cells_sim),PI{8}(cells_sim));
clear mdl

%% Figures all sessions

figure
boxplot([mvlsim(:,2),mvldiff(:,2)]);
ylabel('Locking degree');
xticklabels({'Rel. frequency ~ 1','Rel. frequency > 1'});
set(gca,'fontsize',16,'XColor','k','YColor','k');
axis([0.5 2.5 0 1]);
box off
title('All sessions');
[p,h,stats]=ranksum(mvlsim(:,2),mvldiff(:,2));

figure
boxplot([pisim(:,2),pidiff(:,2)]);
ylabel('Participation index');
xticklabels({'Rel. frequency ~ 1','Rel. frequency > 1'});
set(gca,'fontsize',16,'XColor','k','YColor','k');
axis([0.5 2.5 0 1]);
box off
title('All sessions');
[p,h,stats]=ranksum(pisim(:,2),pidiff(:,2));

for i=1:size(mvlsim,2)
    [p_mvl(i),h_mvl(i)]=ranksum(mvlsim(:,i),mvldiff(:,i));
    [p_pi(i),h_pi(i)]=ranksum(pisim(:,i),pidiff(:,i));
end

figure
plot([5,10,20,30,40,50],p_mvl,'k-*','linewidth',2);
ylabel('p-value - Locking degree');
xlabel('Fraction of cells %');
axis([0 55 0 0.004]);
set(gca,'fontsize',16,'XColor','k','YColor','k');
box off

figure
plot([5,10,20,30,40,50],p_pi,'k-*','linewidth',2);
ylabel('p-value - Participation index');
xlabel('Fraction of cells %');
hold on
yline(0.05,'--');
axis([0 55 0 1]);
set(gca,'fontsize',16,'XColor','k','YColor','k');
box off

for j=1:6
    table_all=[];
    for i=1:15
        table_all=[table_all;table{i,j}];
    end

    if j==2
        figure
        scatter(table_all(:,1),table_all(:,2),80,'k','filled');
        ylabel('Participation index');
        xlabel('Event rate (Hz)');
        set(gca,'fontsize',16,'XColor','k','YColor','k');
        axis([0 inf 0 1])
        box off
        [p_lf,S_lf] = polyfit(table_all(:,1),table_all(:,2),1);
        [f,delta] = polyval(p_lf,table_all(:,1),S_lf);
        hold on
        plot(table_all(:,1),f,'--m','linewidth',2);

    end

    mdl = fitlm(table_all(:,1),table_all(:,2));

    R(j)=mdl.Rsquared.Ordinary;
    p_l(j)=coefTest(mdl);
    [rho_s(j),p_s(j)] = corr(table_all(:,1),table_all(:,2),'Type','Spearman');

    clear table_all mdl
end


figure
plot([5,10,20,30,40,50],p_l,'k-*','linewidth',2);
ylabel('p-value - Line fit');
xlabel('Fraction of cells %');
axis([0 55 0 0.001]);
set(gca,'fontsize',16);
box off


