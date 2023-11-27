%% Computes sequence score by calculating the probability and asking for it to be larger than 3 ensembles
%Without significance
clear all
close all

rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

clusters=10;
mice_number=12;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];
ml=[2.64,3.7,2.7,2.64,3.2,3.3,2.52,3.3,3.5,-10,-10,-10]; %ML coordinates

index_c=0;
index_t=1;
clus=10;
N_sh=500;
sf=7.73;

for m=1:mice_number
    disp(m)
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

    for day=1:dates.daysnum
        for s=1:dates.sesnum(day)

            munit=dates.ses{day}(s);
            file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

            if exist(file_name)==2
                load(file_name,'-mat');
                spikes=full(spikes_d_s);
                [N,T]=size(spikes);

                dt=WS_stat.dt(day,s);
                if isinteger(dt)
                else
                    dt=floor(dt);
                end

                if dt==0
                    dt=66;
                end

                if N>150 && T>7000

                    if (WS_stat.WS(day,s)==1)


                        num_clus_discr=10;
                        make_fig=0;
                        [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);

                        spikes_w=[];
                        for i=1:size(table_u,1)
                            spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
                        end

                        [prob_data,~] = comp_wave_prob_f(spikes_w,floor(dt),clus);
                        Prob=sum(prob_data(3:end));

                      

                        optimal_dt(day,s)=dt;
                        seq_score_prob(day,s)=Prob;

                    else

                        prob_data = comp_wave_prob_f(spikes,dt,clus);
                        Prob=sum(prob_data(3:end));

                      

                        optimal_dt(day,s)=dt;
                        seq_score_prob(day,s)=Prob;


                    end
                else
                    optimal_dt(day,s)=NaN;
                    seq_score_prob(day,s)=NaN;
                  
                end
                clear spikes prob_data spikes_d_s cells_d spikes_w table_u N T Prob_sh prob_data_sh mat_sh
            end
        end
    end

    WS_stat.optimal_dt=optimal_dt;
    WS_stat.seq_score_prob=seq_score_prob;


    save([save_data_path,'WS_Prob_for more than 3 ensembles_dt66_sf7p73_',mice(m,:)],'WS_stat');
    clear optimal_dt wave_score_prob dates mouse thr Entr_up WS_stat wave_score_wave WS_stat
end

%% With significance
clear all
close all

rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

clusters=10;
mice_number=12;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];
ml=[2.64,3.7,2.7,2.64,3.2,3.3,2.52,3.3,3.5,-10,-10,-10]; %ML coordinates

index_c=0;
index_t=1;
clus=10;
sf=7.73;
N_sh=500;

for m=[2,5,6,9]
    disp(m)
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

    for day=1:dates.daysnum
%         if day<15
%             N_sh=2;
%         else
%             N_sh=500;
%         end
        for s=1:dates.sesnum(day)

            munit=dates.ses{day}(s);
            file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];

            if exist(file_name)==2
                load(file_name,'-mat');
                spikes=full(spikes_d_s);
                [N,T]=size(spikes);

                dt=WS_stat.dt(day,s);
                if isinteger(dt)
                else
                    dt=floor(dt);
                end

                if dt==0
                    dt=66;
                end

                if N>150 && T>7000

                    if (WS_stat.WS(day,s)==1)


                        num_clus_discr=10;
                        make_fig=0;
                        [table_u,N,T]=identify_waves_latestversion_6_f(mouse,day,num_clus_discr,dt,make_fig,spikes);

                        spikes_w=[];
                        for i=1:size(table_u,1)
                            spikes_w = horzcat(spikes_w,spikes(:,table_u(i,1):table_u(i,2)));
                        end

                        [prob_data,~] = comp_wave_prob_f(spikes_w,floor(dt),clus);
                        Prob=sum(prob_data(3:end));

                        for i=1:N_sh
                            mat_sh=shuffle(spikes_w')';
                            [prob_data_sh,~] =comp_wave_prob_f(mat_sh,floor(dt),clus);
                            Prob_sh(i)=sum(prob_data_sh(3:end));
                            clear mat_sh
                        end
                        thr_w=prctile(Prob_sh',99);

                        if(Prob>thr_w)
                            seq_score_prob_sig(day,s)=1;
                        else
                            seq_score_prob_sig(day,s)=0;
                        end

                        optimal_dt(day,s)=dt;
                        seq_score_prob(day,s)=Prob;

                    else

                        prob_data = comp_wave_prob_f(spikes,dt,clus);
                        Prob=sum(prob_data(3:end));

                        for i=1:N_sh
                            mat_sh=shuffle(spikes')';
                            prob_data_sh=comp_wave_prob_f(mat_sh,floor(dt),clus);
                            Prob_sh(i)=sum(prob_data_sh(3:end));
                            clear mat_sh
                        end
                        thr_w=prctile(Prob_sh',99);

                        if(Prob>thr_w)
                            seq_score_prob_sig(day,s)=1;
                        else
                            seq_score_prob_sig(day,s)=0;
                        end


                        optimal_dt(day,s)=dt;
                        seq_score_prob(day,s)=Prob;


                    end
                else
  
                    optimal_dt(day,s)=NaN;
                    seq_score_prob(day,s)=NaN;
                    seq_score_prob_sig(day,s)=NaN;

                end
                clear spikes prob_data spikes_d_s cells_d spikes_w table_u N T Prob_sh prob_data_sh mat_sh
            end
        end
    end

    WS_stat.optimal_dt=optimal_dt;
    WS_stat.seq_score_prob=seq_score_prob;
    WS_stat.seq_score_prob_sig=seq_score_prob_sig;

    save([save_data_path,'WS_Prob_for more than 3 ensembles_dt66_sf7p73_',mice(m,:),'with_significance'],'WS_stat');
    clear optimal_dt wave_score_prob dates mouse thr Entr_up WS_stat wave_score_wave WS_stat
end


