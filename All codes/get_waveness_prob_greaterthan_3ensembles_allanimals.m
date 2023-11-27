%% Computes waveness using the probability and asking for it to be larger than 3 ensembles

clear all
close all

rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

clusters=10;
mice_number=12;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];
ml=[2.64,3.7,2.7,2.64,3.2,3.3,2.52,3.3,3.5,-10,-10,-10]; %ML coordinates
% times=[24,48,72,96,120];

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
                        
                        for i=1:N_sh
                            mat_sh=shuffle(spikes_w')';
                            [prob_data_sh,~] =comp_wave_prob_f(mat_sh,floor(dt),clus);
                            Prob_sh(i)=sum(prob_data_sh(3:end));
                            clear mat_sh
                        end
                        thr_w=prctile(Prob_sh',99);
                        
                        if(Prob>thr_w)
                            wave_score_prob_sig(day,s)=1;
                        else
                            wave_score_prob_sig(day,s)=0;
                        end
                        
                        optimal_dt(day,s)=dt;
                        wave_score_prob(day,s)=Prob;
                        wave_score_wave(day,s)=1;
                        
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
                            wave_score_prob_sig(day,s)=1;
                        else
                            wave_score_prob_sig(day,s)=0;
                        end
                        
                        
                        optimal_dt(day,s)=dt;
                        wave_score_prob(day,s)=Prob;
                        wave_score_wave(day,s)=0;
                        
                        
                    end
                else
                    optimal_dt(day,s)=NaN;
                    wave_score_prob(day,s)=NaN;
                    wave_score_wave(day,s)=NaN;
                    wave_score_prob_sig(day,s)=NaN;
                    
                end
                clear spikes prob_data spikes_d_s cells_d spikes_w table_u N T Prob_sh prob_data_sh mat_sh
            end
        end
    end
    
    WS_stat.optimal_dt=optimal_dt;
    WS_stat.wave_score_prob=wave_score_prob;
    WS_stat.wave_score_wave=wave_score_wave;
    WS_stat.wave_score_prob_sig=wave_score_prob_sig;
    
    save([save_data_path,'WS_Prob_for more than 3 ensembles_dt66_sf7p73_',mice(m,:),'with_significance'],'WS_stat');
    clear optimal_dt wave_score_prob dates mouse thr Entr_up WS_stat wave_score_wave WS_stat
end






%% Wave score First half - VS Second Half
clear all
close all

rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

clusters=10;
mice_number=14;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L10M1';'L10M2';'L11M1';'L11M2';'L11M3';'L12M1';'L12M2';'L12M3'];

times=[24,48,72,96];

index_c=0;
index_t=0;

for m=1:mice_number
        
if mice(m,2)=='0'
    mouse=[mice(m,1),mice(m,3:5)];
    mouse_l=mouse(2);
    mouse_a=mouse(4);
else
    mouse=mice(m,:);
    mouse_l=mouse(2:3);
    mouse_a=mouse(5);
end

    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    
    for day=1:dates.daysnum
        
        for s=1:dates.sesnum(day)
            disp(s)
            munit=dates.ses{day}(s);
            file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
                                  
            if exist(file_name)==2
                load(file_name,'-mat');
                spikes=full(spikes_d_s);
                [N,T]=size(spikes);
                
                if N>150 && T>7000
                    for clus=10%[5,10,15]
                        index_t=0;
                        for dt=times    %Loop of different temporal scales
                            index_t=index_t+1;
                            prob_data = comp_wave_prob(spikes,dt,clus);
                            Entr_up(index_t,day,s)=nansum(-prob_data.*log2(prob_data));
                        end
                    end
                    
                    [aup,bup]=max(Entr_up(:,day,s));                    
                    optimal_dt(day,s)=times(bup);
                    wave_score_ent(day,s)=aup;
                    thr(day,s)=2;%th;
                    
                    T_h=floor(T/2);
                    
                    mat1=spikes(:,1:T_h);
                    mat2=spikes(:,T_h+1:T);
                    
                    prob_data_1=comp_wave_prob_f(mat1,times(bup),clus);
                    prob_data_2=comp_wave_prob_f(mat2,times(bup),clus);
                    
                    WS1(day,s)=nansum(-prob_data_1.*log2(prob_data_1));
                    WS2(day,s)=nansum(-prob_data_2.*log2(prob_data_2));

                else
                    Entr_up(1:4,day,s)=NaN;                    
                    optimal_dt(day,s)=NaN;
                    wave_score_ent(day,s)=NaN;
                    thr(day,s)=NaN;%th;
                    WS1(day,s)=NaN;
                    WS2(day,s)=NaN;
                end
                clear spikes prob_data spikes_d_s cells_d prob_data_1 prob_data_2 mat1 mat2
            end
        end
    end
    
    WS_stat.optimal_dt=optimal_dt;
    WS_stat.wave_score_ent=wave_score_ent;
    WS_stat.Entr_up=Entr_up;
    WS_stat.th=thr;
    WS_stat.half1=WS1;
    WS_stat.half2=WS2;
 
    save([save_data_path,'WS2_Ent_',mice(m,:)],'WS_stat');
    clear optimal_dt wave_score_ent dates mouse thr Entr_up WS_stat WS1 WS2
end


%% Wave score First half - VS Second Half + Downsampling of cells
clear all
close all

rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

clusters=10;
mice_number=14;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L10M1';'L10M2';'L11M1';'L11M2';'L11M3';'L12M1';'L12M2';'L12M3'];

times=[24,48,72,96];

index_c=0;
index_t=0;
ds_fraction=[0.5,0.7,0.9];
N_sh=50;

for m=1:8%mice_number
        disp(m);
if mice(m,2)=='0'
    mouse=[mice(m,1),mice(m,3:5)];
    mouse_l=mouse(2);
    mouse_a=mouse(4);
else
    mouse=mice(m,:);
    mouse_l=mouse(2:3);
    mouse_a=mouse(5);
end

    load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
    
    for day=1:dates.daysnum
        
        for s=1:dates.sesnum(day)
%             disp(s)
            munit=dates.ses{day}(s);
            file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
                                  
            if exist(file_name)==2
                load(file_name,'-mat');
                spikes=full(spikes_d_s);
                [N,T]=size(spikes);
                
                if N>150 && T>7000
                    for clus=10%[5,10,15]
                        index_t=0;
                        for dt=times    %Loop of different temporal scales
                            index_t=index_t+1;
                            prob_data = comp_wave_prob(spikes,dt,clus);
                            Entr_up(index_t,day,s)=nansum(-prob_data.*log2(prob_data));
                        end
                    end                    
                    [aup,bup]=max(Entr_up(:,day,s));
                    optimal_dt(day,s)=times(bup);
                    wave_score_ent(day,s)=aup;
                    thr(day,s)=2;%th;
                    
                    index=0;
                    for ds_c=ds_fraction
                        index=index+1;
                        for sh=1:N_sh
                             cell_id=randperm(N,ceil(ds_c*N));
                             prob_data_ds=comp_wave_prob(spikes(cell_id,:),times(bup),clus);
                             WS_ds_cells(sh,index) = nansum(-prob_data_ds.*log2(prob_data_ds));
                             clear cell_id
                        end
                    end                    
                    WD_ds{day,s}=WS_ds_cells;            
                    
                    T_h=floor(T/2);                    
                    mat1=spikes(:,1:T_h);
                    mat2=spikes(:,T_h+1:T);       
                    prob_data_1=comp_wave_prob(mat1,times(bup),clus);                    
                    prob_data_2=comp_wave_prob(mat2,times(bup),clus);                    
                    WS1(day,s)=nansum(-prob_data_1.*log2(prob_data_1));
                    WS2(day,s)=nansum(-prob_data_2.*log2(prob_data_2));

                else
                    Entr_up(1:4,day,s)=NaN;                    
                    optimal_dt(day,s)=NaN;
                    wave_score_ent(day,s)=NaN;
                    thr(day,s)=NaN;%th;
                    WS1(day,s)=NaN;
                    WS2(day,s)=NaN;
                    WD_ds{day,s}=[];
                end
                clear spikes prob_data spikes_d_s cells_d WS_ds_cells
            end
        end
    end
    
    WS_stat.optimal_dt=optimal_dt;
    WS_stat.wave_score_ent=wave_score_ent;
    WS_stat.Entr_up=Entr_up;
    WS_stat.th=thr;
    WS_stat.half1=WS1;
    WS_stat.half2=WS2;
    WS_stat.ds=WD_ds;
    WS_stat.ds_fraction=ds_fraction;
 
    save([save_data_path,'WS3_Ent_',mice(m,:)],'WS_stat');
    clear optimal_dt wave_score_ent dates mouse thr Entr_up WS_stat WD_ds WS_ds_cells WD_ds WS_stat_ds
end


%% KL wave score 

clear all
close all
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
save_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

clusters=10;
mice_number=12;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];

times=[24,48,72,96,120];

index_c=0;
index_t=0;
N_sh=200;

for m=8:mice_number
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
    KL=nan(22,N_sh,4);
    
    for day=1:dates.daysnum
        for s=1:dates.sesnum(day)
            munit=dates.ses{day}(s);
            file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
            
            if exist(file_name)==2
                load(file_name,'-mat');
                spikes=full(spikes_d_s);
                [N,T]=size(spikes);
                
                if N>150 && T>7000
                    for clus=10%[5,10,15]
                        index_t=0;
                        for dt=times    %Loop of different temporal scales
                            index_t=index_t+1;
                            prob_data = comp_wave_prob(spikes,dt,clus);
                            Entr_up(index_t,day,s)=nansum(-prob_data.*log2(prob_data));
                        end
                    end
                    
                    [aup,bup]=max(Entr_up(:,day,s));
                    optimal_dt(day,s)=times(bup);
                    wave_score_ent(day,s)=aup;
                    
                    for dt=times(bup)
                        for sh=1:N_sh
                            mat_sh=shuffle(spikes')';
                            prob_data_sh = comp_wave_prob(mat_sh,dt,clus);
                            alfa=nonzeros(prob_data);
                            ind=find(ismember(prob_data,alfa));
                            KL(day,sh,s) = kldiv(prob_data_sh(ind),prob_data(ind));
                            AUC(day,sh,s) = sum(prob_data - prob_data_sh);
                            AUC_3(day,sh,s) = sum(prob_data(3:end) - prob_data_sh(3:end));

                            clear prob_data_sh alfa ind
                        end
                    end
                    
                    KL_99prc(day,s)= prctile( KL(day,:,s),99);
                    AUC_99prc(day,s)= prctile( AUC(day,:,s),99);
                    AUC_3_99prc(day,s)= prctile( AUC_3(day,:,s),99);

                    
                else
                    Entr_up(1:4,day,s)=NaN;
                    optimal_dt(day,s)=NaN;
                    wave_score_ent(day,s)=NaN;
                    KL_99prc(day,s)=NaN;
                    AUC_99prc(day,s)= NaN;
                    AUC_3_99prc(day,s)= NaN;

                end
                clear spikes prob_data spikes_d_s cells_d
            end
        end
    end
    
    WS_stat.optimal_dt = optimal_dt;
    WS_stat.wave_score_ent = wave_score_ent;
    WS_stat.Entr_up = Entr_up;
    WS_stat.KL = KL;
    WS_stat.AUC = AUC;
    WS_stat.AUC_3 = AUC_3;
    WS_stat.KL_prctile= KL_99prc;
    WS_stat.AUC_prctile= AUC_99prc;
    WS_stat.AUC_3_prctile= AUC_3_99prc;
    
    save([save_data_path,'WS_KL_Ent4_',mice(m,:)],'WS_stat');
    
    clear optimal_dt wave_score_ent dates mouse thr Entr_up WS_stat KL KL_99prc AUC AUC_3  AUC_3_99prc AUC_99prc
end
