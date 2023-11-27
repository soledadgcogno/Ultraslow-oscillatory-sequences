%Generates the spike matrix
clear all
dini_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
snr_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\SNR_H\SNR_flavio_new_L8\';
spath = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

mouse='92227';  %Name of the animal

thr_SNR=4;
res=1.5;
count=0;
new_bin=1;
f=1;
% for m=1:4%1:2%:mice_number

load([dini_path,strcat('recording_dates_',mouse,'.mat')]);

days_num=dates.daysnum;
FOV_n=1;

for day=1:days_num
    
    disp(strcat('Day ',num2str(day),' of ',num2str(days_num)));
    ses_num=dates.sesnum(day);
    
    for s=1:ses_num
        munit=dates.ses{day}(s);
        %             if (ses==dates.FOV{1,f}(day))
        if (length(num2str(dates.days(day)))==8)
            
            dpath=strcat('Z:\2P\V1\Spring 2020\',mouse,'\',num2str(dates.days(day)),'\MUnit_',num2str(munit),'\',dates.folder_name{s,day},...
                '\flavio-develop\plane0\');            
            file_name=strcat(dpath,'Fall.mat');
            file_name_snr= [snr_path strcat(mouse,'\Flavio_2P_V1_Spring 2020_',mouse,'_',num2str(dates.days(day)),'_Munit_',num2str(munit),'.csv')];

        else
            
             dpath=strcat('Z:\2P\V1\Spring 2020\',mouse,'\0',num2str(dates.days(day)),'\MUnit_',num2str(munit),'\',dates.folder_name{s,day},...
                '\flavio-develop\plane0\');            
            file_name=strcat(dpath,'Fall.mat');
            file_name_snr= [snr_path strcat(mouse,'\Flavio_2P_V1_Spring 2020_',mouse,'_0',num2str(dates.days(day)),'_Munit_',num2str(munit),'.csv')];
            
%             disp('Error');
%             break
%             dpath=strcat('Z:\2P\',mouse,'\0',num2str(dates.days(day)),'\MUnit_',num2str(munit),'\suite2p\plane0\');            
%             file_name=strcat(dpath,'Fall.mat');
%             file_name_snr= [snr_path strcat(mouse,'\Flavio_2P_',mouse,'_0',num2str(dates.days(day)),'_Munit_',num2str(munit),'.csv')];
        end
        
        
        if(exist(file_name_snr)==2)
            
            dat=load(file_name,'-mat');
            D = importdata(file_name_snr);
            a=size(D.data,1);
            b=size(D.textdata,1)-1;
            
            if(a==b)
                
               iscell_fall=find(dat.iscell(:,1)>0);
 
                iscell_aux=(D.textdata(2:end,3));     %is cell
                N=size(iscell_aux);
                for i=1:N
                    iscell(i,1) = str2num(cell2mat(iscell_aux(i))) + 1;
                end
                
                sp_nfil=dat.spks;            %spikes (not filtered)
                sp_fil=sp_nfil(iscell,:);   %spikes fitlered
                
                
                [N,T] = size(sp_fil);
                SNR=D.data(:,1);
                
                %Check this part!
                new_bin_num = floor(size(sp_fil,2)/new_bin);
                for i=1:new_bin_num
                    sp_fil_do(:,i)=sum(sp_fil(:,(i-1)*new_bin+1:i*new_bin),2)/new_bin;
                end
                                
                spikes=[];
                count=1;
                for i=1:N
                    if(SNR(i)>thr_SNR && isnan(SNR(i))==0 )
                        
                        threshold(i) = mean(nonzeros(sp_fil_do(i,:))) + res*std(nonzeros(sp_fil_do(i,:))) ;
                        
                        sp_fil_do(i,find(sp_fil_do(i,:)<=threshold(i)))=0;    %Binarize spikes
                        sp_fil_do(i,find(sp_fil_do(i,:)>threshold(i)))=1;     %Binarize spikes
                        spikes(count,:)=sp_fil_do(i,:);
                        
                        cells_d(count)=iscell(i);
                        count=count+1;
                    end
                end
                
                spikes_d_s=sparse(double(spikes));
                
                if (isempty(spikes_d_s)==1)
                    
                else
%                     save([spath,'spikes_120ms_Do_THR',num2str(res),'_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit)],'spikes_d_s','cells_d');
                                    save([spath,'spikes_30ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit)],'spikes_d_s','cells_d');
                    
                end
                
                clear snr_t y_ signal spikes spikes_d spikes_d_s  sp_fil_do SNR T_value cells cells_aux cells_file cells_d spikes_d_s cells_d ...
                    sp_fil_do D
                
                clear prob_data prob sp_fil iscell dat tf_ifd mat_simil U V S Time_edos Cells_edos activity_sh activity_sh_2 one neurons...
                    activity_data H_index H_indexb S_indexp scut_sv corte repetitions stat edos_pks_num edos_temp1...
                    frames_peaks mean_act tot_neurons thr_Simil fac_cut_d min_rep cells_ens corte scut_svd prob_rep repetitions N_sh...
                    cells_per_ens C ia ib neurons empty_cells_2 tf_idf Rep_edos Hamming_edos Norm_act mat_cell sim sim_o count frames_ens cd...
                    sis_query Cellsc Cellsb h_temp states_corr sequences_corr prob_data threshold spikes neuropilCoefficient FcellNeu Fcell N T D dat...
                    SNR cells_d threshold spikes_d_s dat sp_fil_do sp_nfil sp_fil iscell iscell_aux sp_nfil dat spikes_d_s cells_d SNR iscell D;
                
            else
                disp([day , s]);
            end
        end
        
        %             end
        clear snr_t y_ signal spikes spikes_d spikes_d_s  sp_fil_do SNR T_value cells cells_aux cells_file cells_d spikes_d_s cells_d ...
            sp_fil_do D
        
        clear prob_data prob sp_fil iscell dat tf_ifd mat_simil U V S Time_edos Cells_edos activity_sh activity_sh_2 one neurons...
            activity_data H_index H_indexb S_indexp scut_sv corte repetitions stat edos_pks_num edos_temp1...
            frames_peaks mean_act tot_neurons thr_Simil fac_cut_d min_rep cells_ens corte scut_svd prob_rep repetitions N_sh...
            cells_per_ens C ia ib neurons empty_cells_2 tf_idf Rep_edos Hamming_edos Norm_act mat_cell sim sim_o count frames_ens cd...
            sis_query Cellsc Cellsb h_temp states_corr sequences_corr prob_data threshold spikes neuropilCoefficient FcellNeu Fcell N T D dat...
            SNR cells_d threshold spikes_d_s dat sp_fil_do sp_nfil sp_fil iscell iscell_aux sp_nfil dat spikes_d_s cells_d SNR iscell;
    end
end
    
    %
% end
% end

