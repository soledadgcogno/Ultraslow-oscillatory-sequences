clear all
close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

% clusters=53;
clusters=10;

mice_number=14;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L10M1';'L10M2';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];

for m=6
    
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
    
    for day=17
        
        for s=1:dates.sesnum(day)
            disp(s)
            munit=dates.ses{day}(s);
            file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
            
            load(file_name,'-mat');
            spikes=full(spikes_d_s);
            [N,T]=size(spikes);
            
            [~,sorting,~]=get_sorting_smoothed(spikes,117);

            downsampling_factor=8;
            mat=spikes(sorting,:);
            
            %Sorted Non-Threshold
%             FRp = spikes_downsample(mat,clusters,downsampling_factor);            
%             for i=1:clusters     
%                 FRp(i,:)=FRp(i,:)./max(FRp(i,:));                
%             end            
%             
%             figure
%             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
%             imagesc(FRp)
%             yticks([1 floor(clusters/2)+1 clusters])
% %             yticklabels({'50','25','1'});
%             yticklabels({'10','5','1'});
%             ylabel('Ensembles');
%             xlabel('Time (s)');
%             set(gca,'fontsize',16)
%             colormap jet
%             co=colorbar;
%             co.Label.String='Normalized Firing Rate';
            
            %Sorted Threshold

            clear FRp
            FRp = spikes_downsample(mat,clusters,downsampling_factor);
            for i=1:clusters
                FRp(i,:)=FRp(i,:)./max(FRp(i,:));
            end
            
            for i=1:clusters
                thr=mean(FRp(i,:))+std(FRp(i,:));
                FRp(i,find(FRp(i,:)<0.7*thr))=0;
            end
            
            figure
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
            imagesc(FRp)
            yticks([1 floor(clusters/2)+1 clusters])
%             yticklabels({'50','25','1'});
            yticklabels({'10','5','1'});
            ylabel('Ensembles');
            xlabel('Time (s)');
            set(gca,'fontsize',16)
            colormap jet
            co=colorbar;
            co.Label.String='Normalized Firing Rate';

            %Sorted Threshold- keeping the one with the maximum firing rate

            clear FRp
            FRp = spikes_downsample(mat,clusters,downsampling_factor);
            for i=1:clusters
                FRp(i,:)=FRp(i,:)./max(FRp(i,:));
            end
            
            for i=1:clusters
                thr=mean(FRp(i,:))+std(FRp(i,:));
                FRp(i,find(FRp(i,:)<0.7*thr))=0;
            end

            [~,signal]=max(FRp);

            FRmax=zeros(10,size(FRp,2));
            for i=1:size(FRp,2)
                FRmax(signal(i),i)=1;
            end
           
            figure
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
            imagesc(FRmax)
            yticks([1 floor(clusters/2)+1 clusters])
%             yticklabels({'50','25','1'});
            yticklabels({'10','5','1'});
            ylabel('Ensemble #');
            xlabel('Time (s)');
            set(gca,'fontsize',16,'ycolor','k','xcolor','k');
%             colormap jet
%             co=colorbar;
%             co.Label.String='Normalized Firing Rate';



            
            
            
            %Unsorted No-Threshold
% 
%                FRp = spikes_downsample(spikes,clusters,downsampling_factor);            
%             for i=1:clusters     
%                 FRp(i,:)=FRp(i,:)./max(FRp(i,:));                
%             end            
%             
%             figure
%             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
%             imagesc(FRp)
%             yticks([1 5 10])
%             yticklabels({'10','5','1'});
%             ylabel('Ensembles');
%             xlabel('Time (s)');
%             set(gca,'fontsize',16)
%             
%             
%             %Unsorted Threshold
% 
%             clear FRp
%             FRp = spikes_downsample(spikes,clusters,downsampling_factor);
%             
%             for i=1:clusters
%                 thr=mean(FRp(i,:))+std(FRp(i,:));
%                 FRp(i,find(FRp(i,:)<thr))=0;
%                 %                 FRp(i,find(FRp(i,:)>thr))=FRp();
%                 FRp(i,:)=FRp(i,:)./max(FRp(i,:));                
%             end
%             
%             figure
%             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
%             imagesc(FRp)
%             yticks([1 5 10])
%             yticklabels({'10','5','1'});
%             ylabel('Ensembles');
%             xlabel('Time (s)');
%             set(gca,'fontsize',16)
%             
%             
%             %Circular shuffling No-Threshold
%             spikes_sh=circular_shuffling(spikes);
%             [coeffsh]=pca(spikes_sh');
%             anglesh=atan2(coeffsh(:,2),coeffsh(:,1));
%             [a1,sorting_sh_pca]=sort(anglesh,'ascend');
%             mat=spikes_sh(sorting_sh_pca,:);
%             
%             
%             FRp = spikes_downsample(mat,clusters,downsampling_factor);
%             for i=1:clusters
%                 FRp(i,:)=FRp(i,:)./max(FRp(i,:));
%             end
%             
%             figure
%             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
%             imagesc(FRp)
%             yticks([1 5 10])
%             yticklabels({'10','5','1'});
%             ylabel('Ensemble #');
%             xlabel('Time (s)');
%             set(gca,'fontsize',20)
%             
%             
%             %Circular shuffling Threshold
% 
%             clear FRp
%             FRp = spikes_downsample(spikes_sh(sorting_sh_pca,:),clusters,downsampling_factor);
%             
%             for i=1:clusters
%                 thr=mean(FRp(i,:))+std(FRp(i,:));
%                 FRp(i,find(FRp(i,:)<thr))=0;
%                 %                 FRp(i,find(FRp(i,:)>thr))=FRp();
%                 FRp(i,:)=FRp(i,:)./max(FRp(i,:));                
%             end
%             
%             figure
%             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
%             imagesc(FRp)
%             yticks([1 5 10])
%             yticklabels({'10','5','1'});
%             ylabel('Ensembles');
%             xlabel('Time (s)');
%             set(gca,'fontsize',16)
            
            
        end
    end
end
            
            
pbaspect([24 1 1])
