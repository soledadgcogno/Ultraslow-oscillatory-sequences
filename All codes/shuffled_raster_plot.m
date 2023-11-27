clear all
close all

dbeh_path='C:\Users\xscogno\MATLAB\Flavio2\Waves\Tracking data\';
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';


mice_number=6;
mice=['L08M1';'L08M2';'L08M3';'L08M4';'L09M1';'L09M4';'L05M2';'L05M3';'L05M5';'92227';'92229';'60961'];
% mice=['L8M1';'L8M2';'L8M3';'L8M4';'L9M1';'L9M4'];


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
        for s=1%:dates.sesnum(day)
            disp(s)
            munit=dates.ses{day}(s);
            file_name=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];         

            if (exist(file_name) == 2)
                disp(day)
                load(file_name,'-mat');
                spikes_d=full(spikes_d_s);
                [N,T]=size(spikes_d);

                % Angle PCA
                [coeff,score]=pca(spikes_d');
                [sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes_d);
               
                sorting_pca=sorting_descend;
                figure
                spy(spikes_d(sorting_pca,:),'k')
                pbaspect([23,2,1])

                offfset=2*60*8;
                for i=1:N
                    spk_sh(i,:)=circshift(spikes_d(i,:),offfset+floor(rand*(size(spikes_d,2)-offfset)));
                end
                % Angle PCA - shuffle
                [coeff,score]=pca(spk_sh');
                [sorting_ascend,sorting_descend,sorting_0]=get_sorting(spk_sh);
%                 [~,sorting_w,~]=get_sorting_smoothed(spikes_d,117);
                sorting_pca=sorting_descend;
                figure
                spy(spk_sh(sorting_pca,:),'k')
                pbaspect([23,2,1])


            end
        end
    end
end

%%
sf=7.73;
% PCA
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.4]);
sorting_pca=flip(sorting_pca);
hold on
for i=1:size(spikes_d,1)
    scatter((1:size(spk_sh,2))./8,i*spk_sh((sorting_descend(i)),:),5,'k','filled')
    alpha 0.2
end
axis([-inf inf 1 inf])
ylabel('Neurons #');
xlabel('Time (s)');
set(gca,'fontsize',18)
yticks([100 200])
ini=sorting_ascend(1);
set(gcf, 'Renderer', 'opengl');
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Revision\New figures\raster_plot\Rasterplot L9M4Day17 - Circular_shuffling - PCA sorting.fig');
saveas(gcf,'C:\Users\xscogno\Dropbox\Waves\Revision\New figures\raster_plot\Rasterplot L9M4Day17 - Circular_shuffling - PCA sorting.svg');


close all