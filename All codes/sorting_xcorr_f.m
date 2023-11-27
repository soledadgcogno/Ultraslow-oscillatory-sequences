function [sorting] = sorting_xcorr_f(spikes,maxlag,downsampling_factor,dt)

[N,~]=size(spikes);

FRp = spikes_downsample(spikes,N,downsampling_factor);


for i=1:N
    FR(i,:)=full(fire_rate(FRp(i,:),dt,'g')); 
end


% FR=FRp;
for i=1:N
    for j=i+1:N
        [val,time]=(xcorr(FR(i,:),FR(j,:),maxlag));
        [v,in]=max((val));
        if time(in)>=0
            Corr_mat(i,j)=v;
            Corr_mat(j,i)=-v;
        else
            Corr_mat(i,j)=-v;
            Corr_mat(j,i)=v;
        end
    end
end

%Based on one neuron
[r,c]=max(max(Corr_mat));
[~,sorting]=sort(Corr_mat(c,:),'descend');
% % % 
% % % figure
% % % spy(spikes(sorting,:),'k',5)
% % % pbaspect([23,2,1])
% % % axis off




end




%% Old version


% % 
% % [N,~]=size(spikes);
% % 
% % FRp = spikes_downsample(spikes,N,downsampling_factor);
% % 
% % for i=1:N
% %     FR(i,:)=full(fire_rate(FRp(i,:),dt,'g')); %Smoothing in about 10 seconds
% % end
% % 
% % for i=1:N
% %     for j=i+1:N
% %         rho=corr(FR(i,:)',FR(j,:)');
% %         Corr_mat(i,j)=rho;
% %         Corr_mat(j,i)=-rho;
% %     end
% % end
% % 
% % %Based on one neuron
% % [r,c]=max(max(Corr_mat));
% % [~,sorting]=sort(Corr_mat(c,:),'descend');
% % % % %
% % % % % figure
% % % % % spy(spikes(sorting,:),'k',5)
% % % % % pbaspect([23,2,1])
% % % % % axis off
% % 

