function [table_u_r,N,T]=identify_waves_latestversion_6_f(mouse_num,day,num_clus_discr,dt,make_fig,spikes)
sf=7.73; 
min_amp=5; %used to be 2
[N,T]=size(spikes);
[~,sorting,~] = get_sorting(spikes);


for i=1:N
%     FRp(i,:)=full(fire_rate(spikes(i,:),1*dt,'g')); %smooth in the dt chosen for each session
      FRp(i,:)=smoothdata(spikes(i,:),'gaussian',1*dt); %smooth in the dt chosen for each session
end

[~,sorting_FRp,~] = get_sorting(FRp); %gets sorting based on smoothed spike matrix


% [~,scoret,~] = pca(spikes');
% phase_f=(atan2(smooth(scoret(:,2),1*dt),smooth(scoret(:,1),1*dt)));
% phase_r=(atan2(scoret(:,2),scoret(:,1)));

[~,scoret,~] = pca(FRp');
phase_r=(atan2(scoret(:,2),scoret(:,1)));
phase_f=(atan2(smooth(scoret(:,2),floor(1*dt)),smooth(scoret(:,1),floor(1*dt))));


%Check whether waves go upwards or inwards
aux=(diff(phase_f));
pos=length(find(aux>5));
neg=length(find(aux<-5));

if pos<neg
    sign_phase=1;
else
    phase_f=-phase_f;
    sign_phase=1;
    sorting=flip(sorting);
    sorting_FRp=flip(sorting_FRp);
end

%Discretization of the phase and phase difference
alfa= discretize(phase_f,-pi:2*pi/num_clus_discr:pi);
alfac=alfa;

thr=0.3*num_clus_discr;
aux=(diff(alfac));
aux_ini=find(abs(aux)>thr); %Intervals for which the phase jumps

phase2_s_mod=phase_f;

%Some hyperparameters
tol_up=5;
tol_up_in=3;
tol_down=0.7;

count=0;
for i=1:length(aux_ini)-1 %Loop on all intervals
    start=aux_ini(i)+1;
    finish=aux_ini(i+1);
    
    %Removes sustained coactivity at the beggining or at the end of the
    %interval
    range=alfac(start: finish);
    [min_en,mi] = min(range);
    [max_en,ma_m] = max(range(mi:end));
    aux_diff=diff(range(mi:end));
    if isempty(find(aux_diff(ma_m:end),1))
        ma=mi+ma_m-1;
    else
        ma=mi+ma_m-1 + find(aux_diff(ma_m:end),1) - 1;
    end
    
    if find(aux_diff,1)~=1
        mi=find(aux_diff,1);
    end
    
    %If an interval begins with activity at the top that then switches to
    %the bottom to let a wave begin, we delete the first part and only keep
    %the wave
    if(length(range)>3)
        [~,loc_max_range]=findpeaks(range);  %local maxima in range
        [~,loc_min_range]=findpeaks(-range);  %local minima in range
        
        if (isempty(loc_max_range)==0 && isempty(loc_min_range)==0 )
            %         if (isempty(loc_min_range)==0 )
            if(loc_min_range(1)<loc_max_range(1))
                range=range(loc_min_range(1):end);
                start= start + loc_min_range(1);
                [min_en,mi] = min(range);
                [max_en,ma_m] = max(range(mi:end));
                aux_diff=diff(range(mi:end));
                if isempty(find(aux_diff(ma_m:end),1))
                    ma=mi+ma_m-1;
                else
                    ma=mi+ma_m-1 + find(aux_diff(ma_m:end),1) - 1;
                end
            end
        end
        
        if (isempty(loc_max_range)==1 && isempty(loc_min_range)==0 )
            range=range(loc_min_range(1):end);
            start= start + loc_min_range(1);
            [min_en,mi] = min(range);
            [max_en,ma_m] = max(range(mi:end));
            aux_diff=diff(range(mi:end));
            if isempty(find(aux_diff(ma_m:end),1))
                ma=mi+ma_m-1;
            else
                ma=mi+ma_m-1 + find(aux_diff(ma_m:end),1) - 1;
            end
        end
    end
    
    clear loc_max_range loc_min_range range_u
    
    if ((max_en-min_en>min_amp) ) %Condition on minimum amplitude for the waves
                %eliminates 1 ensemble fluctuations in range
%         range_sm=range;
%         trans=find(aux_diff);
%         for m=1:length(trans)-3
%             if( range(trans(m)) == range(trans(m+1)+1) )
%                 disp(m)
%                 range_sm(trans(m):trans(m+1)+1)=range(trans(m));
%                 m=m+1;
%             end
%         end
        
        
        range_sub=range(mi:ma);        
        range_u(1)=range_sub(1);
        range_u(2:length(find(diff(range_sub)~=0))+1) = range_sub(find(diff(range_sub)~=0) + 1);
        
        %eliminates 1 ensemble fluctuations in range_u
        rev=find(diff(range_u)==-1);
        aux=diff(range_u);
        for l=1:length(rev)
            asc=find(aux(1:rev(l))==1,1,'last');
            range_u(asc:rev(l))=range_u(asc);
        end
        ini=find(diff(range_u),1);
        range_u=range_u(ini:end);
        
        template=min_en:max_en;
        l_template=length(template);
        l_rangeu=length(range_u);
        
        if (l_rangeu<=tol_up*l_template && l_rangeu>tol_down*l_template && length(find(diff(range_u)))/length(range_u)>0.5)
                                 
            %Splits into smaller waves
            
            [maxs,loc_max]=findpeaks(range_u);  %local maxima in range_u
            [mins,loc_min]=findpeaks(-range_u); %local minima in range_u
            
            [maxs_range,loc_max_range]=findpeaks(range);  %local maxima in range
            [mins_range,loc_min_range]=findpeaks(-range);  %local minima in range
            
                        
            %This matches the number of local maxima and minima in range
            %and range_u
            if (length(loc_max_range) > length(loc_max) && ~isempty(loc_max))
                for m=1:length(loc_min_range)
                    if ( maxs_range(m)+mins_range(m) > 1)
                        flag_m(m)=1;
                    else
                        flag_m(m)=0;
                    end
                end
                loc_max_range(find(flag_m==0))=[];
                loc_min_range(find(flag_m==0))=[];
            end
            
            
            
            %This matches the number of local maxima and minima in range
            %and range_u
%             if (length(loc_max_range) > length(loc_max) && ~isempty(loc_max))
%                 aux=find(diff(range));
%                 for m=1:length(loc_max_range)
%                     [~,b]= min(abs(aux-loc_max_range(m)));
%                                         
%                     if ((b+3)<=length(aux) && (b-2)>0)
%                         af1 = range(aux(b+2));
%                         br1 = range(aux(b));
%                         af2 = range(aux(b+3));
%                         br2 = range(aux(b-1));
%                         
%                         if (af1>af2 && br1>br2)
%                             flag_m(m)=1;
%                         else
%                             flag_m(m)=0;
%                         end                        
%                     else
%                         flag_m(m)=0;
%                     end                    
%                     clear af1 af2 br1 br2
%                 end         
%                 
%                 count_flag=0;
%                 for k=find(flag_m==0)   
%                     min_from_max=loc_min_range-loc_max_range(k);
%                     min_from_max(min_from_max<=0)=nan;
%                     [~,min_from_max_in]=min(min_from_max);
%                     count_flag=count_flag+1;
%                     flag_min(count_flag)=min_from_max_in;
%                     clear min_from_max  min_from_max_in
%                 end
%                 
%                 loc_max_range(find(flag_m==0))=[];
%                 loc_min_range(flag_min)=[];
%                     
%                 clear flag_m flag_min
%             end
            
            
            
%             if (length(loc_min_range) > length(loc_min) && ~isempty(loc_min))
%                 aux=find(diff(range));
%                 for m=1:length(loc_min_range)
%                     [~,b]= min(abs(aux-loc_min_range(m)));
%                                         
%                     if ((b+3)<=length(aux) && (b-2)>0)
%                         af1 = range(aux(b+2));
%                         br1 = range(aux(b));
%                         af2 = range(aux(b+3));
%                         br2 = range(aux(b-1));
%                         
%                         if (af1<af2 && br1<br2)
%                             flag_m(m)=1;
%                         else
%                             flag_m(m)=0;
%                         end                        
%                     else
%                         flag_m(m)=0;
%                     end                    
%                     clear af1 af2 br1 br2
%                 end                
%                 
%                  count_flag=0;
%                 for k=find(flag_m==0)   
%                     max_from_min=loc_max_range-loc_min_range(k);
%                     max_from_min(max_from_min>=0)=nan;
%                     [~,min_from_max_in]=max(max_from_min);
%                     count_flag=count_flag+1;
%                     flag_max(count_flag)=min_from_max_in;
%                     clear min_from_max  min_from_max_in
%                 end
%                 
%                 loc_max_range(flag_max)=[];
%                 loc_min_range(find(flag_m==0))=[];
%                 clear flag_m
%             end
            
            %%Asigns the interval within intervals that are ocnsistent with
            %%taves to table_u
            if isempty(loc_max)==1
                count=count+1;
                table_u(count,1)=start + mi-1;
                table_u(count,2)=start + ma-1 +1 +1;
            else
                for l=1:length(loc_max)
                    if -mins(l)>=2
                        if l==1
                            len=length(1:loc_min(1));
                            if ((maxs(l)-range_u(1))>=3 && len<tol_up_in*(maxs(l)-range_u(1)))
                                count=count+1;
                                table_u(count,1)=start + mi-1;
                                end_point=find(diff(range(loc_max_range(l):loc_min_range(l)))~=0,1);
                                table_u(count,2)=start + loc_max_range(l)+end_point-1;
                            end
                        else
                            len=length(loc_min(l-1):loc_min(l));
                            if ((maxs(l)+mins(l-1))>=3 && len<tol_up_in*(maxs(l)-range_u(1)))
                                count=count+1;                             

                                table_u(count,1)=start + loc_min_range(l-1) + 1;
                                end_point=find(diff(range(loc_max_range(l):loc_min_range(l)))~=0,1);
                                table_u(count,2)=start + loc_max_range(l) + end_point -1 -1;
                            end
                        end
                    end
                end
                
                if((range_u(end)+mins(end))>=3)
                    count=count+1;
%                     if length(loc_min_range)>length(loc_min)
%                         table_u(count,1)=start + loc_min_range(end-1) + 1;
%                     else
                        table_u(count,1)=start + loc_min_range(end) + 1;
%                     end
                    table_u(count,2)=start + ma-1 +1 +1;
                end
            end
        end
    end
    
    clear range slope b0 range_sub template range_u range_u template loc_max loc_min maxs mins loc_max_range loc_min_range maxs_range mins_range flag_m flag
end


%Simplified version of the previous code for border_beggining
start=1;
finish=aux_ini(1)-1;

range=alfac(start: finish);
[min_en,mi] = min(range);
[max_en,ma_m] = max(range(mi:end));
aux_diff=diff(range(mi:end));
if isempty(find(aux_diff(ma_m:end),1))
    ma=mi+ma_m-1;
else
    ma=mi+ma_m-1 + find(aux_diff(ma_m:end),1) - 1;
end

if ((max_en-min_en>min_amp )) %Should be >= because a bit upward when we have combination os cycles we allow for 5 ensembles cycles
    range_sub=range(mi:ma);
    
    range_u(1)=range_sub(1);
    range_u(2:length(find(diff(range_sub)~=0))+1) = range_sub(find(diff(range_sub)~=0) + 1);
    
    %eliminates 1 ensemble fluctuations
    rev=find(diff(range_u)==-1);
    aux=diff(range_u);
    for l=1:length(rev)
        asc=find(aux(1:rev(l))==1,1,'last');
        range_u(asc:rev(l))=range_u(asc);
    end
    ini=find(diff(range_u),1);
    range_u=range_u(ini:end);
    
    
    template=min_en:max_en;
    l_template=length(template);
    l_rangeu=length(range_u);
    
    if (l_rangeu<=tol_up*l_template && l_rangeu>tol_down*l_template && length(find(diff(range_u)))/length(range_u)>0.5)
       
        
        %Splits into smaller waves
        
        [maxs,loc_max]=findpeaks(range_u);  %local maxima in range_u
        [mins,loc_min]=findpeaks(-range_u); %local minima in range_u
        
        [maxs_range,loc_max_range]=findpeaks(range);  %local maxima in range
        [mins_range,loc_min_range]=findpeaks(-range);  %local minima in range
        
        
        if isempty(loc_max)==1
            count=count+1;
            table_u(count,1)=start + mi-1;
            table_u(count,2)=start + ma-1 +1 +1;
            table_u=circshift(table_u,1,1);

        else
            for l=1:length(loc_max)
                if -mins(l)>=2
                    if l==1
                        len=length(1:loc_min(1));
                        if ((maxs(l)-range_u(1))>=3 && len<tol_up_in*(maxs(l)-range_u(1)))
                            count=count+1;
                            table_u(count,1)=start + mi-1;
                            table_u(count,2)=start + loc_max_range(l);
                        end
                    else
                        len=length(loc_min(l-1):loc_min(l));
                        if ((maxs(l)+mins(l-1))>=3 && len<tol_up_in*(maxs(l)-range_u(1)))
                            count=count+1;
                            table_u(count,1)=start + loc_min_range(l-1) + 1;
                            table_u(count,2)=start + loc_max_range(l)-1;
                        end
                    end
                end
            end
            
            if((range_u(end)+mins(end))>=3)
                
                count=count+1;
                table_u(count,1)=start + loc_min_range(end) + 1;
                table_u(count,2)=start + ma-1 +1 +1;                
                table_u=circ_shift(table_u,1,1);
            end
        end
    end
end
clear range slope b0 range_u  range_sub template range_u range_u template loc_max loc_min maxs mins loc_max_range loc_min_range maxs_range mins_range


%Simplified version of the previous code for border_end
start=aux_ini(end)+1;
finish=length(alfac);

range=alfac(start: finish);
[min_en,mi] = min(range);
[max_en,ma_m] = max(range(mi:end));
aux_diff=diff(range(mi:end));
if isempty(find(aux_diff(ma_m:end),1))
    ma=mi+ma_m-1;
else
    ma=mi+ma_m-1 + find(aux_diff(ma_m:end),1) - 1;
end

if ((max_en-min_en>min_amp) )
    range_sub=range(mi:ma);
    
    range_u(1)=range_sub(1);
    range_u(2:length(find(diff(range_sub)~=0))+1) = range_sub(find(diff(range_sub)~=0) + 1);
    
    %eliminates 1 ensemble fluctuations
    rev=find(diff(range_u)==-1);
    aux=diff(range_u);
    for l=1:length(rev)
        asc=find(aux(1:rev(l))==1,1,'last');
        range_u(asc:rev(l))=range_u(asc);
    end
    ini=find(diff(range_u),1);
    range_u=range_u(ini:end);
    
    
    template=min_en:max_en;
    l_template=length(template);
    l_rangeu=length(range_u);
    
    if (l_rangeu<=tol_up*l_template && l_rangeu>tol_down*l_template && length(find(diff(range_u)))/length(range_u)>0.5)
        %             hold on
        %             plot(range,'r')
        
        %             disp(i)
        %             count=count+1;
        %             table_u(count,1)=start + mi-1;
        %             table_u(count,2)=start + ma-1 +1 +1;
        
        
        %Splits into smaller waves
        
        [maxs,loc_max]=findpeaks(range_u);  %local maxima in range_u
        [mins,loc_min]=findpeaks(-range_u); %local minima in range_u
        
        [maxs_range,loc_max_range]=findpeaks(range);  %local maxima in range
        [mins_range,loc_min_range]=findpeaks(-range);  %local minima in range
        
        
        if isempty(loc_max)==1
            count=count+1;
            table_u(count,1)=start + mi-1;
            table_u(count,2)=start + ma-1 +1 +1;
        else
            for l=1:length(loc_max)
                if -mins(l)>=2
                    if l==1
                        len=length(1:loc_min(1));
                        if ((maxs(l)-range_u(1))>=3 && len<tol_up_in*(maxs(l)-range_u(1)))
                            count=count+1;
                            table_u(count,1)=start + mi-1;
                            table_u(count,2)=start + loc_max_range(l);
                        end
                    else
                        len=length(loc_min(l-1):loc_min(l));
                        if ((maxs(l)+mins(l-1))>=3 && len<tol_up_in*(maxs(l)-range_u(1)))
                            count=count+1;
                            table_u(count,1)=start + loc_min_range(l-1) + 1;
                            table_u(count,2)=start + loc_max_range(l)-1;
                        end
                    end
                end
            end
            
            if((range_u(end)+mins(end))>=3)
                count=count+1;
                table_u(count,1)=start + loc_min_range(end) + 1;
                table_u(count,2)=start + ma-1 +1 +1;
            end
        end
        if(table_u(end,2)>length(phase_f))
            table_u(end,2)=length(phase_f);
        end
    end
end
clear range slope b0 range_u  range_sub template range_u range_u template loc_max loc_min maxs mins loc_max_range loc_min_range maxs_range mins_range

%Removes putative waves that are too short when compared with the time
%scales of the waves
iwi=(table_u(:,2)-table_u(:,1))/sf;
scale=2*dt/sf;
too_short=find(iwi<scale);
table_u(too_short,:)=[];


for wa=1:size(table_u,1)
    cont_min=identification_wave_expand_min(phase_f,table_u,wa);
    cont_max=identification_wave_expand_max(phase_f,table_u,wa);
    table_u_corrected(wa,1)=table_u(wa,1)-cont_min;
    table_u_corrected(wa,2)=table_u(wa,2)+cont_max;
end

table_u_old=table_u;
table_u=table_u_corrected;

%Figure
spikes_d_art=zeros(size(spikes));

table_u_r=table_u*1;
for i=1:size(table_u,1)
    spikes_d_art(:,table_u_r(i,1):table_u_r(i,2))=spikes(:,table_u_r(i,1):table_u_r(i,2));
end

sorting_FRp=flip(sorting_FRp);
if make_fig==1
    phase=phase_f;
    figure
    subplot(2,1,1)
    %     spy(flip(spikes(sorting_FRp,:)),'k')
    %     pbaspect([24 1 4])
    %     %         title(['SNR = ',num2str(snrtt_full)])
    %     subplot(3,1,2)
    %     spy(flip(spikes(sorting_FRp,:)),'k')
    %     pbaspect([24 1 8])
    %     hold on
    %     spy(flip(spikes_d_art(sorting_FRp,:)),'r')
    %     pbaspect([24 1 8])
    %     title([mouse_num,' ', num2str(day), ' dt = ' ,num2str(dt)]);
    
    %       fig=figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.3]);
    hold on
    for i=1:size(spikes,1)
        %                     scatter((1:size(spikes_d,2))./8,i*spikes_d(sorting_descend(i),:),5,'k','filled')
        scatter((1:size(spikes,2))./sf,i*spikes(sorting_FRp(i),:),5,'k','filled')
        alpha 0.2
    end
    axis([-inf inf 1 inf]);
    %title([mouse,' Day',num2str(day)]);
    ylabel('Neurons #');
    set(gca,'fontsize',18);
    xticks([]);
    yticks([100 400])
    cc=plasma(size(table_u,1));
    subplot(2,1,2)
    plot((1:size(spikes,2))./sf,phase,'k','linewidth',1.5);
    hold on
    for i=1:size(table_u,1)
        plot((table_u(i,1):table_u(i,2))./sf,phase(table_u(i,1):table_u(i,2)),'linewidth',2,'color',cc(i,:));
    end
    pbaspect([24 1 5])
    %     ylabel('Phase (rad)');
    %     xlabel('Time (s)');
    axis([-inf inf -inf inf])
    yticks([-3.14,0,3.14])
    yticklabels({'-\pi','0','\pi'})
    set(gca,'fontsize',18)
    xlabel('Time (s)');
    ylabel('Phase (rad)')
    box off
end


end


