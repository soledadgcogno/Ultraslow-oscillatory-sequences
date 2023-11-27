function [table_motion,table_still,motion_2]=get_table_motion_updated(speed_d,threshold_still,gap_for_motion,gap_for_still)

% threshold_still=3;
% gap_for_motion=1; %in seconds
speed_d(speed_d<threshold_still)=0;

T=length(speed_d)+1;
motion=ones(1,T);
motion(speed_d==0)=0;
d_motion=diff(motion);
motion_start=zeros(1,T);
motion_fin=zeros(1,T);
motion_start(d_motion==1)=1;
motion_fin(d_motion==-1)=1;


% This is in case we want to define a minimum gap of immobility that is
% considers motion, and the other way around
immobil_bet_motion=8*gap_for_motion;
motion_2=1*motion;
aux=(find(motion==1));
aux_d=diff(aux); %inter-motion-intervals
for a=1:length(aux_d)
    if aux_d(a)>1 && aux_d(a)<immobil_bet_motion 
        motion_2(aux(a):aux(a+1))=1;
    end
end

clear aux aux_d
immobil_bet_still=8*gap_for_still;
motion_3=motion_2;
aux=(find(motion_2==0));
aux_d=diff(aux);%inter-still-intervals
for a=1:length(aux_d)
    if aux_d(a)>1 && aux_d(a)<immobil_bet_still
        motion_3(aux(a):aux(a+1))=0;
    end
end

% Now we create the table
clear motion_2
motion_2=motion_3;
d_motion_2=diff(motion_2);
motion_start_2=zeros(1,T);
motion_fin_2=zeros(1,T);
motion_start_2(d_motion_2==1)=1;
motion_fin_2(d_motion_2==-1)=1;
motion_start_2_frames=find(motion_start_2>0)+1;
motion_fin_2_frames=find(motion_fin_2>0);

if motion_2(1)==1 %In case the animal is moving from the first frame
    motion_start_2_frames(2:length(motion_start_2_frames)+1)=motion_start_2_frames;
    motion_start_2_frames(1)=1;
end

if motion_start_2_frames(1)==1 && motion_2(end)==1 %Beggining = Run & End = Run
    
    for a=1:length(motion_start_2_frames)-1
        table_motion(a,1)=motion_start_2_frames(a);
        table_motion(a,2)=motion_fin_2_frames(a);
        
        table_still(a,1)=motion_fin_2_frames(a)+1;
        table_still(a,2)=motion_start_2_frames(a+1)-1;
        
        if a==length(motion_start_2_frames)-1
            table_motion(a+1,1)=motion_start_2_frames(a+1);
            table_motion(a+1,2)=T;
        end
    end
    
    
elseif motion_2(1)==1 && motion_2(end)==0 %Beggining = Run & End = Still
    
    for a=1:length(motion_start_2_frames)
        table_motion(a,1)=motion_start_2_frames(a);
        table_motion(a,2)=motion_fin_2_frames(a);
        
        if a==length(motion_start_2_frames)
            table_still(a,1)=motion_fin_2_frames(a)+1;
            table_still(a,2)=T;
        else
            table_still(a,1)=motion_fin_2_frames(a)+1;
            table_still(a,2)=motion_start_2_frames(a+1)-1;
        end
    end
    
elseif motion_start_2_frames(1)>1 && motion_2(end)==1 %Beggining = Still & End = Run
    
    for a=1:length(motion_start_2_frames)-1
        table_motion(a,1)=motion_start_2_frames(a);
        table_motion(a,2)=motion_fin_2_frames(a);
        
        table_still(a+1,1)=motion_fin_2_frames(a)+1;
        table_still(a+1,2)=motion_start_2_frames(a+1)-1;
        
        if a==1
            table_still(a,1)=1;
            table_still(a,2)=motion_start_2_frames(a)-1;
        end
                
        if a==length(motion_start_2_frames)-1
            table_motion(a+1,1)=motion_start_2_frames(a+1);
            table_motion(a+1,2)=T;
        end
        
    end
    
elseif motion_start_2_frames(1)>1 && motion_2(end)==0 %Beggining = Still & End = Still
    
    for a=1:length(motion_start_2_frames)
        table_motion(a,1)=motion_start_2_frames(a);
        table_motion(a,2)=motion_fin_2_frames(a);
        
        if a==length(motion_start_2_frames)
            table_still(a,1)=motion_fin_2_frames(a)+1;
            table_still(a,2)=T;
        else
            table_still(a,1)=motion_fin_2_frames(a)+1;
            table_still(a,2)=motion_start_2_frames(a+1)-1;
        end
        
        if a==1
            table_still(a,1)=1;
            table_still(a,2)=motion_start_2_frames(a)-1;
        end
        
    end
    
end


end
