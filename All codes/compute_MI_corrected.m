function MI = compute_MI_corrected(table,disc_var1,disc_var2)

table(:,1)=discretize(table(:,1),0:disc_var1+1);
table(:,2)=discretize(table(:,2),disc_var2);

%Joint probability
for l=1:disc_var1+1
    for j=1:disc_var2
        aux=find(table(:,1)==l);
        aux2=find(table(aux,2)==j);
        joint_p(l,j)=length(aux2);%./size(table,1);
        
        clear aux aux2
    end
end

if (sum(sum(joint_p)))~=size(table,1)
    disp('Error: Joint prob');
end

 joint_p= joint_p./size(table,1);
 
%Marginal probabilities
prob_var1=sum(joint_p,2);
prob_var2=sum(joint_p,1);

H_joint=-sum(nonzeros(joint_p(:)).*log2(nonzeros(joint_p(:))));
H_var2=-sum(nonzeros(prob_var2).*log2(nonzeros(prob_var2)));
H_var1=-sum(nonzeros(prob_var1).*log2(nonzeros(prob_var1)));

MI = H_var1 + H_var2 - H_joint;
end