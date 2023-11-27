function [C, idx] = bin_cells_position(r_i,bin_size,max_length)

%r_i            position in (x,y) of each cell
%max_length     is the maximum distance in one direction (um)
%bin_size       is bin size in one diretion (um)

N=size(r_i,1);
n_bins=max_length/bin_size;

posx=discretize(r_i(:,1),0:bin_size:max_length);
posy=discretize(r_i(:,2),0:bin_size:max_length);

C1=[];
C2=[];
for l=1:n_bins
    C1 = [C1;(l*ones(n_bins,1))];
    C2 = [C2;(1:n_bins)'];
end
C=[C1,C2];

for n=1:N
    a=find(posx(n)==C(:,1));
    b=find(posy(n)==C(:,2));
    idx(n)=intersect(a,b);
    clear a b
end