function [coac_vec_1,s_vec] = rcoact_dist(spike_train_m)
%Returns a vector with the coactivity distribution of a matrix of spike
%trains. The matrix should have n rows, corresponding to the neurons, and p
%columns, corresponding to the time frames. The output is an n dimensional
%vector containing the number of counts for each possible amount of
%coactive neurons, from 1 to n.

s_vec = sum(spike_train_m);
coac_vec = zeros(size(spike_train_m,1)+1,1);
for i=1:size(coac_vec)
    coac_vec_1(i) = sum(s_vec(:)==i-1)/size(spike_train_m,2);
end
end