function responsive_matrix = isresponsive(F4,odor,base_start,base_end,stim_start,stim_end)
%responsive_matrix = isresponsive(F4,odor,base_start,base_end,stim_start,stim_end)
%
%F4 is a 4d matrix of shape [k, frames, trials, odors]
%[base_start, base_end] defines timebin of baseline F
%[stim_start, stim_end] defines timebin for stim F
%odor specifies the 4th dimension of F4 to look at
%
%outputs:
%a matrix of length k and width length(stim_start:stim_end) that is all the
%pvalues for ranksum tests
%%
%1. pool baseline DF/F to get reference distribution
%first, we get the size of F4 and give them useful names:
k = size(F4,1); frames = size(F4,2); trials = size(F4,3); classes = size(F4,4);
base_length = length(base_start:base_end);
%then we pool all the baseline frames
Fbase = reshape(F4(:,base_start:base_end,:,:),[k,base_length*trials*classes]);
%%
%2. reduce the dataset to the stim frames and specified odor
Fstim = F4(:,stim_start:stim_end,:,odor);
%%
%3. for each neuron, for each odor, for each frame after odor onset, test against baseline to get a p-value
%   shape of third matrix is [k,10,6]
responsive_matrix = zeros(k,size(Fstim,2));
for n = 1:k
    x = Fbase(n,:);
    for t = 1:size(Fstim,2)
        y = Fstim(n,t,:);
        responsive_matrix(n,t) = ranksum(x(:),y(:));
    end
end

end