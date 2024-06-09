% current path:
currentScriptPath = fileparts(mfilename('fullpath'));
% generate string for all subdirectories:
pathString = genpath(currentScriptPath);
% add all subdirectories to path
addpath(pathString);

% Load F4 data:
load('F4.mat')

% set default options for clustering using clustheat_defaults() fn.
defaultOptions = clustheat_defaults();
options = defaultOptions; % and use default options
%%
for a = 1:length(unique(neuron_to_animal))
    od_abytrials(a,:) = pooled_olist(trial_to_animal==a);
end

for k = 1:size(pooled_F3,1)
    for o = 1:6
        selectF3(k,:,o) = mean(pooled_F3(k,:,...
            od_abytrials(neuron_to_animal(k),:)==options.odord(o)),3);
    end
end

Fsiz = size(selectF3);
selectF2 = reshape(selectF3,[Fsiz(1),Fsiz(2)*Fsiz(3)]);
normF2 = normalize(selectF2,2);
normF3 = reshape(normF2,[Fsiz(1),Fsiz(2),Fsiz(3)]);
%%
odbins = 71:75;
alist = unique(animid);

for a = 1:length(alist)
    tempFsnap = permute(mean(normF3(animid==alist(a),odbins,[1,5,3,6,2,4]),...
        2),[1,3,2]);
    temp2d = squareform(pdist(tempFsnap','euclidean')./...
        max(pdist(tempFsnap','euclidean'),[],'all'));
    temp2d(eye(size(temp2d))==1) = nan;
    temp2d(tril(ones(size(temp2d,1)))==1) = nan;
    pair_odist(:,:,a) = temp2d;
    clear tempFsnap;
end
%%
intra_or_inter = categorical(repmat("inter",size(pair_odist))); %initialize zero matrix
for i = 1:size(intra_or_inter,3)
    temp_mat = intra_or_inter(:,:,i);
    temp_mat(tril(ones(size(intra_or_inter,1))==1)) = "exclude";
    intra_or_inter(:,:,i) = categorical(temp_mat);
end
for i = [1,3,5]
    intra_or_inter(i,i+1,:)="intra"; %set intra-valence elements to 'intra'
end
%%
%to get the average intra or inter odist for each animal:
for a = 1:size(pair_odist,3)
    tempdist = pair_odist(:,:,a);
    intra_by_anim(a) = mean(tempdist(intra_or_inter(:,:,a)=="intra"),'all');
    inter_by_anim(a) = mean(tempdist(intra_or_inter(:,:,a)=="inter"),'all');
    clear tempdist
end

% calculate mean and sem for each
inter_m = mean(inter_by_anim); 
intra_m = mean(intra_by_anim);
inter_s = std(inter_by_anim,[],2)/sqrt(length(alist));
intra_s = std(intra_by_anim,[],2)/sqrt(length(alist));
%%
%to run stats on this we just extract the upper triangle minus the diag:
Y_dist = pair_odist(:);
Y_dist = Y_dist(~(isnan(Y_dist)));
X_label = intra_or_inter(~(intra_or_inter=="exclude"));
A_label = categorical(reshape(repmat(alist',15,1),[],1));
%%
% the average pairdist heatmap. for this we want to leave the
% lower triangle:
for a = 1:length(alist)
    tempFsnap = permute(mean(normF3(animid==alist(a),odbins,[1,5,3,6,2,4]),...
        2),[1,3,2]);
    temp2d = squareform(pdist(tempFsnap','euclidean')./...
        max(pdist(tempFsnap','euclidean'),[],'all'));
    temp2d(eye(size(temp2d))==1) = nan;
    pair_odist(:,:,a) = temp2d;
    clear tempFsnap;
end

mu_pdist_2d = (mean(pair_odist,3));
mu_pdist_2d(mu_pdist_2d==0) = nan;

panelS3E = figure('Position',[0,0,520,500]);
h2 = heatmap(mu_pdist_2d); 
colormap parula
h2.XData = {"4MT","TMT","IAA","Heptanol","Peanut Oil","2PE"};
h2.YData = {"4MT","TMT","IAA","Heptanol","Peanut Oil","2PE"};

h2.CellLabelFormat = '%.2f';
set(gca,'FontName','Arial','FontSize',15);