% current path:
currentScriptPath = fileparts(mfilename('fullpath'));
% generate string for all subdirectories:
pathString = genpath(currentScriptPath);
% add all subdirectories to path
addpath(pathString);

% Load F4 data:
load('F4.mat')
%%
% grab the trial order from each animal
for a = 1:length(unique(neuron_to_animal))
    od_abytrials(a,:) = pooled_olist(trial_to_animal==a);
end

% standardize the ordering of trials
for k = 1:size(pooled_F3,1)
    for o = 1:6
        F3(k,:,o) = mean(pooled_F3(k,:,...
            od_abytrials(neuron_to_animal(k),:)==options.odord(o)),3);
    end
end
%%
% concatenate trial-averaged responses
Fsiz = size(F3);
F2 = reshape(F3,[Fsiz(1),Fsiz(2)*Fsiz(3)]); 
% normalize across time
normF2 = normalize(F2,2);
% restack the normalized matrix so that the 3rd dimension is different
% odors
normF3 = reshape(normF2,[Fsiz(1),Fsiz(2),Fsiz(3)]);
%%
% check which animals were anterior or posterior
check_ant = @(x) strcmp(aloc{x}, 'a');
check_post = @(x) strcmp(aloc{x}, 'p');

% apply cellfun from above to animid:
isant = cellfun(check_ant, num2cell(animid));
ispost = cellfun(check_post, num2cell(animid));
%%
%do the analysis on the last second of the odor:
odbins = 71:75;
%now calculate the valence score defined as abs(Fnut-F2pe) - abs(Fpos-Fneg)
Fsnap = mean(normF3(:,odbins,:),2);

inter_idx = [1,2;1,4;5,2;5,4];
intra_idx = [1,5;2,4;3,6];

for i = 1:size(inter_idx,1)
    inter_diff(:,i) = (Fsnap(:,:,inter_idx(i,2)) - Fsnap(:,:,inter_idx(i,1)));
end

for j = 1:size(intra_idx,1)
    intra_diff(:,j) = abs(Fsnap(:,:,intra_idx(j,1)) - Fsnap(:,:,intra_idx(j,2)));
end

valsco = mean(inter_diff,2)./mean(intra_diff,2);

idx_perm = perms(1:6);
for s = 1:size(idx_perm,1)
    Shfsnap = Fsnap(:,:,idx_perm(s,:));
    for i = 1:size(inter_idx,1)
        sinter_diff(:,i,s) = (Shfsnap(:,:,inter_idx(i,2)) - Shfsnap(:,:,inter_idx(i,1)));
    end

    for j = 1:size(intra_idx,1)
        sintra_diff(:,j,s) = abs(Shfsnap(:,:,intra_idx(j,1)) - Shfsnap(:,:,intra_idx(j,2)));
    end

    valshf(:,s) = mean(sinter_diff(:,:,s),2)./mean(sintra_diff(:,:,s),2);
end

%%
panelI = figure('Position',[0,0,450,500]);
% if visualizing ant vs post separately:
bandwidth= 0.5;
subsiz = 100;
subshf = valshf(:,randperm(size(valshf,2),subsiz));
V = violinplot([valsco(isant,:);valsco(ispost,:);valsco(:);subshf(:)], ...
    [ones(sum(isant),1);2*ones(sum(ispost),1);...
    3*ones(size(valsco,1),1);4*ones(size(subshf,1)*size(subshf,2),1)], ...
    'BandWidth',bandwidth);


V(1).ViolinPlot.FaceColor = [255,0,114]./256;
V(1).ScatterPlot.MarkerFaceColor = [255,0,114]./256;
V(2).ViolinPlot.FaceColor =[0,185,246]./256;
V(2).ScatterPlot.MarkerFaceColor = [0,185,246]./256;
V(3).ViolinPlot.FaceColor = [75,0,130]./256;
V(3).ScatterPlot.MarkerFaceColor = [75,0,130]./256;
V(4).ViolinPlot.FaceColor = 'none';
V(4).ScatterPlot.MarkerFaceColor = [0.5,0.5,0.5];


for i = 1:4
    V(i).ScatterPlot.MarkerFaceAlpha = 0.6;
    V(i).ViolinPlot.EdgeColor = 'none';
    V(i).BoxPlot.FaceColor = [0.2,0.2,0.2];
    V(i).BoxPlot.EdgeColor = [0.2,0.2,0.2];
    V(i).BoxPlot.LineWidth = 5;
end
V(4).ScatterPlot.MarkerFaceAlpha = 0.05;

ylim([-6 6]);  ylabel('signed valence score')
xlim([0,5]); xtickangle(-45);
set(gca,'xticklabel',{'ant.','post.','pooled','shuffled'},'box','off','tickdir','out',...
    'FontName','Arial','FontSize',15,'ytick',-8:2:8,'TickLength',[0.02, 0.02])