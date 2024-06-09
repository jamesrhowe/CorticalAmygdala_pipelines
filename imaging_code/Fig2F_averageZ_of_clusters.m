% current path:
currentScriptPath = fileparts(mfilename('fullpath'));
% generate string for all subdirectories:
pathString = genpath(currentScriptPath);
% add all subdirectories to path
addpath(pathString);

% Load F3 data:
load('F3.mat')
%%
% set frame bins to average in to generate vectors for clustering
% (corresponds to 6s pre-odor, then 6x 2s bins from the onset of odor):
startbins = [21,51:10:101];
endbins = [50,55:10:105];

% set default options for clustering using clustheat_defaults() fn.
defaultOptions = clustheat_defaults();
options = defaultOptions; % and use default options
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
% pre-allocate matrix tosort which will have average F's from the
% pre-defined timebins:
tosort = zeros(Fsiz(1),length(startbins),Fsiz(3));

% loop through each odor
for o = 1:Fsiz(3)
    % and loop through each timebin
    for t = 1:length(startbins)
        % to find time-average during said timebin for said odor:
        tosort(:,t,o) = mean(normF3(:,startbins(t):endbins(t),o),2);
    end
end

% flatten tosort matrix into 2d before clustering:
tosort = reshape(tosort,[size(tosort,1),size(tosort,2)*size(tosort,3)]);

% then perform the clustering/reordering:
% 1. find pairwise distance using the distance metric specified in options.metric
D = pdist(tosort,options.metric); 
% 2. perform hierarchical clustering using method specified in options.method
tree = linkage(D,options.method);
% 3. order the fluorescence matrix such that adjacent neurons are most similar
leafOrder = optimalleaforder(tree,D);
% 4. find colorthreshold that matches the number of maximum clusters as
% defined by options.numclust:
cthresh = tree(end-options.numclust+2,3)-eps;
% 5. define clusters according to options.numclust:
T = cluster(tree,"maxclust",options.numclust);
%%
panelF = figure('Position',[0,0,400,700]);
ocolor = [1,0,0;1,0,0;0,0,1;0,0,1;0,1,0;0,1,0];
ymin = -3.5; ymax = 3.5;
xmin = 0; xmax = size(normF2,2);
ralpha = 0.3;

subplot(1,options.numclust,1);
ylabel({'avg Z-score';'of cluster'});
for c = 1:options.numclust
    subplot(options.numclust,1,c); hold on
    M = movmean(mean(normF2(T==c,:),1),[3,0],2);
    S = movmean(std(normF2(T==c,:),[],1),[3,0],2);

    plot(M,'k-','LineWidth',1.5);
    patch([1:length(normF2),fliplr(1:length(normF2))],...
        [M-S,fliplr(M+S)],[0.6,0.6,0.6],'FaceAlpha',0.5,'EdgeColor','none');

    ylim([ymin,ymax]); xlim([xmin,xmax]);
    for o = 1:6
        h=fill([(o-1)*Fsiz(2)+50.5,(o-1)*Fsiz(2)+75.5,...
            (o-1)*Fsiz(2)+75.5,(o-1)*Fsiz(2)+50.5],...
            [ymin,ymin,ymax,ymax],ocolor(o,:));
        set(h,'FaceAlpha',ralpha,'EdgeColor','none');
    end
    
    set(gca,'xtick',[],'ytick',-3:3:3,'box','off','tickdir','out',...
        'FontName','Arial','TickLength',[0.02, 0.02]);
end