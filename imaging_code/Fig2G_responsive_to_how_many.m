% current path:
currentScriptPath = fileparts(mfilename('fullpath'));
% generate string for all subdirectories:
pathString = genpath(currentScriptPath);
% add all subdirectories to path
addpath(pathString);

% Load F4 data:
load('F4.mat')

% set default options for responsive period using clustheat_defaults() fn.
defaultOptions = clustheat_defaults();
options = defaultOptions; % and use default options
%%
%check responsiveness to each odor using fn "isresponsive"
for o = 1:6
    response_pval(:,:,o) = isresponsive(pooled_F4,o,options.pre_start,...
        options.pre_end,options.od_start,options.od_end);
end

%now we have to switch pval to null rejection with Holm-Bonferroni method
alf = 0.05/length(response_pval); %define cutoff
response_rej = response_pval<alf;
consecutive_rej = (diff(response_rej,1,2)==0 & response_rej(:,2:end,:)==1);

final_resp = reshape(sum(consecutive_rej,2) > 2,...
    [size(consecutive_rej,1),size(consecutive_rej,3)]);

for a = 1:length(unique(animid))
    responsive_by_animal(a,:) = mean(final_resp(animid==a,:),1);
end
%%
%count how many neurons are repsonsive to 0,1,2,..., etc. odors
responsive_to_how_many = sum(final_resp,2);

for a = 1:length(unique(animid))
    resp_distr_by_animal(a,:) = histcounts(sum(final_resp(animid==a,:),2),-0.5:1:6.5);
end
resp_distr_by_animal = resp_distr_by_animal./sum(resp_distr_by_animal,2);
%%
panelG = figure('Position',[0,0,350,500]);
b2 = bar(0:6,100*mean(resp_distr_by_animal,1));
b2.FaceColor = 'flat';
hold on
sem = std(resp_distr_by_animal,[],1)./...
    sqrt(size(resp_distr_by_animal,1));
e2n = errorbar(b2.XData,b2.YData,100*sem,[],'Color','w','LineStyle','none');
e2p = errorbar(b2.XData,b2.YData,[],100*sem,'Color','k','LineStyle','none');
ylabel("% of neurons"); xlabel("responsive to how many odors")

for i = 1:size(b2.CData,1)
    if mod(i,2) == 0
        b2.CData(i,:) = [0.1,0.1,0.1];
    else
        b2.CData(i,:) = [0.4,0.4,0.4];
    end
end
ylim([-1 60])
set(gca,'ytick',0:10:80,'xtick',0:6,'tickdir','out',...
    'box','off','FontName','Arial','FontSize',15,'TickLength',[0.02, 0.02])