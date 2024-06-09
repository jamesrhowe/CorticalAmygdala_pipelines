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
% check which animals were anterior or posterior
check_ant = @(x) strcmp(aloc{x}, 'a');
check_post = @(x) strcmp(aloc{x}, 'p');

% apply cellfun from above to animid:
isant = cellfun(check_ant, num2cell(animid));
ispost = cellfun(check_post, num2cell(animid));
%%
% average responsive_by_animal grouped by anterior vs posterior animals
mu_ap = [mean(responsive_by_animal(aloc=="a",:),1);...
    mean(responsive_by_animal(aloc=="p",:),1)];
% also calculate SEM
sem_ap = [std(responsive_by_animal(aloc=="a",:),[],1)./sqrt(sum(aloc=="a"));...
    std(responsive_by_animal(aloc=="p",:),[],1)./sqrt(sum(aloc=="p"))];
%%
aloc_transp = repmat(permute(aloc,[2,1]),1,6);
% specify valence (V for aversive, P for appetitive, N for neutral
val_transp = repmat(["V","P","N","P","V","N"],13,1);

% 2-way anova on lens location and valence
[~,table1] = anovan(responsive_by_animal(:),{aloc_transp(:),val_transp(:)},'varnames',{'a vs p','valence'},...
    'model','interaction');
%%
panelH = figure('Position',[0,0,450,500]);
hold on
b1 = bar([1,5,3,6,2,4]*3,100*mu_ap(1,:)); b1.BarWidth = 0.25;
b2 = bar([1,5,3,6,2,4]*3+1,100*mu_ap(2,:)); b2.BarWidth = 0.25;
e1n = errorbar(b1.XData,b1.YData,100*sem_ap(1,:),[],'Color','w','LineStyle','none');
e2n = errorbar(b2.XData,b2.YData,100*sem_ap(2,:),[],'Color','w','LineStyle','none');
e1p = errorbar(b1.XData,b1.YData,[],100*sem_ap(1,:),'Color','k','LineStyle','none');
e2p = errorbar(b2.XData,b2.YData,[],100*sem_ap(2,:),'Color','k','LineStyle','none');
title(strcat('F=',num2str(table1{4,6}),',p=',num2str(table1{4,7})))

xlim([1,21]); xtickangle(-45)
ylabel("% of neurons responsive");
b1.FaceColor = [255,0,114]./256;
b2.FaceColor = [0,185,246]./256;

set(gca,'ytick',0:10:80,'xtick',(1:6)*3+0.5,...
    'xticklabel',{"4MT","TMT","IAA","Heptanol",...
    "Peanut Oil","2PE"},'box','off','FontName','Arial','FontSize',15,...
    'tickdir','out','TickLength',[0.02, 0.02])

leg = legend([b1,b2],{"anterior","posterior"},'box','off');
x1 = 5; x2= 4;
leg.ItemTokenSize = [x1,x2]; leg.FontSize = 10;