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
panelS3D = figure('Position',[0,0,350,500]);
b1 = bar(1:6,100*mean(responsive_by_animal(:,[1,5,3,6,2,4],:),1));
hold on
sem = std(responsive_by_animal(:,[1,5,3,6,2,4],:),[],1)./...
    sqrt(size(responsive_by_animal,1));
e2n = errorbar(b1.XData,b1.YData,100*sem,[],'Color','w','LineStyle','none');
e2p = errorbar(b1.XData,b1.YData,[],100*sem,'Color','k','LineStyle','none');
ylabel("% of neurons");
b1.FaceColor = 'flat';
b1.CData = [0.33,0,0;0.33,0,0;0.5,0.5,0;0.5,0.5,0;0,0,0.5;0,0,0.5];

set(gca,'ytick',0:10:80,'xticklabel',{"4MT","TMT","IAA","Heptanol",...
    "Peanut Oil","2PE"},'box','off','FontName','Arial','FontSize',15,...
    'tickdir','out','TickLength',[0.02, 0.02])

Y = responsive_by_animal(:,[1,5,3,6,2,4]);
X = repmat(["V","V","N","N","P","P"],size(responsive_by_animal,1),1);
[~,table0] = anovan(Y(:),{X(:)});

title(strcat('F=',num2str(table0{2,6}),',p=',num2str(table0{2,7})))