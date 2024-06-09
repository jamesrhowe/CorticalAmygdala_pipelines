% current path:
currentScriptPath = fileparts(mfilename('fullpath'));
% generate string for all subdirectories:
pathString = genpath(currentScriptPath);
% add all subdirectories to path
addpath(pathString);

% Load F4 data:
load('F4.mat')
%%
k=10; totneur = size(pooled_F3,1);
resp = reshape(repmat(1:6,20,1),[],1);
% train classifiers on the last second of odor
tstart = 71; tend = 75;

numneur = 200; %for an example confusion matrix, train a svm ecoc on 200 neurons     
subF3 = pooled_F3(randperm(totneur,numneur),:,:);       
pred = permute(mean(subF3(:,tstart:tend,:),2),[3,1,2]);
        
% Train multinomial SVM classifier
svm_model = fitcecoc(pred, resp);
shf_model = fitcecoc(pred, resp(randperm(length(resp))));
        
% Cross-validation example (optional)
% Perform k-fold cross-validation to evaluate the model
cv_model = crossval(svm_model,'KFold', k);
cv_shf = crossval(shf_model,'KFold',k);

oofLabel = kfoldPredict(cv_model);
ConfMat = confusionchart(resp,oofLabel,'RowSummary','total-normalized');
%%
panelO = figure('Position',[0,0,520,500]);
% reorder the confusion matrix and normalize to total number of trials for
% each odor:
h = heatmap(ConfMat.NormalizedValues([1,5,3,6,2,4],[1,5,3,6,2,4])./20); 
colormap parula
h.XData = {"4MT","TMT","IAA","Heptanol","Peanut Oil","2PE"};
h.YData = {"4MT","TMT","IAA","Heptanol","Peanut Oil","2PE"};

h.CellLabelFormat = '%.2f';
set(gca,'FontName','Arial','FontSize',15);