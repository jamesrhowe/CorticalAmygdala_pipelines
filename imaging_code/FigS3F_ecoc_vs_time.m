% current path:
currentScriptPath = fileparts(mfilename('fullpath'));
% generate string for all subdirectories:
pathString = genpath(currentScriptPath);
% add all subdirectories to path
addpath(pathString);

% Load F4 data:
load('F4.mat')
%%
% as a default, the code will analyze the results of ecoc accuracies and
% confusions calculated previously. if you would like to rerun the ecoc
% training and calculations, set rerun_ecoc = true;
rerun_ecoc = false;
%%
if rerun_ecoc == true
    k=10; totneur = size(pooled_F3,1);
    resp = reshape(repmat(1:6,20,1),[],1);
    numneur = 200;
    
    tic
    for t = 3:size(pooled_F3,2)
        for s = 1:100
            subF3 = pooled_F3(randperm(totneur,numneur),:,:);        
            pred = permute(mean(subF3(:,(t-2):t,:),2),[3,1,2]);
            
            % Train multinomial SVM classifier
            svm_model = fitcecoc(pred, resp);
            shf_model = fitcecoc(pred, resp(randperm(length(resp))));
            
            % Cross-validation example (optional)
            % Perform k-fold cross-validation to evaluate the model
            cv_model = crossval(svm_model,'KFold', k);
            cv_shf = crossval(shf_model,'KFold',k);
            
            % Calculate the k-fold loss (misclassification rate)
            acc(t,s) = 1-kfoldLoss(cv_model);
            shfacc(t,s) = 1-kfoldLoss(cv_shf);
        end
        disp([num2str(t),' frames out of ',num2str(size(pooled_F3,2)),' completed.'])
        toc

        save(strcat(currentScriptPath,"\data\ecoc_vs_time_new.mat"),...
            'acc','shfacc','neurlist');
    end
else
    load('ecoc_vs_time.mat','acc','shfacc','neurlist');
end
%%
panelS3F = figure('Position',[0,0,400,500]);
hold on
m1 = movmean(mean(acc,2),[2,0],1);
s1 = movmean(std(acc,[],2),[2,0],1);
p1 = plot(m1,'-','LineWidth',1.5);
patch([1:length(m1),fliplr(1:length(m1))],...
    [(m1-s1)',fliplr((m1+s1)')],p1.Color,...
    'FaceAlpha',0.5,'EdgeColor','none');

m2 = movmean(mean(shfacc,2),[2,0],1);
s2 = movmean(std(shfacc,[],2),[2,0],1);
p2 = plot(m2,'k--','LineWidth',1.5);
patch([1:length(m2),fliplr(1:length(m2))],...
    [(m2-s2)',fliplr((m2+s2)')],[0.5,0.5,0.5],'FaceAlpha',0.5,'EdgeColor','none');

xlim([25 125.5]); ylim([0 1]);
xlabel('time to odor (s)'); ylabel('CV Accuracy of SVM (k=200)')
set(gca,'box','off','tickdir','out','FontName','Arial',...
    'FontSize',15,'xtick',25.5:25:125.5,'xticklabel',-5:5:20,...
    'TickLength',[0.02,0.02])