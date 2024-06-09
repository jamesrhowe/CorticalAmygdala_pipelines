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
    % train classifiers on the last second of odor
    tstart = 71; tend = 75;
    
    neurlist = [1,5,10,20,30,50,70,100,150,200,250,300];
    tic
    for n = 1:length(neurlist)
        numneur = neurlist(n);
        for s = 1:100        
            subF3 = pooled_F3(randperm(totneur,numneur),:,:);
            
            pred = permute(mean(subF3(:,tstart:tend,:),2),[3,1,2]);
            
            % Train multinomial SVM classifier
            svm_model = fitcecoc(pred, resp);
            shf_model = fitcecoc(pred, resp(randperm(length(resp))));
            
            % Cross-validation example (optional)
            % Perform k-fold cross-validation to evaluate the model
            cv_model = crossval(svm_model,'KFold', k);
            cv_shf = crossval(shf_model,'KFold',k);
            
            % Calculate the k-fold loss (misclassification rate)
            acc(n,s) = 1-kfoldLoss(cv_model);
            shfacc(n,s) = 1-kfoldLoss(cv_shf);
        end
        disp([num2str(n),' numneur values out of ',num2str(length(neurlist)),' completed.'])
        toc

        save(strcat(currentScriptPath,"\data\ecoc_vs_neur_new.mat"),...
            'acc','shfacc','neurlist');
    end
else
    load('ecoc_vs_neur.mat','acc','shfacc','neurlist');
end
%%
panelN = figure('Position',[0,0,400,500]);
p1 = plot(neurlist(1:end),mean(acc(1:end,:),2),'ko-');
hold on
p2 = plot(neurlist(1:end),mean(shfacc(1:end,:),2),'o--','Color',[0.5,0.5,0.5]);
p1.LineWidth=1.5; p2.LineWidth=1.5; 
p1.MarkerFaceColor = p1.Color; p2.MarkerFaceColor = p2.Color;
e1 = errorbar(neurlist(1:end),mean(acc(1:end,:),2),std(acc(1:end,:),[],2));
e2 = errorbar(neurlist(1:end),mean(shfacc(1:end,:),2),std(shfacc(1:end,:),[],2));
e1.Color = p1.Color; e2.Color = p2.Color;
e1.LineStyle = 'none'; e2.LineStyle = 'none';
e1.CapSize=8; e2.CapSize=8;

leg = legend([p1,p2],{"data","shuffled"},'box','off');
x1 = 10; x2= 4;
leg.ItemTokenSize = [x1,x2]; leg.FontSize = 13;

xlim([-10 305]); ylim([0 1]);
xlabel('Number of Neurons'); ylabel('CV Accuracy of SVM during odor')
set(gca,'box','off','tickdir','out','FontName','Arial','FontSize',15,...
    'xtick',[1,50:50:300],'TickLength',[0.02,0.02])