% current path:
currentScriptPath = fileparts(mfilename('fullpath'));
% generate string for all subdirectories:
pathString = genpath(currentScriptPath);
% add all subdirectories to path
addpath(pathString);

% Load F4 data:
load('F4.mat')
%%
% as a default, the code will analyze the results of MNR accuracies and
% confusions calculated previously. if you would like to rerun the MNR
% training and calculations, set rerun_mnr = true;
rerun_mnr = false;

if rerun_mnr == true
    warning('off', 'all');
    resp = reshape(repmat(1:6,20,1),[],1);
    k = 5;
    % train classifier on the last second of odor
    tstart = 71; tend = 75;
    % as a control, compare against 100 trials shuffled data:
    shf_n = 100;
    
    % initialize some matrices:
    accuracy = zeros(size(pooled_F3,1),k);
    shf_acc = zeros(size(pooled_F3,1),shf_n,k);
    
        for n = 1:size(pooled_F3,1)
            tic
            pred = permute(mean(pooled_F3(n,tstart:tend,:),2),[3,1,2]);
            % Create a cross-validation partition
            c = cvpartition(resp, 'KFold', k);
        
            for i = 1:k
                trainingPredictors = pred(c.training(i),:);
                trainingResponse = categorical(resp(c.training(i)));
        
                testPredictors = pred(c.test(i),:);
                testResponse = categorical(resp(c.test(i),:));
                   
        
                for j = 1:size(pred,2)
                    % Predict on the test set
                    B{n,i} = mnrfit(trainingPredictors(:,j), categorical(trainingResponse),...
                            'model', 'nominal', 'interactions', 'on');
                    pihat = mnrval(B{n,i},testPredictors(:,j));
                    [~,predicted] = max(pihat,[],2);
        
                    accuracy(n,i) = mean(categorical(predicted) == testResponse);
        
                    confusion{n,i} = confusionmat(testResponse, categorical(predicted));
                        
                    for f = 1:shf_n
                        warning('off', 'all');
                        shf_train = trainingResponse(randperm(length(trainingResponse)));
                        shf_test = testResponse(randperm(length(testResponse)));
        
                        Bshf{n,f,i} = mnrfit(trainingPredictors(:,j), categorical(shf_train),...
                            'model', 'nominal', 'interactions', 'on');
        
                        % Predict on a shuffled set
                        sfhat = mnrval(Bshf{n,f,i}, testPredictors(:,j));
                        [~, sf_pred] = max(sfhat, [], 2);            
        
                        % Compute the shuffled accuracy for this fold
                        shf_acc(n,f,i) = mean(categorical(sf_pred) == shf_test);
        
                        % Shuffled confusion matrix:
                        shf_conf{n,f,i} = confusionmat(testResponse, categorical(sf_pred));
                    end
        
                end
            end
            toc
            save(strcat(currentScriptPath,"\data\single_mnr_new.mat"),...
                'B','Bshf','accuracy','shf_acc','shf_conf','confusion');
        end
else
    load("single_mnr.mat");
end    
%%
shf_conf_L = reshape(shf_conf,[size(shf_conf,1)*size(shf_conf,2),5]);
for N = 1:length(shf_conf_L)
    tempC = zeros(size(shf_conf_L{N,1}));
    for i = 1:size(shf_conf_L,2)
        tempC = tempC+shf_conf_L{N,i};
    end
    kfold_acc(N,:) = diag(tempC)./sum(tempC,2); %get accuracies for individual fold
    conf(:,:,N) = tempC;
    
    for a = 1:size(tempC,1) %also get accuracies by odor
        shf_subacc(N,a) = tempC(a,a)/sum(tempC(a,:),2);
    end

    sorted_shf_subacc(N,:) = sort(shf_subacc(N,:),'descend'); %also look at Bshf
end
%%
clear conf
for N = 1:length(confusion)
    tempC = zeros(size(confusion{N,1}));
    for i = 1:size(confusion,2)
        tempC = tempC+confusion{N,i};
    end
    kfold_acc(N,:) = diag(tempC)./sum(tempC,2); %get accuracies for individual fold
    conf(:,:,N) = tempC;
    
    for a = 1:size(tempC,1) %also get accuracies by odor
        subacc(N,a) = tempC(a,a)/sum(tempC(a,:),2);
    end

    sorted_subacc(N,:) = sort(subacc(N,:),'descend'); %also look at Bshf
end
%%
alist = unique(animid);
for a = 1:length(alist)
    for i = 1:size(sorted_subacc,2)
        sig_subacc(a,i) = mean(sorted_subacc(animid==a,i)>...
            prctile(sorted_shf_subacc(:,i),95));
    end
end

x = repmat(1:6,size(alist,1),1); x = x(:);
y = sig_subacc; y= y(:);

% Fit linear model
linearModel = fitlm(x, y);

% Fit exponential decay model
% The model form is y = a * exp(b * x)
expDecayModel = @(b, x) b(1) * exp(b(2) * x);
initialGuess = [y(1), -0.1]; % Initial guess for parameters [a, b]
expDecayParams = nlinfit(x, y, expDecayModel, initialGuess);

% Calculate predicted values and residuals for exponential model
y_pred_exp = expDecayModel(expDecayParams, x);
residuals_exp = y - y_pred_exp;
RSS_exp = sum(residuals_exp.^2);

% Number of data points and number of parameters
n = length(y); % Number of data points
p_linear = 2; % Number of parameters in linear model
p_exp = 2; % Number of parameters in exponential model

% AIC and BIC for linear model
linearAIC = linearModel.ModelCriterion.AIC;
linearBIC = linearModel.ModelCriterion.BIC;

% AIC and BIC for exponential decay model
expAIC = n * log(RSS_exp/n) + 2 * p_exp;
expBIC = n * log(RSS_exp/n) + log(n) * p_exp;

% Generate predictions
y_pred_linear = predict(linearModel, (1:6)');
y_pred_exp = expDecayModel(expDecayParams, (1:6)');
%%
panelM = figure('Position',[0,0,400,500]);
scatter(x,y,45,'k','filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
hold on
pl=plot(1:6,y_pred_linear,'x-','LineWidth',1.5,'MarkerSize',15);
pe=plot(1:6,y_pred_exp,'x-','LineWidth',1.5,'MarkerSize',15);
pd=plot(1:6,mean(sig_subacc,1),'-','LineWidth',1.5,'Color',[0.7,0.7,0.7]);

pl.Color=[255,0,114]./256;
pe.Color=[0,185,246]./256;

yline(0,'k-'); 

leg = legend([pd,pl,pe],{"data average",strcat("lin. fit (BIC= ",...
    num2str(linearBIC,'%.1f'),")"),strcat("exp. fit (BIC= ",...
    num2str(expBIC,'%.f'),")")},'box','off');
x1 = 10; x2= 4;
leg.ItemTokenSize = [x1,x2]; leg.FontSize = 13;

xlim([0 7]); ylim([-0.1 1])

ylabel('% subaccuracies > 95^t^h of shuffled','Interpreter','tex'); 
xlabel('sorted rank')
set(gca,'xtick',1:6,'box','off','tickdir','out',...
    'FontName','Arial','FontSize',15,'ytick',0:0.2:1,...
    'yticklabel',abs(0:20:100),'TickLength',[0.02, 0.02])
%%