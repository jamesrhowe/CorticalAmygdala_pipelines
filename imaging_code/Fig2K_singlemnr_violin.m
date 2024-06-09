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
bandwidth = 0.05; malpha1 = 0.25; malpha2 = 0.02;
msize = 25;

panelK1 = figure('Position',[0,0,300,500]);

v0 = violinplot([mean(subacc,2);mean(shf_subacc,2)],...
    [ones(1,length(subacc),1),2*ones(1,length(shf_subacc),1)],'BandWidth',bandwidth);
xtickangle(-45);
xlim([0,3]); ylim([0,0.5]); ylabel('MNR accuracies');
set(gca,'xtick',[1,2],'xticklabel',["data","shuffled"],...
    'ytick',[0,0.25,0.5,0.75,1],'FontName','Arial','FontSize',15,...
    'TickLength',[0.02, 0.02],'tickdir','out','box','off')

for a = 1:2
    v0(a).ScatterPlot.MarkerFaceAlpha = malpha1; 
    v0(a).ScatterPlot.SizeData = msize;
end

v0(1).ScatterPlot.MarkerFaceColor = [0.01,0.01,0.01];
v0(2).ScatterPlot.MarkerFaceColor = [0.8,0.8,0.8];
v0(1).ViolinPlot.FaceColor = 'none';
v0(2).ViolinPlot.FaceColor = 'none';

[h,p] = kstest2(mean(subacc,2),mean(shf_subacc,2));

title(strcat("kstest2 p=",num2str(p)))
%%
%now look across biological replicates
alist = unique(animid);
for a = 1:length(alist)
    sig_acc(a) = mean(mean(subacc(animid==a,:),2)>...
    prctile(mean(shf_subacc,2),95));
end
panelK2 = figure('Position',[0,0,300,500]);

b1 = bar(1,mean(sig_acc,2)); hold on
b1.FaceColor = 'none'; b1.LineWidth = 1.5;
e1 = errorbar(b1.XData,b1.YData,std(sig_acc,[],2)/sqrt(size(sig_acc,2)),...
    'LineStyle','none','Color','k','CapSize',8);
s1 = scatter(ones(size(sig_acc)),sig_acc,50,'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',0.7,'MarkerFaceColor',[0.25,0.25,0.25]);

ylim([-0.05 1])
ylabel('% of MNR acc. > 95^t^h shuffled','Interpreter','tex');
set(gca,'xtick',[],'box','off','tickdir','out',...
    'FontName','Arial','FontSize',15,'ytick',0:0.2:1,...
    'yticklabel',abs(0:20:100),'TickLength',[0.02, 0.02])