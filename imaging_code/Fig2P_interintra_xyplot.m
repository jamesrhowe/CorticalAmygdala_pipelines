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
    
    tstart = 71; tend = 75;
    
    neurlist = [1,5,10,20,30,50,70,100,150,200,250,300];
    
    for n = 1:length(neurlist)
        nStart = tic; 
        numneur = neurlist(n);
        for s = 1:100        
            tic
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
    
            oofLabel = kfoldPredict(cv_model);
            ConfMat = confusionchart(resp,oofLabel,'RowSummary','total-normalized');
    
            conf_intrav(n,s)=mean([ConfMat.NormalizedValues(1,5),ConfMat.NormalizedValues(5,1),...
                ConfMat.NormalizedValues(2,4),ConfMat.NormalizedValues(4,2)])/20;
            conf_interv(n,s)=mean([...
                ConfMat.NormalizedValues(1,2),ConfMat.NormalizedValues(2,1),...
                ConfMat.NormalizedValues(5,2),ConfMat.NormalizedValues(2,5),...
                ConfMat.NormalizedValues(1,4),ConfMat.NormalizedValues(4,1),...
                ConfMat.NormalizedValues(5,4),ConfMat.NormalizedValues(4,5)])/20;
    
            shf_oofLabel = kfoldPredict(cv_shf);
            shf_ConfMat = confusionchart(resp,shf_oofLabel,'RowSummary','total-normalized');
    
            shf_conf_intrav(n,s)=mean([shf_ConfMat.NormalizedValues(1,5),shf_ConfMat.NormalizedValues(5,1),...
                shf_ConfMat.NormalizedValues(2,4),shf_ConfMat.NormalizedValues(4,2)])/20;
            shf_conf_interv(n,s)=mean([...
                shf_ConfMat.NormalizedValues(1,2),shf_ConfMat.NormalizedValues(2,1),...
                shf_ConfMat.NormalizedValues(5,2),shf_ConfMat.NormalizedValues(2,5),...
                shf_ConfMat.NormalizedValues(1,4),shf_ConfMat.NormalizedValues(4,1),...
                shf_ConfMat.NormalizedValues(5,4),shf_ConfMat.NormalizedValues(4,5)])/20;
            toc
        end
        disp([num2str(n),' numneur values out of ',num2str(length(neurlist)),' completed.'])
        save(strcat(currentScriptPath,"\data\intrainter_vs_neur_new.mat"),...
            'acc','shfacc','neurlist','shf_conf_intrav','shf_conf_interv',...
            'conf_intrav','conf_interv');
        tEnd = toc(nStart)
    end
else
    load('intrainter_vs_neur.mat',...
    'acc','shfacc','neurlist','shf_conf_intrav','shf_conf_interv',...
    'conf_intrav','conf_interv');
end
%%
lw1 = 1.2; lw2 = 1; cs = 8; cs2 = 6;
cintra = [0.5,0.5,0.5]; cinter = [1;0.05;0.05];
cshf = [0;0;0];

panelP = figure('Position',[0,0,450,500]); hold on

p1 = plot(neurlist+1,mean(conf_interv,2),'o-',...
    'LineWidth',lw1,'Color',cinter,'MarkerFaceColor',cinter);
p2 = plot(neurlist,mean(conf_intrav,2),'o-',...
    'LineWidth',lw1,'Color',cintra,'MarkerFaceColor',cintra);
p3 = plot(neurlist+1,mean(shf_conf_interv,2),'x-',...
    'LineWidth',lw1,'Color',cshf,'MarkerSize',18);
p4 = plot(neurlist,mean(shf_conf_intrav,2),'o--',...
    'LineWidth',lw1,'Color',cshf);

e1 = errorbar(p1.XData,p1.YData,std(conf_interv,[],2),...
    'LineStyle','none','CapSize',cs,'Color',cinter,'LineWidth',lw2);
e2 = errorbar(p2.XData,p2.YData,std(conf_intrav,[],2),...
    'LineStyle','none','CapSize',cs,'Color',cintra,'LineWidth',lw2);
e3 = errorbar(p3.XData,p3.YData,std(shf_conf_interv,[],2),...
    'LineStyle','none','CapSize',cs2,'Color',cshf,'LineWidth',lw2);
e4 = errorbar(p4.XData,p4.YData,std(shf_conf_intrav,[],2),...
    'LineStyle','none','CapSize',cs2,'Color',cshf,'LineWidth',lw2);

leg = legend([p1,p2,p3,p4],...
    {"inter-valence","intra-valence","inter shuffle","intra shuffle"},'box','off');
x1 = 20; x2= 4;
leg.ItemTokenSize = [x1,x2]; leg.FontSize = 13;


xlim([-10 310]); ylim([0 0.3]);
xlabel('number of neurons'); ylabel('average confusion')
set(gca,'box','off','tickdir','out','FontName','Arial','FontSize',15,...
    'xtick',[1,50:50:300],'TickLength',[0.02,0.02])