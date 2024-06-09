% current path:
currentScriptPath = fileparts(mfilename('fullpath'));
% generate string for all subdirectories:
pathString = genpath(currentScriptPath);
% add all subdirectories to path
addpath(pathString);

% Load F4 data:
load('F4.mat')
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
% check which animals were anterior or posterior
check_ant = @(x) strcmp(aloc{x}, 'a');
check_post = @(x) strcmp(aloc{x}, 'p');

% apply cellfun from above to animid:
isant = cellfun(check_ant, num2cell(animid));
ispost = cellfun(check_post, num2cell(animid));
%%
trialbins = [1:16;5:20];
inter_idx = [1,2;1,4;5,2;5,4];
intra_idx = [1,5;2,4;3,6];

%do the analysis on the last second of the odor:
odbins = 71:75;

%now calculate the valence score defined as abs(Fnut-F2pe) - abs(Fpos-Fneg)
%across different bins for trials:

%we should look at this across animals!!
alist = unique(neuron_to_animal);

%loop through every animal
for a = 1:length(alist)
    %and every timebin
    for t = 1:size(trialbins,2)
        clear Fsnap; clear inter_diff; clear intra_diff; clear valsco;
        clear sinter_diff; clear sintra_diff; clear valshf;
    
        Fsnap = permute(mean(mean(pooled_F4(neuron_to_animal==alist(a),odbins,...
            trialbins(1,t):trialbins(2,t),:),2),3),[1,4,2,3]);
    
        for i = 1:size(inter_idx,1)
            inter_diff(:,i) = (Fsnap(:,inter_idx(i,2)) - Fsnap(:,inter_idx(i,1)));
        end
        
        for j = 1:size(intra_idx,1)
            intra_diff(:,j) = abs(Fsnap(:,intra_idx(j,1)) - Fsnap(:,intra_idx(j,2)));
        end
        
        valsco = mean(inter_diff,2)./mean(intra_diff,2);
        
        idx_perm = perms(1:6); %there's 720 combinations total
        for s = 1:size(idx_perm,1)
            Shfsnap = Fsnap(:,idx_perm(s,:));
            for i = 1:size(inter_idx,1)
                sinter_diff(:,i,s) = (Shfsnap(:,inter_idx(i,2)) - Shfsnap(:,inter_idx(i,1)));
            end
        
            for j = 1:size(intra_idx,1)
                sintra_diff(:,j,s) = abs(Shfsnap(:,intra_idx(j,1)) - Shfsnap(:,intra_idx(j,2)));
            end
        
            valshf(:,s) = mean(sinter_diff(:,:,s),2)./mean(sintra_diff(:,:,s),2);
        end
    
    %3rd index = 1 is significant negative, = 2 is significant positive
    sign_val(a,t,1) = mean(valsco(:)<prctile(valshf(:),5));
    sign_val(a,t,2) = mean(valsco(:)>prctile(valshf(:),95));
    end
end

%structure the data for rmANOVA
table_neg = array2table(sign_val(:,:,1),'VariableNames', ...
    cellstr(strcat('T', string(1:size(sign_val, 2)))));
table_neg = addvars(table_neg,cellstr(aloc'),'Before','T1',...
    'NewVariableNames','lensloc');
trialtime = trialbins(1,:);
rm_neg = fitrm(table_neg,'T1-T16 ~ lensloc','WithinDesign',trialtime);
ranovatbl_neg = ranova(rm_neg);
%the rmANOVA fails to show significant interaction of time
%and lensloc
%%
%now to plot, for every point along t, we will graph mean significant
%positive valence neurons above 0 and mean significant negative valence
%neurons below 0 with sem across biological replicates. for anterior
%imaging animals, we will plot in magenta and posterior imaging animals we
%will plot in cyan:
panelJ = figure('Position',[0,0,500,500]);

%pinkish color for anterior, blueish color for posterior:
aC = [255,0,114]./256;
pC = [0,185,246]./256;
hold on
p1 = plot(trialbins(1,:),100*mean(sign_val(aloc=="a",:,2),1));
p2 = plot(trialbins(1,:),-100*mean(sign_val(aloc=="a",:,1),1));
p3 = plot(trialbins(1,:)+0.2,100*mean(sign_val(aloc=="p",:,2),1));
p4 = plot(trialbins(1,:)+0.2,-100*mean(sign_val(aloc=="p",:,1),1));

p1.Color = aC; p2.Color = aC; p3.Color = pC; p4.Color = pC;
p1.Marker = 'o'; p2.Marker = 'o'; p3.Marker = 'o'; p4.Marker = 'o';
p1.MarkerFaceColor = aC; p2.MarkerFaceColor = aC;
p3.MarkerFaceColor = pC; p4.MarkerFaceColor = pC;
p1.LineWidth = 1.5; p2.LineWidth = 1.8; 
p3.LineWidth = 1.5; p4.LineWidth = 1.8;
p2.LineStyle = '--'; p4.LineStyle = '--';

s1 = errorbar(p1.XData,p1.YData,100*std(sign_val(aloc=="a",:,2),[],1)./...
    sqrt(sum(aloc=="a")),'LineStyle','none','Color',aC,...
    'CapSize',12,'LineWidth',1.5);
s2 = errorbar(p2.XData,p2.YData,100*std(sign_val(aloc=="a",:,1),[],1)./...
    sqrt(sum(aloc=="a")),'LineStyle','none','Color',aC,...
    'CapSize',12,'LineWidth',1.5);
s3 = errorbar(p3.XData,p3.YData,100*std(sign_val(aloc=="p",:,2),[],1)./...
    sqrt(sum(aloc=="p")),'LineStyle','none','Color',pC,...
    'CapSize',12,'LineWidth',1.5);
s4 = errorbar(p4.XData,p4.YData,100*std(sign_val(aloc=="p",:,1),[],1)./...
    sqrt(sum(aloc=="p")),'LineStyle','none','Color',pC,...
    'CapSize',12,'LineWidth',1.5);
yline(0,'k-'); xlim([0 17]);
ylabel('% with significant valence score'); xlabel('trial#')
set(gca,'xtick',[1,5,10,16],'box','off','tickdir','out',...
    'FontName','Arial','FontSize',15,'ytick',-20:5:20,...
    'yticklabel',abs(-20:5:20),'TickLength',[0.02, 0.02])

leg = legend([p1,p3],{"anterior","posterior"},'box','off');
x1 = 5; x2= 4;
leg.ItemTokenSize = [x1,x2]; leg.FontSize = 10;

title(strcat('F=',num2str(ranovatbl_neg{2,4}),',p=',num2str(ranovatbl_neg{2,5})))