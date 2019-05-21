function [canbesolved deadendmets] = MissingTransDeadEnd(model)

% This function is to detect whether the deadend metabolites can be solved by adding a transport reaction.
% Output: canbesolved is a list of deadend metabolites that can be solved by adding a transport rxn.
%         deadendmets is a list of deadend metaolites
%         output format: deadentMetsName MNXID deadendMetsInanotherComps Possible solution
% Feiran Li 2019-02-01

changeCobraSolver('gurobi', 'LP');
exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = -1000;
model.ub(exchangeRxns) = 1000;
outputMets_Idx = detectDeadEnds(model);
outputMets = model.metNames(outputMets_Idx);
outputMetsMNX = model.metMetaNetXID(outputMets_Idx);

model_r = ravenCobraWrapper(model);

deadendmets = {'deadentMetsName','MNXID','deadendMetsInanotherComps','Possible solution'};
% deadendmets is a list contains all deadend mets in the model
deadendmets1 = {'deadentMetsName','MNXID','deadendMetsInanotherComps','Possible solution'};
% This one refers to deadend mets that may be solved by adding a transport rxn or change lb, which will be further checked in the later part of this code.
for i = 1:length(outputMets)
    mets_Idx = find(strcmp(model_r.metNames,model_r.metNames(outputMets_Idx(i))));
    metsinOtherComps_Idx = setdiff(mets_Idx,outputMets_Idx(i));
    metinOtherExE_Idx = metsinOtherComps_Idx(model_r.metComps(metsinOtherComps_Idx) ~= 3);% 3 refers to extracelluar compartment in model_r.metComps
    if ~isempty(metinOtherExE_Idx)
        for j = 1:length(metinOtherExE_Idx)
            if isempty(find(ismember(metinOtherExE_Idx(j),outputMets_Idx), 1))
                transrxn = intersect(find(model.S(outputMets_Idx,:)~=0),find(model.S(metinOtherExE_Idx(j),:)~=0));
                if ~isempty(transrxn)
                    for m = 1:length(transrxn)
                        if model.lb(transrxn(m)) == 0
                            deadendmets1 = [deadendmets1;outputMets(i),outputMetsMNX(i),model.metNames(metinOtherExE_Idx(j)), ['change lb for rxn ',model.rxns{transrxn(m)}]];
                        elseif model.lb(transrxn) == -1000
                            deadendmets = [deadendmets;outputMets(i),outputMetsMNX(i),model.metNames(metinOtherExE_Idx(j)), ['has transport reaction ', model.rxns{transrxn(m)},' for ', model.metNames{metinOtherExE_Idx(j)}, ' but not work']];
                        end
                    end
                else
                    if ~isempty(find(model_r.metComps(metinOtherExE_Idx(j)) == 1))
                        deadendmets1 = [deadendmets1;outputMets(i),outputMetsMNX(i),model.metNames(metinOtherExE_Idx(j)), ['add a transport reaction from cytosol for ',model.metNames{metinOtherExE_Idx(j)}]];
                    else
                        deadendmets1 = [deadendmets1;outputMets(i),outputMetsMNX(i),model.metNames(metinOtherExE_Idx(j)), ['add a transport reaction for ',model.metNames{metinOtherExE_Idx(j)}]];
                    end
                end
            else
                deadendmets = [deadendmets;outputMets(i),outputMetsMNX(i),model.metNames(metinOtherExE_Idx(j)),' mets in other compartment is also deadend'];
            end
        end
    else
        deadendmets = [deadendmets;outputMets(i),outputMetsMNX(i),'noMetsInOtherComps','mets only appear in one compartment'];
    end
end

% trying to add a transport reaction to see whether the deadend metaolites
% canbesolved is a list that deadend mets can be linked into the model by adding a transport rxns
canbesolved = deadendmets1(1,:); % title line
for i = 2:length(deadendmets1(:,1))% There is a title line
    if ~isempty(cell2mat(deadendmets1(i,1)))
        if strncmpi(deadendmets1(i,4),'add a transport reaction from cytosol',37)
            mets = [deadendmets1(i,1),deadendmets1(i,3)];
            [~,metindex] = ismember(mets,model.metNames);
            metsID = model.mets(metindex);
            cd ../otherChanges
            newID    = getNewIndex(model.rxns);
            cd ../modelCuration/
            TransRxn  = ['r_' newID];
            newModel = addReaction(model,TransRxn, ...
                'reactionName', [deadendmets1{i,2}, ' transport'], ...
                'metaboliteList', metsID, 'stoichCoeffList', [-1 1], ...
                'lowerBound', -1000, 'upperBound', 1000, 'subSystem', '', ...
                'checkDuplicate', false);
            outputMets_Idx_temp = detectDeadEnds(newModel);
            if ~ismember(deadendmets1{i,1},model.metNames(outputMets_Idx_temp)) % deadmet is not in the list
                canbesolved = [canbesolved;deadendmets1(i,:)];
                model = newModel;
            end
        elseif strncmpi(deadendmets1(i,4),'change lb',9)
            rxn_temp = deadendmets1{i,4};
            rxnID = rxn_temp(end-5:end);
            [~,rxnindex] = ismember({rxnID},model.rxns);
            newModel = model;
            newModel.lb(rxnindex) = -1000;
            outputMets_Idx_temp = detectDeadEnds(newModel);
           if ~ismember(deadendmets1{i,1},model.metNames(outputMets_Idx_temp))
                canbesolved = [canbesolved;deadendmets1(i,:)];
                model = newModel;
           end
        end
    end
end

end
