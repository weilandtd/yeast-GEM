function [canbesolved]= MissingTransDeadEnd(model)

%This function is to detect whether the deadend metabolites can be solved
%by adding a transport reaction.
% Feiran Li 2019-02-01



solverOK = changeCobraSolver ('gurobi', 'LP');
exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = -1000;
model.ub(exchangeRxns) = 1000;
outputMets_index = detectDeadEnds(model);
outputMets = model.metNames(outputMets_index);
outputMetsMNX = model.metMetaNetXID(outputMets_index);
% * Detect reactions contained deadEnds metabolites
%[rxnList, rxnFormulaList] = findRxnsFromMets(model, outputMets);

model_r = ravenCobraWrapper(model);
for j = 1:length(model_r.metComps)
    metE = find(model_r.metComps == 3);
end
for j = 1:length(model_r.metComps)
    metC = find(model_r.metComps == 1);
end
deadendmets = [];
deadendmets1 = [];
for i = 1:length(outputMets)
    [metsindiffcomps] = metsinComps(model,model_r,outputMets(i));
    metsinOtherComps = setdiff(metsindiffcomps,outputMets(i));
    metinOtherExE = setdiff(metsinOtherComps,model.metNames(metE));
    if ~isempty(metinOtherExE)
        for j = 1:length(metinOtherExE)
            if isempty(find(strcmpi(metinOtherExE(j),outputMets)))
                [~,metOtherindex] = ismember(metinOtherExE(j),model.metNames);
                [~,metindex] = ismember(outputMets(i),model.metNames);
                transrxn = intersect(find(model.S(metindex,:)~=0),find(model.S(metOtherindex,:)~=0));
                if ~isempty(transrxn)
                    for m = 1:length(transrxn)
                        if model.lb(transrxn(m)) == 0
                            deadendmets1 = [deadendmets1;outputMets_index(i),outputMets(i),outputMetsMNX(i),metinOtherExE{j}, ['change lb for rxn ',model.rxns{transrxn(m)}]];
                        elseif model.lb(transrxn) == -1000
                            deadendmets = [deadendmets;outputMets_index(i),outputMets(i),outputMetsMNX(i),metinOtherExE{j}, ['has transport reaction ', model.rxns{transrxn(m)},' for ', metinOtherExE{j}, ' but not work']];
                        end
                    end
                else
                    if ~isempty(find(strcmp(metinOtherExE(j),model.metNames(metC))))
                        deadendmets1 = [deadendmets1;outputMets_index(i),outputMets(i),outputMetsMNX(i),metinOtherExE{j}, ['add a transport reaction from cytosol for ',metinOtherExE{j}]];
                    else
                        deadendmets1 = [deadendmets1;outputMets_index(i),outputMets(i),outputMetsMNX(i),metinOtherExE{j}, ['add a transport reaction for ',metinOtherExE{j}]];
                    end
                end
            else
                deadendmets = [deadendmets;outputMets_index(i),outputMets(i),outputMetsMNX(i),metinOtherExE{j},[' mets in other compartment is also deadend']];
            end
        end
    else
        deadendmets = [deadendmets;outputMets_index(i),outputMets(i),outputMetsMNX(i),'noMetsInOtherComps','mets only appear in one compartment'];     
    end
end

% trying to add a transport reaction to see whether the deadend metaolites
% can be solved
canbesolved = [];
for i = 1:length(deadendmets1(:,1))
        if ~isempty(cell2mat(deadendmets1(i,3)))
            if strncmpi(deadendmets1(i,5),'add a transport reaction from cytosol',37)
                mets = [deadendmets1(i,2),deadendmets1(i,4)];
                [~,metindex] = ismember(mets,model.metNames);
                metsID = model.mets(metindex);
                cd ../otherChanges
                newID    = getNewIndex(model.rxns);
                TransRxn  = ['r_' newID];
                newModel = addReaction(model,TransRxn, ...
                'reactionName', [deadendmets1{i,2}, ' transport'], ...
                'metaboliteList', metsID, 'stoichCoeffList', [-1 1], ...
                'lowerBound', -1000, 'upperBound', 1000, 'subSystem', '', ...
                'checkDuplicate', false);
                outputMets_index_temp = detectDeadEnds(newModel);
                %if ~isempty(setdiff(outputMets_index,outputMets_index_temp))
                if ~any(outputMets_index_temp(:)==deadendmets1{i,1})
                    canbesolved = [canbesolved;deadendmets1(i,:)];
                    model = newModel;
                end
            elseif strncmpi(deadendmets1(i,5),'change lb',9)
                rxn_temp = deadendmets1{i,5};
                rxnID = rxn_temp(end-5:end);
                [~,rxnindex] = ismember({rxnID},model.rxns);
                newModel = model;
                newModel.lb(rxnindex) = -1000;
                outputMets_index_temp = detectDeadEnds(newModel);
                if ~any(outputMets_index_temp(:)==deadendmets1{i,1})
                    canbesolved = [canbesolved;deadendmets1(i,:)];
                    model = newModel;
                end 
            end
        end
end

                