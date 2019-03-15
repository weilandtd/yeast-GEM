function modelOut = translateYeastGEMtoBiGGIdentifiers(model, option)
% This Function is to translate the model to BIGG identifiers. In
% particular, it replaces the content of model.mets with BiGG IDs. 
% Inputs:  
%           model: yeast model
%           option: 1 is to translate just metabolites which are currently
%                   in the BiGG database. 2 is to also add suggested BiGG
%                   identifiers (not in BiGG) database.
%
% author: Sebastián N. Mendoza 15/03/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelOut = model;
switch option
    
    case 1
        %just translates those metabolites which are in the BiGG database
        if isfield(modelOut, 'metBiGGID')
            modelOut.mets(~cellfun(@isempty, modelOut.metBiGGID)) = modelOut.metBiGGID(~cellfun(@isempty, modelOut.metBiGGID));
        end
        
        if isfield(modelOut, 'rxnBiGGID')
            modelOut.rxns(~cellfun(@isempty, modelOut.rxnBiGGID)) = modelOut.rxnBiGGID(~cellfun(@isempty, modelOut.rxnBiGGID));
        end
    case 2
        %also assigns BiGG-like IDs to unmapped metabolites
        [~,s] = xlsread('../ComplementaryData/modelCuration/globalMetDictionary.csv');
        for i = 1:size(s,1)
            dict = splitString(s{i},',');
            modelOut.mets{strcmp(model.metSIDs,dict{1})==1} = dict{2};
        end
        
        [~,s] = xlsread('../ComplementaryData/modelCuration/globalRxnDictionary.csv');
        for i = 1:size(s,1)
            dict = splitString(s{i},',');
            modelOut.rxns{strcmp(model.rxnRIDs,dict{1})==1} = dict{2};
        end
end
end