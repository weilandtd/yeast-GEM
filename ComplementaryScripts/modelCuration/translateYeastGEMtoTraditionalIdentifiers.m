function modelOut = translateYeastGEMtoTraditionalIdentifiers(model)
% This function is to translate the model to tradionally used identifiers 
% in the consensus model (s_xxxx[compartment] for metabolites and r_xxxx 
% for reactions)
% Inputs:  
%           model: yeast model
%
% author: Sebastián N. Mendoza 15/03/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelOut = model;
modelOut.mets = modelOut.metSIDs;
modelOut.rxns = modelOut.rxnRIDs;

end