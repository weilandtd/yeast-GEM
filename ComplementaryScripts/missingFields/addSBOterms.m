%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addSBOterms(model)
%
% Benjamín J. Sánchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addSBOterms(model)

%Get RAVEN model for matching names & compartments
model_r = ravenCobraWrapper(model);

%Add SBO terms for mets:
model.metSBOTerms = cell(size(model.mets));
for i = 1:length(model.mets)
    metName = model_r.metNames{i};
    if ismember(metName,{'biomass','DNA','RNA','protein','carbohydrate','lipid'}) ...
            || endsWith(metName,' backbone') || endsWith(metName,' chain')
        model.metSBOTerms{i} = 'SBO:0000649';     %Biomass
    else
        model.metSBOTerms{i} = 'SBO:0000247';     %Simple chemical
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
