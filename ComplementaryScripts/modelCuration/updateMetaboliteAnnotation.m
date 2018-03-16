%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = updateMetaboliteAnnotation(model)
% update the metabolite annotation information in the model
%
% Hongzhong Lu & Benjamín J. Sánchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = updateMetaboliteAnnotation(model)

%Load data:
metData = tdfread('../../ComplementaryData/metabolite_manual_curation.tsv','\t');

for i = 1:length(metData)
    for j = 1:length(model.mets)
        if startsWith(model.metNames{j},[metData.name_original(i,:) ' ['])	%old name
            
            %find corresponding compartment:
            metName = model.metNames{j};
            comp    = metName(strfind(metName,' ['):strfind(metName,']'));
            metName = [metData.name_new(i,:) comp];
            
            %Update other fields:
            model.metNames{j}   = metName;
            model.metChEBIID{j} = metData.CHEBI_new(i,:);
            model.metKEGGID{j}  = metData.KEGG_new(i,:);
            model.metCharges(j) = metData.charge_new(i);
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%