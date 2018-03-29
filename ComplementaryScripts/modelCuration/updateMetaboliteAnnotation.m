%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = updateMetaboliteAnnotation(model)
% update the metabolite annotation information in the model
%
% Hongzhong Lu & Benjamín J. Sánchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = updateMetaboliteAnnotation(model)

%Load data:
fid = fopen('../../ComplementaryData/metabolite_manual_curation.tsv','r');
metaboliteData = textscan(fid,'%s %s %s %s %f32 %s %s %s %f32 %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

for i = 1:length(metaboliteData{1})
    for j = 1:length(model.mets)
        if startsWith(model.metNames{j},[metaboliteData{2}{i} ' ['])	%old name
            %find corresponding compartment:
            metName = model.metNames{j};
            comp    = metName(strfind(metName,' ['):strfind(metName,']'));
            metName = [metaboliteData{6}{i} comp];
            
            %Update other fields:
            model.metNames{j}   = metName;              %new name
            model.metChEBIID{j} = metaboliteData{7}{i};	%new CHEBI
            model.metKEGGID{j}  = metaboliteData{8}{i};	%new KEGG
            model.metCharges(j) = metaboliteData{9}(i);	%new charge
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%