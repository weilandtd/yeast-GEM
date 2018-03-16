%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = updateMetaboliteAnnotation(model)
% update the metabolite annoation information in the model
%
% Hongzhong Lu & Benjamín J. Sánchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = updateMetaboliteAnnotation(model)

%Load data:
fid = fopen('../../ComplementaryData/metabolite_manual_curation.csv');
metaboliteData = textscan(fid,'%s %s %s %s %f32 %s %s %s %f32 %s','Delimiter',',','HeaderLines',1);
fclose(fid);

for i = 1:length(metaboliteData)
    for j = 1:length(model.mets)
        if contains(model.mets{j},metaboliteData{1}{i})	%metID
            metName = model.metNames{j};
            comp    = metName(strfind(metName,' ['):strfind(metName,']'));
            metName = [metaboliteData{6}{i} comp];
            model.metNames{j}   = metName;              %new name
            model.metChEBIID{j} = metaboliteData{7}{i};	%new CHEBI
            model.metKEGGID{j}  = metaboliteData{8}{i};	%new KEGG
            model.metCharges(j) = metaboliteData{9}(i);	%new charge
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%