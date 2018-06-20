%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addDBNewGeneAnnotation
% Add changes from the database new anootation for new genes + manual curation on those changes
% Input: model, databasenewGPR.tsv,SGDgeneNames.tsv.
As for the reference of new GPR, please find detailed information in the /ComplementaryData/databases/DBnewGeneAnnotation.tsv

% NOTE: changeGeneAssociation.m is a function from cobra
%
% Feiran Li & Hongzhong Lu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load model
model = readCbModel('../../ModelFiles/xml/yeastGEM.xml');
%Change GPR relations

fid3 = fopen('../../ComplementaryData/modelCuration/databasenewGPR.tsv');
changegpr = textscan(fid3,'%s %s %s','Delimiter','\t','HeaderLines',1);
newGPR.ID = changegpr{1};
newGPR.oldGPR = changegpr{2};
newGPR.GPR = changegpr{3};
fclose(fid3);

for i = 1:length(newGPR.ID)
RxnIndex = find(strcmp(model.rxns, newGPR.ID(i)));
model = changeGeneAssociation(model, model.rxns{RxnIndex}, newGPR.GPR{i});% for new reactions we add new GPR. For existing reactions, we update the GPR correlations.
end

% add gene standard name for new genes
%fid2 = fopen('../../ComplementaryData/yeast_gene_annotation_SGD.tsv');
fid2 = fopen('../../ComplementaryData/databases/SGDgeneNames.tsv');
yeast_gene_annotation = textscan(fid2,'%s %s','Delimiter','\t','HeaderLines',1);
fclose(fid2);


geneIndex = zeros(1,1);
for i = 1: length(model.genes)
    if ~isempty(find(strcmp(yeast_gene_annotation{1}, model.genes{i})))
        geneIndex(i) = find(strcmp(yeast_gene_annotation{1}, model.genes{i}));
	 			model.geneNames{i} = yeast_gene_annotation{2}{geneIndex(i)};
    else
        geneIndex(i) = nan;
				model.geneNames{i} = ' ';
    end
end

% add protein name for genes
ss4 = length(model.genes);
proteinName = cell(ss4,1);
for i = 1:ss4
proteinName{i} = strcat('COBRARProtein',num2str(i));
model.proteins{i} = proteinName{i};
end

model = rmfield(model,'grRules');

cd ..
saveYeastModel(model)
cd modelCuration
