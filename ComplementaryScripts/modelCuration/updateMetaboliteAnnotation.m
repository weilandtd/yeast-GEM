%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update the metabolite annoation information in the model
% updateMetaboliteAnnotation.m is a function from cobra
% March 11, 2018 by Hongzhong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = readCbModel('yeastGEM.xml');
update_annotation = load('updated metabolite annotation.mat')
% update_annotation
% column 1: Abbreviation
% column 2: Description
% column 3: Charged formula
% column 4: Charge
% column 5: Compartment
% column 6: KEGG ID
% column 7: ChEBI ID
function model = updateMetaboliteAnnotation(model,update_annotation)
for  i = 1 : size(model.mets, 1)
model.metNames{i} = update_annotation.metabolite{i,2};
model.metCharges(i,1) = str2double(update_annotation.metabolite{i,4});
model.metKEGGID{i} = update_annotation.metabolite{i,6};
model.metChEBIID{i} = update_annotation.metabolite{i,7};
end

%model.metCharges(isnan(model.metCharges))= []

end



