%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveYeastModel(model)
% Saves model as a .mat, .xml and .txt file.
%
% Benjamín J. Sánchez. Last edited: 2018-01-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveYeastModel(model)

model.description = 'yeastGEM';

%Remove any space in rxnECNumbers:
model.rxnECNumbers = strrep(model.rxnECNumbers,' ','');

%Save changes to current model:
save('../ModelFiles/mat/yeastGEM.mat','model');
writeCbModel(model,'sbml','../ModelFiles/xml/yeastGEM.xml');
writeCbModel(model,'text','../ModelFiles/txt/yeastGEM.txt');

%Detect boundary metabolites and save them in a .txt file:
fid = fopen('boundaryMets.txt','wt');
for i = 1:length(model.rxns)
    pos = find(model.S(:,i) ~= 0);
    if length(pos) == 1 %Exchange rxn
        fprintf(fid,[model.mets{pos} '\t' model.metNames{pos} '\n']);
    end
end
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%