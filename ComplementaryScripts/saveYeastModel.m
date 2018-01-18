%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveYeastModel(model,version)
% Saves model as a .mat, .xml and .txt file.
%
% Benjamín J. Sánchez. Last edited: 2018-01-18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveYeastModel(model,version)

model.description = ['yeastGEM_' version];

%Remove any space in rxnECNumbers:
model.rxnECNumbers = strrep(model.rxnECNumbers,' ','');

%Save changes to current model:
cd ../ModelFiles/mat
save('yeastGEM.mat','model');
cd ../xml
writeCbModel(model,'sbml','yeastGEM.xml');
cd ../txt
writeCbModel(model,'text','yeastGEM.txt');

%Detect boundary metabolites and save them in a .txt file:
cd ../../ComplementaryScripts
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