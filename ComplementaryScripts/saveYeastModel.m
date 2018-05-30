%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveYeastModel(model)
% Saves model as a .xml, .txt and .yml file. Also updates complementary
% files (boundaryMets.txt and dependencies.txt).
%
% Benjamín J. Sánchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveYeastModel(model)

%Remove any space in rxnECNumbers:
model.rxnECNumbers = strrep(model.rxnECNumbers,' ','');

%Save changes to current model:
writeCbModel(model,'sbml','../ModelFiles/xml/yeastGEM.xml');
writeCbModel(model,'text','../ModelFiles/txt/yeastGEM.txt');
exportForGit(model,'yeastGEM','..',{'yml'});

%Detect boundary metabolites and save them in a .txt file:
fid = fopen('../ModelFiles/boundaryMets.txt','wt');
for i = 1:length(model.rxns)
    pos = find(model.S(:,i) ~= 0);
    if length(pos) == 1 %Exchange rxn
        fprintf(fid,[model.mets{pos} '\t' model.metNames{pos} '\n']);
    end
end
fclose(fid);

%Convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
%inconsistencies between Windows and MAC:
copyfile('../ModelFiles/xml/yeastGEM.xml','backup.xml')
fin  = fopen('backup.xml','r');
fout = fopen('../ModelFiles/xml/yeastGEM.xml','w');
still_reading = true;
while still_reading
    inline = fgets(fin);
    if ~ischar(inline)
        still_reading = false;
    else
        if length(regexp(inline,'[0-9]e-00[0-9]')) == 1 && length(regexp(inline,'e-00')) == 1
            inline = strrep(inline,'e-00','e-0');
        end
        fwrite(fout,inline);
    end
end
fclose('all');
delete('backup.xml');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
