%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveYeastModel(model)
% Saves model as a .xml and .txt file. Also updates complementary files
% (boundaryMets.txt and dependencies.txt).
%
% Benjamín J. Sánchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveYeastModel(model)

%Remove any space in rxnECNumbers:
model.rxnECNumbers = strrep(model.rxnECNumbers,' ','');

%Remove any empty notes (created automatically by addReaction):
model.rxnNotes = strrep(model.rxnNotes,'NOTES: NA','');

%Save changes to current model:
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

%Retrieve RAVEN version:
RAVENver = getVersion('checkInstallation.m','version.txt');

%Retrieve latest COBRA commit:
COBRApath   = which('initCobraToolbox.m');
slashPos    = getSlashPos(COBRApath);
COBRApath   = COBRApath(1:slashPos(end)-1);
currentPath = pwd;
cd(COBRApath)
COBRAcommit = git('log -n 1 --format=%H');
cd(currentPath)

%Save file with versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['RAVEN_toolbox\tv' RAVENver '\n']);
fprintf(fid,['COBRA_toolbox\tcommit ' COBRAcommit(1:7) '\n']);
fields = fieldnames(model.modelVersion);
for i = 1:length(fields)
    value = model.modelVersion.(fields{i});
    fprintf(fid,[fields{i} '\t' num2str(value) '\n']);
end
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function version = getVersion(IDfileName,VERfileName)

try
    path     = which(IDfileName);
    slashPos = getSlashPos(path);
    path     = path(1:slashPos(end-1));
    fid      = fopen([path VERfileName],'r');
    version  = fscanf(fid,'%s');
    fclose(fid);
catch
    version = '?';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slashPos = getSlashPos(path)

slashPos = strfind(path,'\');       %Windows
if isempty(slashPos)
    slashPos = strfind(path,'/');   %MAC/Linux
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%