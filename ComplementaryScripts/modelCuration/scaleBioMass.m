%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scaleBioMass
% Corrects the stoichiometry coefficients of all pseudo-rxns in an
% iterative fashion:
% 1. Switch back to original model's abundance values (Forster et al. 2003)
% 
% Function adapted from SLIMEr: https://github.com/SysBioChalmers/SLIMEr
%
% Benjamin Sanchez. Last update: 2018-09-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scaleBioMass

%Load model:
initCobraToolbox
cd ..
model = loadYeastModel;
cd modelCuration

%Load data from Forster 2003:
fid = fopen('../../ComplementaryData/physiology/biomassComposition_Forster2003.tsv');
Forster2003     = textscan(fid,'%s %s %f32 %f32 %s','Delimiter','\t','HeaderLines',1);
data.mets       = Forster2003{1};
data.abundances = double(Forster2003{3});
data.MWs        = double(Forster2003{4});
data.groups     = Forster2003{5};
fclose(fid);

%Compute current composition:
sumBioMass(model,data);

%Switch to Forster 2003 data:
for i = 1:length(data.mets)
    %Find positions:
    rxnName = [data.groups{i} ' pseudoreaction'];
    rxnName = strrep(rxnName,'other','biomass');
    rxnPos  = strcmp(model.rxnNames,rxnName);
    metPos  = strcmp(model.mets,data.mets{i});
    %Extra changes in protein pseudoreaction:
    if strcmp(data.groups{i},'protein')
        oldVal = abs(model.S(metPos,rxnPos));
        metPos = abs(model.S(:,rxnPos)) == oldVal;
        if sum(metPos) ~= 2
            error('did not found tRNA(aa) pair')
        end
    end
    %Change stoichiometry:
    model.S(metPos,rxnPos) = sign(model.S(metPos,rxnPos))*data.abundances(i);
end
sumBioMass(model,data);

%Finally, save model:
cd ..
saveYeastModel(model)
cd modelCuration

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
