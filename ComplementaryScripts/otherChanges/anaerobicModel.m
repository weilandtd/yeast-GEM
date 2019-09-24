%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anaerobicModel.m
% Converts model to anaerobic
%
% Benjamín J. Sánchez
% Feiran Li -Last update: 2019-09-24 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1st change: Changes media to anaerobic (no O2 uptake and allows sterol
%            and fatty acid exchanges)
model.lb(strcmp(model.rxns,'r_1992')) = 0;        %O2
model.lb(strcmp(model.rxns,'r_1757')) = -1000;    %ergosterol
model.lb(strcmp(model.rxns,'r_1915')) = -1000;    %lanosterol
model.lb(strcmp(model.rxns,'r_1994')) = -1000;    %palmitoleate
model.lb(strcmp(model.rxns,'r_2106')) = -1000;    %zymosterol
model.lb(strcmp(model.rxns,'r_2134')) = -1000;    %14-demethyllanosterol
model.lb(strcmp(model.rxns,'r_2137')) = -1000;    %ergosta-5,7,22,24(28)-tetraen-3beta-ol
model.lb(strcmp(model.rxns,'r_2189')) = -1000;    %oleate

%2nd change: Removes the requirement of heme a in the biomass equation
%            (not used under aerobic conditions)
mets = {'s_3714[c]','s_1198[c]','s_1203[c]','s_1207[c]','s_1212[c]'};
[~,met_index] = ismember(mets,model.mets);
model.S(met_index,strcmp(model.rxns,'r_4598')) = 0;

%3rd change: Blocked pathways for proper glycerol production
%Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
model.lb(strcmp(model.rxns,'r_0713')) = 0; %Mithocondria
model.lb(strcmp(model.rxns,'r_0714')) = 0; %Cytoplasm
%Block glycerol dehydroginase (only acts in microaerobic conditions)
model.ub(strcmp(model.rxns,'r_0487')) = 0;
%Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
model.ub(strcmp(model.rxns,'r_0472')) = 0;

%4th change: Refit GAM and NGAM to exp. data, change biomass composition
GAM   = 30.49;  %Data from Nissen et al. 1997
P     = 0.461;  %Data from Nissen et al. 1997
NGAM  = 0;      %Refit done in Jouthen et al. 2012
%cd ../modelCuration
model = changeBiomass(model,P,GAM,NGAM);
%cd ../otherChanges
ScaleBiomassAnaerobic(model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ScaleBiomassAnaerobic(model)
fid = fopen('../../ComplementaryData/physiology/biomassComposition_Forster2003.tsv');
Forster2003     = textscan(fid,'%s %s %f32 %f32 %s','Delimiter','\t','HeaderLines',1);
data.mets       = Forster2003{1};
data.abundances = double(Forster2003{3});
data.MWs        = double(Forster2003{4});
data.groups     = Forster2003{5};
fclose(fid);
%Correct with data from biomassComposition_Cofactor_Ion:
data_original = data;
fid = fopen('../../ComplementaryData/physiology/biomassComposition_Cofactor_Ion.tsv');
Cofactors       = textscan(fid,'%s %s %f32 %f32 %s %s','Delimiter','\t','HeaderLines',1);
data.mets       = Cofactors{1};
data.abundances = double(Cofactors{3});
data.MWs        = double(Cofactors{4});
data.groups     = Cofactors{5};
fclose(fid);
for j = 1:length(data_original.mets)
    if ~ismember(data_original.mets(j),data.mets)
        data.mets = [data.mets; data_original.mets(j)];
        data.abundances = [data.abundances; data_original.abundances(j)];
        data.MWs = [data.MWs; data_original.MWs(j)];
        data.groups = [data.groups; data_original.groups(j)];
    end
end

%Balance out mass with carbohydrate content:
[X,P,C,R,~,~,~,~] = sumBioMass(model,data);
delta = X - 1;           %difference to balance
fC    = (C - delta)/C;
model = rescalePseudoReaction(model,'carbohydrate',fC);
sumBioMass(model,data);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = rescalePseudoReaction(model,metName,f)

rxnName = [metName ' pseudoreaction'];
rxnPos  = strcmp(model.rxnNames,rxnName);
for i = 1:length(model.mets)
    S_ir   = model.S(i,rxnPos);
    isProd = strcmp(model.metNames{i},[metName ' [cytoplasm]']);
    if S_ir ~= 0 && ~isProd
        model.S(i,rxnPos) = f*S_ir;
    end
end

end