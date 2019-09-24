function model = NitrogenlimitedBiomass(model)

% This function is to change biomass composition to nitrogen limited
% condition. As we know, the protein content drop drmatically under nitogen
% limited condition, this code is to rescale the protein content.

% Feiran Li 2018-08-25

%Compute current composition:
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
cd ../modelCuration/
[X,P,C,R,~,L,~,~] = sumBioMass(model,data);


fC    = 0.587/C;
fR    = 0.077/R;
fP    = 0.289/P;
fL    = 0.048/L;
model = rescalePseudoReaction(model,'carbohydrate',fC);
model = rescalePseudoReaction(model,'lipid backbone',fL); 
model = rescalePseudoReaction(model,'lipid chain',fL); 
model = rescalePseudoReaction(model,'protein',fP);
model = rescalePseudoReaction(model,'RNA',fR);

[X,P,C,R,~,L,~,~] = sumBioMass(model,data);


%Balance out mass with carbohydrate content:
delta = X - 1;           %difference to balance
fC    = (C - delta)/C;
model = rescalePseudoReaction(model,'carbohydrate',fC);
sumBioMass(model,data);
cd ../otherChanges/
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
end
