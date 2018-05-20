%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addSLIMErxns
%
% Benjamín J. Sánchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Start COBRA & load model:
clear variables
initCobraToolbox
cd ..
model = loadYeastModel;

%Move to SLIMEr and add paths:
cd ../../SLIMEr/models
addpath('../simulations')
addpath('../data')

%Add SLIMErxns:
data = readLahtveeData(1);
[~,model_SLIMEr,~,~] = modelsFromData(model,data,'backbones');

%Remove added paths:
rmpath('../simulations')
rmpath('../data')

%Reinsert previous GAM and correct O.F.:
GAM    = 61.9779;
bioRxn = strcmp(model_SLIMEr.rxnNames,'biomass pseudoreaction');
ATPpos = strcmp(model_SLIMEr.metNames,'ATP [cytoplasm]');
H2Opos = strcmp(model_SLIMEr.metNames,'H2O [cytoplasm]');
ADPpos = strcmp(model_SLIMEr.metNames,'ADP [cytoplasm]');
Hpos   = strcmp(model_SLIMEr.metNames,'H+ [cytoplasm]');
Ppos   = strcmp(model_SLIMEr.metNames,'phosphate [cytoplasm]');
model_SLIMEr.S(ATPpos,bioRxn) = -GAM;
model_SLIMEr.S(H2Opos,bioRxn) = -GAM;
model_SLIMEr.S(ADPpos,bioRxn) = +GAM;
model_SLIMEr.S(Hpos,bioRxn)   = +GAM;
model_SLIMEr.S(Ppos,bioRxn)   = +GAM;

%Save model & return to origin:
cd ../../yeast-GEM/complementaryScripts
saveYeastModel(model_SLIMEr);
cd otherChanges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%