%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeBiomass(model,P,GAM,NGAM)
%
% Benjamín J. Sánchez. Last edited: 2017-12-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeBiomass(model,P,GAM,NGAM)

%Recalculate GAM based on original GAM:
Pbase   = 0.4668;         %Value from original biomass comp. (Förster et al. 2003)
Cbase   = 0.4514;         %Value from original biomass comp. (Förster et al. 2003)
Pfactor = P/Pbase;
Cfactor = (Cbase+Pbase-P)/Cbase;    %Assumption: change in protein is balanced with a change in carbohydrate
fullGAM = GAM + 16.965*Pfactor + 5.210*Cfactor;

%Get current contents and calculate conversion factors for proteins and carbs:
[Pbase,Cbase] = calculateContent(model);
Pfactor = P/Pbase;
Cfactor = (Cbase+Pbase-P)/Cbase;    %Assumption: change in protein is balanced with a change in carbohydrate

%Change biomass composition:
protPos = strcmp(model.metNames,'protein [cytoplasm]');
carbPos = strcmp(model.metNames,'carbohydrate [cytoplasm]');
bioRxn  = model.S(protPos,:) == -1;
protRxn = model.S(protPos,:) == 1;
carbRxn = model.S(carbPos,:) == 1;
for i = 1:length(model.mets)
    Sbio  = model.S(i,bioRxn);
    Sprot = model.S(i,protRxn);
    Scarb = model.S(i,carbRxn);
    if Sbio ~= 0      
        name  = model.metNames{i};
        isATP = strcmpi(name,'ATP [cytoplasm]');
        isADP = strcmpi(name,'ADP [cytoplasm]');
        isH2O = strcmpi(name,'H2O [cytoplasm]');
        isH   = strcmpi(name,'H+ [cytoplasm]');
        isP   = strcmpi(name,'phosphate [cytoplasm]');
                
        %Variable ATP growth related maintenance (GAM):
        if isATP || isADP || isH2O || isH || isP
            model.S(i,bioRxn) = sign(Sbio)*fullGAM;
        end
        
    elseif Sprot ~= 0 && Sprot ~= 1
        %Variable aa content in biomass eq:
        model.S(i,protRxn) = Sprot*Pfactor;
        
    elseif Scarb ~= 0 && Scarb ~= 1
        %Variable carb content in biomass eq:
        model.S(i,carbRxn) = Scarb*Cfactor;
    end
end

%Add/change NGAM reaction:
%            ATP    +    H2O    ->  ADP     +   H+      +  PO4
mets  = {'s_0434[c]','s_0803[c]','s_0394[c]','s_0794[c]','s_1322[c]'};
coefs = [-1,-1,1,1,1];
model = addReaction(model,{'r_4046','non-growth associated maintenance reaction'}, ...
                    mets,coefs,false,NGAM,NGAM);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%