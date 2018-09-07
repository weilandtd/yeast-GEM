%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,P,C,R,D,L,N] = sumBioMass(model)
% Calculates breakdown of biomass
%
% model     metabolic model in COBRA format
% data      structure with at least the following 2 fields:
%   mets    Cell array with metabolite ids
%   MWs     Numeric array with molecular weights for each metabolite
%
% X         Total biomass fraction [gDW/gDW]
% P         Protein fraction [g/gDW]
% C         Carbohydrate fraction [g/gDW]
% R         RNA fraction [g/gDW]
% D         DNA fraction [g/gDW]
% L         Lipid fraction [g/gDW]
% N         Other mets [g/gDW]
%
% Function adapted from SLIMEr: https://github.com/SysBioChalmers/SLIMEr
%
% Benjamin Sanchez. Last update: 2018-09-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,P,C,R,D,L,N] = sumBioMass(model,data)

%Get main fractions:
[P,X] = getFraction(model,data,'P',0);
[C,X] = getFraction(model,data,'C',X);
[R,X] = getFraction(model,data,'R',X);
[D,X] = getFraction(model,data,'D',X);
[L,X] = getFraction(model,data,'L',X);
[N,X] = getFraction(model,data,'N',X);

disp(['X -> ' num2str(X) ' gDW/gDW'])

% Simulate growth:
sol = optimizeCbModel(model);
disp(['Growth = ' num2str(sol.f) ' 1/h'])
disp(' ')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,X] = getFraction(model,data,compType,X)

%Define pseudoreaction name:
rxnName = [compType ' pseudoreaction'];
rxnName = strrep(rxnName,'P','protein');
rxnName = strrep(rxnName,'C','carbohydrate');
rxnName = strrep(rxnName,'N','biomass');
rxnName = strrep(rxnName,'L','lipid backbone');
rxnName = strrep(rxnName,'R','RNA');
rxnName = strrep(rxnName,'D','DNA');

%Add up fraction:
rxnPos = strcmp(model.rxnNames,rxnName);
isSub   = model.S(:,rxnPos) < 0;        %substrates in pseudo-rxn
if strcmp(compType,'L')
    F = -sum(model.S(isSub,rxnPos));   %g/gDW
else
    F = 0;
    %Add up all components:
    for i = 1:length(model.mets)
        pos = strcmp(data.mets,model.mets{i});
        if isSub(i) && sum(pos) == 1
            if strcmp(compType,'N')
                MW = data.MWs(pos);
            else
                MW = data.MWs(pos)-18;
            end
            abundance = -model.S(i,rxnPos)*MW/1000;
            F         = F + abundance;
        end
    end
end
X = X + F;

disp([compType ' -> ' num2str(F) ' g/gDW'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
