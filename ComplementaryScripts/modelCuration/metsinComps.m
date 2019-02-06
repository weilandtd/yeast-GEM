function [metsindiffcomps] = metsinComps(model,modle_r,mets)
% this function is to identify same metabolites occuring in different
% compartements
% usage: [metsindiffcomps] = SpecificModel(model,mets)
%
% model     model structure to save (note: must be in COBRA format)
% mets      should be a cell array of metlist (here we are matching metname)
%           (opt, default all mets in the model)
%
% Feiran Li 2018.10.12


if nargin<2
    mets = model.metNames;
end
%model_r = ravenCobraWrapper(model);
metsindiffcomps = [];
for i = 1:length(mets)
    S = regexp(mets{i}, ' [', 'split');
    met_temp = char(S(1));
    metMapping = strcmpi(met_temp,modle_r.metNames);
    metsindiffcomps = [metsindiffcomps;model.metNames(metMapping)];
end

