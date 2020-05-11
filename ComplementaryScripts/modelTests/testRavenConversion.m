function testRavenConversion
% testRavenConversion
%   Test that changes to the wrapper don't affect the model structure and
%   which fields get lost applying it twice.
%
%   Benjamin J. Sanchez, 2020-05-07
%

% Load model:
cd ..
initCobraToolbox
model_c1 = loadYeastModel;

% Convert to RAVEN, back to COBRA, and back again to RAVEN:
model_r1 = ravenCobraWrapper(model_c1);
model_c2 = ravenCobraWrapper(model_r1);
model_r2 = ravenCobraWrapper(model_c2);

% Compare models:
compareModels(model_c1,model_c2,'COBRA')
compareModels(model_r1,model_r2,'RAVEN')

% Save model:
saveYeastModel(model_c2);
cd missingFields

end

function compareModels(model1,model2,toolbox)

% Sort Miriam fields:
if strcmp(toolbox,'RAVEN')
    for i = 1:length(model1.rxnMiriams)
        [model1.rxnMiriams{i}.name,order] = sort(model1.rxnMiriams{i}.name);
        model1.rxnMiriams{i}.value = model1.rxnMiriams{i}.value(order);
        [model2.rxnMiriams{i}.name,order] = sort(model2.rxnMiriams{i}.name);
        model2.rxnMiriams{i}.value = model2.rxnMiriams{i}.value(order);
    end
    for i = 1:length(model1.metMiriams)
        [model1.metMiriams{i}.name,order] = sort(model1.metMiriams{i}.name);
        model1.metMiriams{i}.value = model1.metMiriams{i}.value(order);
        [model2.metMiriams{i}.name,order] = sort(model2.metMiriams{i}.name);
        model2.metMiriams{i}.value = model2.metMiriams{i}.value(order);
    end
end

fieldnames1 = fieldnames(model1);
fieldnames2 = fieldnames(model2);
fieldnamesAll = union(fieldnames1,fieldnames2);
for i = 1:length(fieldnamesAll)
    if ~isfield(model1,fieldnamesAll{i})
        disp([fieldnamesAll{i} ' is not in ' toolbox ' model before wrapping'])
    elseif ~isfield(model2,fieldnamesAll{i})
        disp([fieldnamesAll{i} ' is not in ' toolbox ' model after wrapping'])
    else
        if ~eval(['isequal(model1.' fieldnamesAll{i} ',model2.' fieldnamesAll{i} ')'])
            disp([fieldnamesAll{i} ' changed after ' toolbox ' wrapping'])
        end
    end
end

end