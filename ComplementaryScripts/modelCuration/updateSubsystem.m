%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update the metabolite annoation information in the model
% updateMetaboliteAnnotation.m is a function from cobra
% March 11, 2018 by Hongzhong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = readCbModel('yeastGEM.xml');
update_annotation = load('subsystem.mat')
% update_annotation
% column 1: reactionID
% column 2: corrected subsystem

function model = updateSubsystem(model,update_annotation)
for  i = 1 : size(model.subSystems, 1)
model.subSystems{i} = update_annotation.subsystem{i,2};
end


end



