%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% increaseVersion(model,version)
%
% Benjamín J. Sánchez. Last edited: 2018-01-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function increaseVersion(model,version)

%Include tag in model:
model.description = ['yeastGEM_v' version];
saveYeastModel(model)

%Update version file:
fid = fopen('../version.txt','wt');
fprintf(fid,version);
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%