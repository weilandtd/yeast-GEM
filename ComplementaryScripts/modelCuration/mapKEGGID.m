% This script maps reaction IDs of different databases
% to find new rxnKEGGID(s)
%
% Inputs: model
%
% getRxnsFromMetaCyc is a function from RAVEN
%
% reactions.tsv from modelSEED database is accessed directly from cloned GitHub repository ModelSEEDDatabase
% link: github.com/ModelSEED/ModelSEEDDatabase
%
% Cheng Wei Quan (Eiden), 2020-05-07

%load model
cd ..
model = loadYeastModel;

%% rxnKEGGID
%map rxnMetaNetXID to rxnKEGGID via function mapIDsViaMNXref
cd ./missingFields
xref_rxnKEGG = mapIDsViaMNXref('rxns',model.rxnMetaNetXID,'MetaNetX','KEGG');
xref_rxnKEGG(:,2) = model.rxns; %match rxnKEGGID with model.rxns
xref_rxnKEGG(:,3) = model.rxnMetaNetXID; %match rxnKEGGID with model.rxnMetaNetXID
xref_rxnKEGG(:,4) = model.rxnKEGGID; %match rxnKEGGID with model.rxnKEGGID

%only retain mapped rxnKEGGID if model.rxnKEGGID is empty
empties = find(cellfun('isempty',xref_rxnKEGG(:,1)));
xref_rxnKEGG(empties,:) = [];
empties = find(~cellfun('isempty',xref_rxnKEGG(:,4)));
xref_rxnKEGG(empties,:) = [];

%add rxnKEGGID into model
[~,idx] = ismember(xref_rxnKEGG(:,2),model.rxns);
idx = idx(idx~=0);

for i = 1:length(idx)
    model.rxnKEGGID(idx(i)) = xref_rxnKEGG(i,1);
end

%correction for MNXref error in MNXR96123/R00124#2
idx = find(ismember(model.rxnMetaNetXID,'MNXR96123'));
model.rxnKEGGID(idx) = {'R00124'};

%map rxnMetaNetXID to rxnMetaCycID via function mapIDsViaMNXref, then map rxnMetaCycID to rxnKEGGID
xref_rxnMetaCyc = mapIDsViaMNXref('rxns',model.rxnMetaNetXID,'MetaNetX','MetaCyc');
xref_rxnMetaCyc(:,2) = model.rxns; %match rxnMetaCycID with model.rxns
xref_rxnMetaCyc(:,3) = model.rxnMetaNetXID; %match rxnMetaCycID with model.rxnMetaNetXID
xref_rxnMetaCyc(:,4) = model.rxnKEGGID; %match rxnMetaCycID with model.rxnKEGGID
empties = find(cellfun('isempty',xref_rxnMetaCyc(:,1)));
xref_rxnMetaCyc(empties,:) = [];

%load metaCycRxns.mat via getRxnsFromMetaCyc from RAVEN directory
MetaCyc_rxnInfo = getRxnsFromMetaCyc('..RAVEN/external/metacyc');

%Input rxnKEGGID into xref array
for i = 1:length(xref_rxnMetaCyc)
    if isempty(xref_rxnMetaCyc{i,4})
        idx = find(ismember(MetaCyc_rxnInfo.rxns,xref_rxnMetaCyc(i,1)));
        if idx~=0
            if ~isempty(MetaCyc_rxnInfo.rxnMiriams{idx})
                temp = cell2mat(MetaCyc_rxnInfo.rxnMiriams(idx));
                idx2 = find(ismember(temp.name,'kegg.reaction'));
                if ~isempty(idx2)
                    xref_rxnMetaCyc(i,5) = temp.value(idx2);
                end
            end
        end
    end
end

if size(xref_rxnMetaCyc,2) == 5
    %only retain mapped rxnKEGGID if model.rxnKEGGID is empty
    empties = find(cellfun('isempty',xref_rxnMetaCyc(:,5)));
    xref_rxnMetaCyc(empties,:) = [];
    empties = find(~cellfun('isempty',xref_rxnMetaCyc(:,4)));
    xref_rxnMetaCyc(empties,:) = [];
    
    %add rxnKEGGID into model
    [~,idx] = ismember(xref_rxnMetaCyc(:,2),model.rxns);
    idx = idx(idx~=0);
    
    for i = 1:length(idx)
        model.rxnKEGGID(idx(i)) = xref_rxnMetaCyc(i,5);
    end
else
    xref_rxnMetaCyc(:) = []; %erase all cell array, since no KEGGID matched and added
end

%map rxnMetaNetXID to rxnSEEDID via function mapIDsViaMNXref, then map rxnSEEDID to rxnKEGGID
xref_rxnSEED = mapIDsViaMNXref('rxns',model.rxnMetaNetXID,'MetaNetX','SEED');
xref_rxnSEED(:,2) = model.rxns; %match rxnSEEDID with model.rxns
xref_rxnSEED(:,3) = model.rxnMetaNetXID; %match rxnSEEDID with model.rxnMetaNetXID
xref_rxnSEED(:,4) = model.rxnKEGGID; %match rxnSEEDID with model.rxnKEGGID
empties = find(cellfun('isempty',xref_rxnSEED(:,1)));
xref_rxnSEED(empties,:) = [];

%load reactions.tsv from modelSEEDDatabase
temp = what('ModelSEEDDatabase\Biochemistry'); %require cloned/downloaded GitHub repo ModelSEEDDatabase to be added to MATLAB path
cd(temp.path);
fid2 = fopen('reactions.tsv');
format = repmat('%s ',1,22);
format = strtrim(format);
rxn_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(rxn_temp)
    SEED_rxnInfo(:,i) = rxn_temp{i};
end
commentLines = startsWith(SEED_rxnInfo(:,1),'#');
SEED_rxnInfo(commentLines,:) = [];
fclose(fid2);

%Input rxnKEGGID into xref array
for i = 1:length(xref_rxnSEED)
    if isempty(xref_rxnSEED{i,4})
        idx = find(ismember(SEED_rxnInfo(:,1),xref_rxnSEED(i,1)));
        if idx~=0
            if ~isempty(SEED_rxnInfo{idx,13}) && contains(SEED_rxnInfo(idx,13),'KEGG')
                temp_ID = extractBetween(SEED_rxnInfo(idx,13),'KEGG: ','|');
                xref_rxnSEED(i,5) = temp_ID;
            end
        end
    end
end

if size(xref_rxnSEED,2) == 5
    %only retain mapped rxnKEGGID if model.rxnKEGGID is empty
    empties = find(cellfun('isempty',xref_rxnSEED(:,5)));
    xref_rxnSEED(empties,:) = [];
    empties = find(~cellfun('isempty',xref_rxnSEED(:,4)));
    xref_rxnSEED(empties,:) = [];
    
    %add rxnKEGGID into model
    [~,idx] = ismember(xref_rxnSEED(:,2),model.rxns);
    idx = idx(idx~=0);
    
    for i = 1:length(idx)
        model.rxnKEGGID(idx(i)) = xref_rxnSEED(i,5);
    end
else
    xref_rxnSEED(:) = []; %erase all cell array, since no KEGGID matched and added
end

%Save model
temp = what('yeast-GEM\ComplementaryScripts');
cd(temp.path);
saveYeastModel(model);
