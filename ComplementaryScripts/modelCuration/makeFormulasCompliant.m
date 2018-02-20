%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeFormulasCompliant.m
% Corrects SBML-incompatible metabolite formulas.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%aa codes:
aas = {'s_0955[c]'	'ala'   	% A     Alanine
       's_0981[c]'	'cys'       % C     Cysteine
       's_0973[c]'	'asp'       % D     Aspartic acid
       's_0991[c]'	'glu'       % E     Glutamic acid
       's_1032[c]'	'phe'   	% F     Phenylalanine
       's_1003[c]'	'gly'       % G     Glycine
       's_1006[c]'	'his'   	% H     Histidine
       's_1016[c]'	'ile'       % I     Isoleucine
       's_1025[c]'	'lys'       % K     Lysine
       's_1021[c]'	'leu'   	% L     Leucine
       's_1029[c]'	'met'       % M     Methionine
       's_0969[c]'	'asn'   	% N     Asparagine
       's_1035[c]'	'pro'       % P     Proline
       's_0999[c]'	'gln'   	% Q     Glutamine
       's_0965[c]'	'arg'       % R     Arginine
       's_1039[c]'	'ser'   	% S     Serine
       's_1045[c]'	'thr'       % T     Threonine
       's_1056[c]'	'val'       % V     Valine
       's_1048[c]'	'trp'       % W     Tryptophan
       's_1051[c]'	'tyr'};     % Y     Tyrosine

%tRNA's are just produced and consumed in single loops, so we can just remove
%the radical part and still maintain mass balances, i.e. aa-tRNA(aa) will have
%the same formula as aa and tRNA(aa) won't have any formula at all. Example:
%     tRNA(Ala): C10H17 O10PR2(C5H8O6PR)n
% Ala-tRNA(Ala): C13H22NO11PR2(C5H8O6PR)n
%    diff (ala): C3 H5 NO
% Cycle in which the 2 are involved:
% r_0157: ATP + L-alanine + tRNA(Ala) -> Ala-tRNA(Ala) + AMP + diphosphate
% r_4047: 0.4193 Ala-tRNA(Ala) + ...  -> 0.4193 tRNA(Ala) + ... + protein
for i = 1:length(model.mets)
    name    = model.metNames{i};
    formula = model.metFormulas{i};
    if contains(model.metNames{i},'-tRNA(')
        %Correct metabolite:
        try
            aaName = lower(name(1:strfind(name,'-')-1));
            aaID   = aas{strcmp(aas(:,2),aaName),1};
            aaPos  = strcmp(model.mets,aaID);
            model.metFormulas{i} = model.metFormulas{aaPos};
        catch
            model.metFormulas{i} = 'NA';    %f-met CHANGE LATER
        end
        
        %Correct associated tRNA pair:
        pairName = name(strfind(name,'-')+1:end);
        pairPos  = strcmp(model.metNames,pairName);
        if sum(pairPos) > 0
            model.metFormulas{pairPos} = 'NA';
        end
    end
end

%Display the remaining problems:
for i = 1:length(model.mets)
    if contains(model.metFormulas{i},')n')
        disp([model.metNames{i} ': ' model.metFormulas{i}]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%