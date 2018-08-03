function [L_fixed]=triSurfLogicSharpFix(F,L,dirOpt)

% function [L_fixed]=triSurfLogicSharpFix(F,L,dirOpt)
%-------------------------------------------------------------------------
%
% 
%-------------------------------------------------------------------------

%%

if size(F,1)~=size(L,1)
    error('size(F,1)~=size(L,1)');
end

switch dirOpt
    case 1 %Add "inward teeth" to the list
        indInLogic=unique(F(L,:));
        L_fixed=all(ismember(F,indInLogic),2);
    case 2 %Remove sharp stand-alone triangles from list
        indNotInLogic=unique(F(~L,:));
        L_fixed=L & ~all(ismember(F,indNotInLogic),2);                
    case 3
        [L]=triSurfLogicSharpFix(F,L,1);
        [L]=triSurfLogicSharpFix(F,L,2);        
        L_fixed=L;
end