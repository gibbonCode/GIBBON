function [P]=tesSmooth(TES,V,IND_V,cPar)

% function [P]=tesSmooth(TES,V,IND_V,cPar)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/06/02
%------------------------------------------------------------------------

%% 

%Get/set method
if isfield(cPar,'Method')
    smoothMethod=cPar.Method;
else
    smoothMethod='LAP'; %DEFAULT
end

%Smooth
switch smoothMethod
    case 'LAP'
        [P]=tesSmooth_LAP(TES,V,IND_V,cPar);
    case 'HC'
        [P]=tesSmooth_HC(TES,V,IND_V,cPar);
    otherwise
        error('Invalid smooth method specified');
end