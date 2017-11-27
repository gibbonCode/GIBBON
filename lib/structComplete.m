function [C]=structComplete(A,B,emptyFixOpt)

% function [C]=structComplete(A,B,emptyFixOpt)
% ------------------------------------------------------------------------
% This function fills in the missing data in the structure A with the
% content from B. The structure B can be seen as a default input structure
% for instance and A might be an incomplete set of inputs. The optional
% parameter emptyFixOpt (0=no, 1=yes) determines whether empty entries in A
% are overwritten by the corresponding default values in B. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2017/11/20
%------------------------------------------------------------------------

%%
fieldNameSet=fieldnames(B);

C=B;
for q=1:1:numel(fieldNameSet)
    fieldNameNow=fieldNameSet{q};
    if isfield(A,fieldNameNow)
        C.(fieldNameNow)=A.(fieldNameNow);
        if emptyFixOpt==1 %if empty field fixing is on
            if isempty(C.(fieldNameNow)) %if the current field is empty
                C.(fieldNameNow)=B.(fieldNameNow); %Replace empty by default
            end
        end        
    end
end