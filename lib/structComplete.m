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

C=A; %Initialize C as the same as A
for q=1:1:numel(fieldNameSet) %Loop over field names
    fieldNameNow=fieldNameSet{q}; %Current field name
    if isfield(A,fieldNameNow) %If A contains the field in B                
        
        if isempty(A.(fieldNameNow)) %if the field in A is empty
            if emptyFixOpt==1 %if empty field fixing is on
                C.(fieldNameNow)=B.(fieldNameNow); %Replace empty in C by default
            end
        end    
        
        if isstruct(A.(fieldNameNow)) && isstruct(B.(fieldNameNow)) %If the field in A is a structure, check structure recursively
            [C.(fieldNameNow)]=structComplete(A.(fieldNameNow),B.(fieldNameNow),emptyFixOpt);
        end
        
    else %If the field is missing, add it
        C.(fieldNameNow)=B.(fieldNameNow);            
    end
end

%%