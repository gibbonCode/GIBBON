function [b]=mcol(a)

% function [b]=mcol(a)
% ------------------------------------------------------------------------
% This function converts the input matrix (or vector) a to a column vector. 
%
%
% Change log: 
% 2023/03/10: Kevin Moerman Created
%
% ------------------------------------------------------------------------
%%

if ismatrix(a)
    if ~isvector(a) %Check if it is a vector already
        b=a(:);
    else %Vector so either row already or a column
        if ~iscolumn(a)
            b=a'; %Transpose
        else
            b=a; %Keep
        end
    end
else
    error('Input is not of matrix type. The mrow function is only defined for matrix arrays, see also ismatrix');
end