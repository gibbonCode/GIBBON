function [varargout]=hexVol(E,V)

% function VE=hexVol(E,V)
% ------------------------------------------------------------------------
% Calculates the volume (VE) of the hexahedral elements specified by the
% element matrix E and the vertices V.
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2018/04/10
%-------------------------------------------------------------------------

%%

[Et,Vt]=hex2tet(E,V,[],5);
[VE,logicPositive]=tetVol(Et,Vt);

if size(E,1)==1
    VE=sum(VE);
    logicPositive=all(logicPositive);
else   
    VE=sum(reshape(VE,[5 size(E,1)])',2);    
    logicPositive=all(reshape(logicPositive,[5 size(E,1)])',2);
end

%%

varargout{1}=abs(VE);
varargout{2}=logicPositive;
