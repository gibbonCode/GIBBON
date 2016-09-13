function M=checkerBoard3D(siz)

% function M=checkerBoard3D(siz)
% ------------------------------------------------------------------------
% This function creates a checkboard image of the size siz whereby elements
% are either black (0) or white (1). The first element is white.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2008/08/15
% ------------------------------------------------------------------------

%%

%Coping with 1D or 2D input
if numel(siz)==2
    siz(3)=1; 
elseif numel(siz)==1
    siz(2:3)=[1 1];
end

%%
[I,J,K]=ndgrid(1:1:siz(1),1:1:siz(2),1:1:siz(3));
logic_ij=((iseven(I)| iseven(J)) & ((iseven(I)~=iseven(J))));
M=false(siz);
M(iseven(K))=logic_ij(iseven(K));
M(~iseven(K))=~logic_ij(~iseven(K));