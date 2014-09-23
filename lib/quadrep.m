function [Mq]=quadrep(M)

% function [Mq]=quadrep(M)
% ------------------------------------------------------------------------
% The input matrix 'M' is mirrored in three directions to create the matrix
% 'Mq'. The size of this matrix is ((2*size(M))-1).
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 15/09/2008
% ------------------------------------------------------------------------

%% Setting up Mq matrix
Mq=zeros((2*size(M))-1);

%% LEFT SIDE FRONT

%Front top left
Mq(1:size(M,1),1:size(M,2),1:size(M,3))=M;

%Front bottom left
Mq(size(M,1):end,1:size(M,2),1:size(M,3))=flipdim(M,1);

%% RIGHT SIDE FRONT
Mq(1:end,size(M,2):end,1:size(M,3))=flipdim(Mq(1:end,1:size(M,2),1:size(M,3)),2);

%% BACK
Mq(1:end,1:end,size(M,3):end)=flipdim(Mq(1:end,1:end,1:size(M,3)),3);

%% END