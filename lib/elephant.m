function [F,V]=elephant

% [F,V]=elephant
% ------------------------------------------------------------------------
%
% This function generates patch data (faces=F and vertices=V) for an
% elephant model.
%
% This adjusted model consists of 892 triangular faces and 448 vertices. 
% 
% The model was constructed based on the model given here: 
% https://www.rocq.inria.fr/gamma/gamma/download/affichage.php?dir=DINOSAUR/&name=Parasaurolophus
% 
%
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
%------------------------------------------------------------------------
%%

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
offPath=fullfile(toolboxPath,'data','OFF');
fileName=fullfile(offPath,'elephant-50kv.off');
[F,V] = import_off(fileName);



