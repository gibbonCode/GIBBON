function [CV]=faceToVertexMeasure(F,V,CF)
% function [CV]=faceToVertexMeasure(F,V,CF)
% ------------------------------------------------------------------------
%
%
% Change log:
% 2009
% 2019/04/22 Minor changes to improve efficiency through improved allocation
% 2019/04/22 Renamed variables
% 2019/04/23 Updated to handle cell input
% ------------------------------------------------------------------------

%%

%Check if F is a cell of face arrays
cellMode=isa(F,'cell');

if cellMode==1
    CV=repmat({zeros(size(V,1),size(CF{1},2))},1,numel(CF));
    for qc=1:1:numel(F)
        % Get vertex-face connectivity matrix
        conStruct=patchConnectivity(F{qc},V,'vf');
        IND_F=conStruct.vertex.face;
        logicValid=IND_F>0;
        CV_now=ones(size(V,1),size(IND_F,2),size(CF{qc},2));
        cv=nan(size(IND_F));
        for q=1:1:size(CF{qc},2) %Loop over data columns
            cf=CF{qc}(:,q);
            cv(logicValid)=cf(IND_F(logicValid));
            CV_now(:,:,q)=cv;
        end
        CV{qc}=CV_now;
    end
    CV=mean(cell2mat(CV),2,'omitnan'); %Take mean across faces
    if size(CF{1},2)>1
        CV=permute(CV,[1 3 2]); %Get rid of singleton dimension
    end
else
    % Get vertex-face connectivity matrix
    conStruct=patchConnectivity(F,V,'vf');
    IND_F=conStruct.vertex.face;    
    logicValid=IND_F>0;    
    CV=ones(size(V,1),size(CF,2));
    cv=nan(size(IND_F));
    for q=1:1:size(CF,2) %Loop over data columns
        cf=CF(:,q);
        cv(logicValid)=cf(IND_F(logicValid));
        CV(:,q)=mean(cv,2,'omitnan');
    end
end

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
