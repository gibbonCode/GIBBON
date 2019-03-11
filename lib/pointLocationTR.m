function [TI,BC]=pointLocationTR(varargin)

% function [TI,BC]=pointLocationTR(TR,QP,distCropOpt,chullCropOpt,waitbarOpt,toleranceMagnitude)
% ------------------------------------------------------------------------
% 
% This function finds the triangles or tetrahedrons in which the query
% points QP are contained. In addition it outputs the barycentric
% coordinates.
%
%
% [TI,BC]=pointLocationTR(TR,QP,distCropOpt,chullCropOpt,waitbarOpt,toleranceMagnitude)
% [TI,BC]=pointLocationTR(TR,QP,optionStruct)
%
% Default options: 
% optionStructDef.distCropOpt=0;
% optionStructDef.chullCropOpt=0;
% optionStructDef.waitbarOpt=1;
% optionStructDef.toleranceMagnitude=0;
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% Change log:
% 2014/04/01 Created
% 2018/08/03 Fixed bug in relation to tolerance level 
% 2018/08/03 Added tolerance level as input
% 2018/08/03 Added input structure to handle options
%------------------------------------------------------------------------
%% Parse input
TR=varargin{1};
QP=varargin{2};

%Default options
optionStructDef.distCropOpt=0;
optionStructDef.chullCropOpt=0;
optionStructDef.waitbarOpt=1;
optionStructDef.toleranceMagnitude=0;

switch nargin
    case 2
        optionStruct=[];
    case 3
        if isa(varargin{3},'struct')
            optionStruct=varargin{3};
        else
            optionStruct.distCropOpt=varargin{3};
        end
    case 4
        optionStruct.distCropOpt=varargin{3};
        optionStruct.chullCropOpt=varargin{4};
        optionStruct.waitbarOpt=[];
        optionStruct.toleranceMagnitude=[];
    case 5
        optionStruct.distCropOpt=varargin{3};
        optionStruct.chullCropOpt=varargin{4};        
        optionStruct.waitbarOpt=varargin{5};
        optionStruct.toleranceMagnitude=[];
    case 6    
        optionStruct.distCropOpt=varargin{3};
        optionStruct.chullCropOpt=varargin{4};
        optionStruct.waitbarOpt=varargin{5};
        optionStruct.toleranceMagnitude=varargin{6};
end

[optionStruct]=structComplete(optionStruct,optionStructDef,1); %Complement provided with default if missing
distCropOpt=optionStruct.distCropOpt;
chullCropOpt=optionStruct.chullCropOpt;
waitbarOpt=optionStruct.waitbarOpt;
toleranceMagnitude=optionStruct.toleranceMagnitude;

%%
E=TR.ConnectivityList;
numPerE=size(E,2);
numElements=size(E,1);

if distCropOpt==1
    [TR_centre_all,TR_rad_all] = circumcenter(TR,(1:numElements)');   
end

if chullCropOpt==1    
    if size(TR.ConnectivityList,2)>3 %ie tetrahedra
        FBtri = freeBoundary(TR);
        DT=delaunayTriangulation(TR.Points(unique(FBtri(:)),:));
        ti_DT= pointLocation(DT,QP);
        logicDT=~isnan(ti_DT);
        QP=QP(logicDT,:);    
    else
        warning('chullCropOpt not available for traingulated meshes');
        logicDT=true(size(QP,1),1);
    end
else
    logicDT=true(size(QP,1),1);
end

%% START LOOP

if waitbarOpt==1
    waitBarString='Finding points inside tesselation elements...';
    hw = waitbar(0,waitBarString);
    nUpdateSteps=25;
    waitbarSet=unique(round(linspace(1,numElements,nUpdateSteps)));
end

ti=nan(size(QP,1),1);
bc=nan(size(QP,1),numPerE);
for q=1:1:numElements
    
    if distCropOpt==1
        %Get current circum centre and radius
        TR_centre=TR_centre_all(q,:);        
        TR_rad_max=TR_rad_all(q);
        
        %Check distance with current points        
        QP_min=QP-TR_centre(ones(size(QP,1),1),:); %Query points with current element centre subtracted
        QP_dist=sqrt(sum(QP_min.^2,2)); %Distances of query points to current element centre        
        logicClose=QP_dist<=(TR_rad_max+toleranceMagnitude); 
        bcClose=nan(size(logicClose,1),numPerE);
        QP_test=QP(logicClose,:);
        if nnz(logicClose)>0
            [logicClose_sub,bcClose_sub]=isInsideTR(TR,QP_test,q*ones(size(QP_test,1),1),toleranceMagnitude);
            bcClose(logicClose,:)=bcClose_sub;            
            logicClose(logicClose)=logicClose_sub;        
        end             
    else
        [logicClose,bcClose]=isInsideTR(TR,QP,q*ones(size(QP,1),1),toleranceMagnitude);        
    end
    
    ti(logicClose)=q;
    
    bc(logicClose,:)=bcClose(logicClose,:);
    
    if waitbarOpt==1 && ismember(q,waitbarSet)
        waitbar(q/numElements,hw,waitBarString);
    end
end

if waitbarOpt==1
    close(hw);
end

TI=nan(size(logicDT));
TI(logicDT)=ti;

BC=nan(size(logicDT,1),numPerE);
BC(logicDT,:)=bc;

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
