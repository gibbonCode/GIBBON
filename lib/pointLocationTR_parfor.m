function [TI,BC]=pointLocationTR_parfor(TR,QP,distCropOpt,chullCropOpt)

% function [TI,BC]=pointLocationTR(TR,QP,distCropOpt,chullCropOpt)
% ------------------------------------------------------------------------
%
% This function finds the triangles or tetrahedrons in which the query
% points QP are contained. In addition it outputs the barycentric
% coordinates.
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/04/01
%------------------------------------------------------------------------

E=TR.ConnectivityList;
numPerE=size(E,2);
numElements=size(E,1);

if distCropOpt==1
    [TR_centre_all,TR_rad_all] = circumcenter(TR,(1:numElements)');   
end

if chullCropOpt==1
    FBtri = freeBoundary(TR);
    DT=delaunayTriangulation(TR.Points(unique(FBtri(:)),:));
    ti_DT= pointLocation(DT,QP);
    
    logicDT=~isnan(ti_DT);
    QP=QP(logicDT,:);    
else
    logicDT=true(size(QP,1),1);
end

%% START LOOP

siz1=size(QP,1);
siz2=numElements;
SP_ti=sparse(siz1,siz2);
SP_bc1=sparse(siz1,siz2);
SP_bc2=sparse(siz1,siz2);
SP_bc3=sparse(siz1,siz2);
SP_bc4=sparse(siz1,siz2);

%%
try
    matlabpool close;
catch exception
    disp('--> No existing pools to close');
end

%Define profile and number of workers
defaultProfile = parallel.defaultClusterProfile;
myCluster = parcluster(defaultProfile);
myCluster.NumWorkers=2; 

%Open matlab pool
matlabpool(myCluster, 'open');

%%
parfor q=1:numElements

    ti=nan(size(QP,1),1);
    bc=nan(size(QP,1),numPerE);

    if distCropOpt==1
        %Get current circum centre and radius
        TR_centre=TR_centre_all(q,:);        
        TR_rad_max=TR_rad_all(q);
        
        %Check distance with current points        
        QP_min=QP-TR_centre(ones(size(QP,1),1),:); %Query points with current element centre subtracted
        QP_dist=sqrt(sum(QP_min.^2,2)); %Distances of query points to current element centre        
        logicClose=QP_dist<=(TR_rad_max+eps(TR_rad_max)); 
        bcClose=nan(size(logicClose,1),numPerE);
        QP_test=QP(logicClose,:);
        if nnz(logicClose)>0
            [logicClose_sub,bcClose_sub]=isInsideTR(TR,QP_test,q*ones(size(QP_test,1),1));
            bcClose(logicClose,:)=bcClose_sub;            
            logicClose(logicClose)=logicClose_sub;        
        end             
    else
        [logicClose,bcClose]=isInsideTR(TR,QP,q*ones(size(QP,1),1));        
    end    
    ti(logicClose)=q;    
    bc(logicClose,:)=bcClose(logicClose,:);    
    
    ti(isnan(ti))=0;
    bc(isnan(bc))=0;
    
    SP_ti(:,q)=ti;
    SP_bc1(:,q)=bc(:,1);
    SP_bc2(:,q)=bc(:,2);
    SP_bc3(:,q)=bc(:,3);
    SP_bc4(:,q)=bc(:,4);
end

[~,J_max]=max(SP_ti,[],2);

J_max=full(J_max);
[indMax]=sub2ind(size(SP_ti),[1:size(SP_ti,1)]',J_max);

ti=full(SP_ti(indMax));
bc=[full(SP_bc1(indMax)) full(SP_bc2(indMax)) full(SP_bc3(indMax)) full(SP_bc4(indMax))];

TI=nan(size(logicDT));
TI(logicDT)=ti;

BC=nan(size(logicDT,1),numPerE);
BC(logicDT,:)=bc;

end
 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
