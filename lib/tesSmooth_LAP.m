function [P]=tesSmooth_LAP(TES,V,IND_V,cPar)

%% CONTROL PARAMETERS

if isfield(cPar,'LambdaSmooth')
    LambdaSmooth=cPar.LambdaSmooth;
else
    LambdaSmooth=0.5; %DEFAULT
end

if isfield(cPar,'n')
    nMax=cPar.n;
else
    nMax=1; %DEFAULT
end

if isfield(cPar,'RigidConstraints')
    indRigid=cPar.RigidConstraints;
else
    indRigid=[]; %DEFAULT
end

if isfield(cPar,'Tolerance')
    SSQD_Tol=cPar.Tolerance;
else
    SSQD_Tol=[]; %DEFAULT
end

%%

if ~isempty(SSQD_Tol)
    SSQD_old=[];
    SSQD_ratio=0;
end

if isempty(IND_V)
    [~,IND_V]=patchIND(TES,V,2);
end
logicValid=IND_V>0;
indNoneValid=find(sum(logicValid,2)==0);

%%

nDims=size(V,2); %Number of dimensions

P=V;
PP=V; 
Q=V;
for qIter=1:nMax 
        
    %% SIMPLE LAPLACIAN SMOOTHENING
    
    %Loop for all dimensions
    for qDim=1:1:nDims
        Xp=NaN(size(IND_V,1),size(IND_V,2));
        Xp(logicValid)=P(IND_V(logicValid),qDim);
        Xp=gnanmean(Xp,2);       
        Xp(indNoneValid)=V(indNoneValid,qDim);
        PP(:,qDim)=Xp;
    end
    P=P+LambdaSmooth.*(PP-P);
    
    %%
        
    %Put back constrained points
    if ~isempty(indRigid)
       P(indRigid,:)=V(indRigid,:);
    end
    
    if ~isempty(SSQD_Tol)
        %Compute sum of squared differences with respect to previous iteration
        SSQD_new=gnansum((P(:)-Q(:)).^2);
        if ~isempty(SSQD_old)
            SSQD_ratio=SSQD_new./SSQD_old;            
        end
        
        %Store current metrics
        Q=P;
        SSQD_old=SSQD_new;        
        if abs(1-SSQD_ratio)<=SSQD_Tol
           break %STOP SMOOTHING LOOP IF TOLERANCE IS REACHED 
        end
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
