function [F,V]=subQuad(F,V,n)

% function [Fs,Vs]=subQuad(F,V,n)
% ------------------------------------------------------------------------
% Sub-devides the quadrilateral faces defined by the patch format data F
% (faces) and V (vertices). Each face is split n times. 
% 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2010/06/01 Created
% 2017/11/29 Fixed single face input related bug (requires transpose of
% arrays)
% ------------------------------------------------------------------------

%%
if size(V,2)==2
    V(:,3)=0;
end
if n>0
    for qIter=1:1:n
        
        numV=size(V,1);
        X=V(:,1); Y=V(:,2); Z=V(:,3);
        
        %% DERIVE EDGE INDICES AND FACE-EDGE INDEX MATRIX
        
        %Format of column index in F
        EColumnInd=[(1:size(F,2)); (1:size(F,2))];
        EColumnInd=[EColumnInd(2:end) EColumnInd(1)];
        
        %Derive edges matrix
        E=F(:,EColumnInd)'; %Use index into F to create edges matrix
        E=reshape(E,2,numel(E)/2)';
        
        E=sort(E,2); %Sort edge order
        [E,~,ind2] = unique(E,'rows'); %Removing double edges, i.e. [1  4] = [4  1]
        
        Fe=reshape(1:numel(F),size(F,2),size(F,1))';
        Fe=ind2(Fe);
        
        if size(F,1)==1
            Fe=Fe';
        end
   
        %% Calculate mid-face vertices
                
        if size(F,1)==1
            XF=X(F)'; YF=Y(F)'; ZF=Z(F)';
        else
            XF=X(F); YF=Y(F); ZF=Z(F);
        end
        V_midFace=[mean(XF,2) mean(YF,2) mean(ZF,2)];
        
        %% Calculate mid-edge vertices
        
        XE=X(E); YE=Y(E); ZE=Z(E);
        V_midEdge=[mean(XE,2) mean(YE,2) mean(ZE,2)];        
        
        %% Create new faces matrix
        
        Vs=[V;V_midEdge;V_midFace];
        startIndMidEdge=numV;
        startIndMidFace=startIndMidEdge+size(V_midEdge,1);
        indMidFace=startIndMidFace+1:size(Vs,1);        
        
        Fs=repmat(F,4,1);
        for q=1:1:4
            startInd=1+(q-1)*size(F,1);
            endInd=startInd-1+size(F,1);
            
            if q==1
                ind4=4;
            else
                ind4=q-1;
            end

            fs=[F(:,q) Fe(:,q)+numV indMidFace(:) Fe(:,ind4)+numV];

            Fs(startInd:endInd,:)=fs;
        end
        
        %% Override input
        F=Fs;
        V=Vs;
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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
