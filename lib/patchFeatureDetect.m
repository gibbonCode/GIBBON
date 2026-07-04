function [G]=patchFeatureDetect(F,V,a)

%%
%Get face normals
N=patchNormal(F,V);

%%

% cFigure; hold on;
% gpatch(F,V,'w','none',0.5,1);
% hp=gpatch(F,V,ones(size(F,1),1),'none',ones(size(F,1),1),1);
% hp.FaceColor='flat';
% axisGeom; camlight headlight; 
% colormap(gjet(4)); %icolorbar;
% axis off; 
% gdrawnow; 

%%

G=zeros(size(F,1),1);

doneFlag=0;
startNew=true;
groupLabel=1;


[C]=patchConnectivity(F,V,'ff');
IND_FF=C.face.face;
% [~,~,IND_FF]=tesIND(F,V,0);
logicValid=IND_FF~=0;
IND_FF(~logicValid)=NaN;
nff=size(IND_FF,2);

logicDone=false(size(F,1),1); 
%%
while doneFlag==0     
    if startNew %Start new group        
        L=false(size(F,1),1); %Reset logic group
        Lp=false(size(F,1),1); %Reset previous logic
        indStart=find(~logicDone,1);                 
        L(indStart)=1; %The current face logic
        startNew=false; %Flip start switch
        Ln=L; %The new face logic        
    else        
        IND_FF_cand=IND_FF(Ln,:); %Entries for candidate faces
        
        %Dot dot product
        NX_cand=IND_FF_cand; 
        
        NX_cand(logicValid(Ln,:));
        IND_FF_cand(logicValid(Ln,:));
        N(IND_FF_cand(logicValid(Ln,:)),1);
        
        
        NX_cand(logicValid(Ln,:))=N(IND_FF_cand(logicValid(Ln,:)),1);
        NY_cand=IND_FF_cand; 
        NY_cand(logicValid(Ln,:))=N(IND_FF_cand(logicValid(Ln,:)),2);
        NZ_cand=IND_FF_cand; 
        NZ_cand(logicValid(Ln,:))=N(IND_FF_cand(logicValid(Ln,:)),3);
        
        %Compute angle
        A_cand=acos((N(Ln,1*ones(nff,1)).*NX_cand)+(N(Ln,2*ones(nff,1)).*NY_cand)+(N(Ln,3*ones(nff,1)).*NZ_cand));
        
        Lp = L; %The previous face logic
        L(IND_FF_cand(A_cand<=a))=1; %The current face logic
        Ln = L & ~Lp; %The new face logic
                
%         C=zeros(size(L));
%         C(L )=1;
%         C(Lp)=2;
%         C(Ln)=3;
%         set(hp,'CData',double(C));        
%         set(hp,'FaceVertexAlphaData',double(C)); 
%         drawnow;     
    end
    
    if nnz(Lp)==nnz(L) %If the group has not grown
        startNew=true;        
        G(L)=groupLabel;     
        logicDone(L)=1;
        groupLabel=groupLabel+1;                
    end
    
    if all(logicDone)
        doneFlag=1;                
    end            
end

%%

% %% Prepare output
% varargout{1}=G;
% varargout{2}=G_iter;
% 
% if nargout==3
%     switch outputType
%         case 'array'
%             groupSize=sum(G,1);
%         case 'label'
%             groupSize=zeros(1,max(G(:)));
%             for q=1:1:max(G(:))
%                 groupSize(q)=nnz(G==q);
%             end
%     end
%     varargout{3}=groupSize;
% end

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
