function [varargout]=quad2tri(varargin)

% function [Ft,Vt,Ct]=quad2tri(Fq,Vq,triType,Cq)
% ------------------------------------------------------------------------
% This function converts quadrilateral patch data to triangular patch data.
% 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/11/03
% 2018/08/21 Added tolerances for edge and angle based slashing
%------------------------------------------------------------------------


%% Parse input
switch nargin
    case 2
        Fq=varargin{1};
        Vq=varargin{2};
        triType=[];
        Cq=[];
    case 3        
        Fq=varargin{1};
        Vq=varargin{2};
        triType=varargin{3};
        Cq=[];
    case 4        
        Fq=varargin{1};
        Vq=varargin{2};
        triType=varargin{3};
        Cq=varargin{4};        
    otherwise
        error('False number of input arguments');
end

if isempty(triType)
    triType='a';
end

%%
switch triType
    case 'b' %backward slash
        Ft=[Fq(:,[1 2 3]);Fq(:,[3 4 1]);];
        Vt=Vq;
        if size(Cq,1)==size(Fq,1) %Face colors
            Ct=repmat(Cq,[2,1]);
        else %Assume vertex colors
            Ct=Cq;
        end
    case 'f' %Forward slash
        Ft=[Fq(:,[1 2 4]);Fq(:,[2 3 4]);];        
        Vt=Vq;        
        if size(Cq,1)==size(Fq,1) %Face colors
            Ct=repmat(Cq,[2,1]);
        else %Assume vertex colors
            Ct=Cq;
        end
    case 'e'        
        D=patchEdgeLengths(Fq,Vq);
        edgeTolerance=min(D)./1000; %Tolerance level
        d1=sum((Vq(Fq(:,1),:)-Vq(Fq(:,3),:)).^2,2);
        d2=sum((Vq(Fq(:,2),:)-Vq(Fq(:,4),:)).^2,2);
        d=d1-d2;                 
        L1=d<-edgeTolerance; %Flip if difference is more than tolerance        
        Ft=[Fq(L1,[1 2 3]);Fq(L1,[3 4 1]); Fq(~L1,[1 2 4]);Fq(~L1,[2 3 4])];
        Vt=Vq;
        
        if size(Cq,1)==size(Fq,1) %Face colors
            Ct=[Cq(L1,:); Cq(L1,:); Cq(~L1,:); Cq(~L1,:)];
        else %Assume vertex colors
            Ct=Cq;
        end        
    case 'a'        
        angleTolerance=(pi/3)/1000; %Tolerance level
        F1=[Fq(:,[1 2 3]);Fq(:,[3 4 1])];
        F2=[Fq(:,[1 2 4]);Fq(:,[2 3 4])];
        [a1]=patchEdgeAngles(F1,Vq);
        a1=reshape(a1,[size(Fq,1),6]);
        d1=sum(abs(a1-(pi/3)),2);
        
        [a2]=patchEdgeAngles(F2,Vq);
        a2=reshape(a2,[size(Fq,1),6]);
        d2=sum(abs(a2-(pi/3)),2);
        d=d1-d2;                 
        L1=d<-angleTolerance; %Flip if difference is more than tolerance
        Ft=[Fq(L1,[1 2 3]);Fq(L1,[3 4 1]); Fq(~L1,[1 2 4]);Fq(~L1,[2 3 4])];          
        Vt=Vq;
        
        if size(Cq,1)==size(Fq,1) %Face colors
            Ct=[Cq(L1,:); Cq(L1,:); Cq(~L1,:); Cq(~L1,:)];
        else %Assume vertex colors
            Ct=Cq;
        end                
    case 'x' %Cross type
        Vm=zeros(size(Fq,1),size(Vq,2));
        for q=1:1:size(Vq,2)
            X=Vq(:,q);
            FX=X(Fq);
            if size(Fq,1)==1 %Treat special case of single face
                FX=FX';
            end
            Vm(:,q)=mean(FX,2);
        end
        %Join point sets
        Vt=[Vq;Vm];
        
        indVm=(size(Vq,1)+1):size(Vt,1);
        %Create faces
        Ft=[Fq(:,1) Fq(:,2) indVm(:);...
            Fq(:,2) Fq(:,3) indVm(:);...
            Fq(:,3) Fq(:,4) indVm(:);...
            Fq(:,4) Fq(:,1) indVm(:)];
        
        if size(Cq,1)==size(Fq,1) %Face colors
            Ct=repmat(Cq,[4,1]);
        else %Assume vertex colors
            Cm=zeros(size(Vm,1),size(Cq,2));
            for q=1:1:size(Cq,2)
                c=Cq(:,q);
                Cm(:,q)=mean(c(Fq),2);
            end
            Ct=[Cq; Cm];
        end
end
    
varargout{1}=Ft;
varargout{2}=Vt;
varargout{3}=Ct;

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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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
