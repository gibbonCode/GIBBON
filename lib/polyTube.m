function [varargout]=polyTube(Vg,optionStruct)

% function [varargout]=polyTube(Vg,optionStruct)
% ------------------------------------------------------------------------
%
%
%
%
% Change log (started 2019/03/28): 
% 2019/03/28 Added spatially varying radii
% 2019/03/29 Added capped ends option (for triangulated patch type) 
% ------------------------------------------------------------------------

%% Parse input

numSteps=size(Vg,1);

%Create default option structure
dPath=pathLength(Vg);
defaultOptionStruct.r=min(diff(dPath))/2; %Half the minimum point spacing
defaultOptionStruct.C=(1:1:numSteps)'; %Point indices
defaultOptionStruct.nr=3;
defaultOptionStruct.patchType='quad';
defaultOptionStruct.closeOpt=0;
defaultOptionStruct.fixOpt=eps;

%Parse input structure (complete with default)
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

%Access options
C=optionStruct.C;
nr=optionStruct.nr;
patchType=optionStruct.patchType;
closeOpt=optionStruct.closeOpt;
r=optionStruct.r;
fixOpt=optionStruct.fixOpt;

if numel(r)==1
    r=r.*ones(numSteps,1);
end
r=r(:);

indRepPath=repmat(1:1:numSteps,nr,1);
indRepPath=indRepPath(:);
colorCirc=repmat((1:1:numSteps),nr,1)';

%% Create coordinates for circles

%Create single circle
t=linspace(0,2*pi,nr+1)'; t=t(1:end-1);
Vc=[sin(t) cos(t) zeros(size(t))];

%Replicate circle coordinates numSteps times
VC=repmat(Vc,[numSteps,1]);

%Assign radii to all circles by scaling
RC=r(indRepPath);
VC=VC.*RC(:,ones(1,size(VC,2)));

%% Construct rotation tensors

%Allong path edge vectors
Ue=Vg(2:end,:)-Vg(1:end-1,:); 
Ue(end+1,:)=Ue(end,:);
Ue=vecnormalize(Ue);

%Create allong path vertex vectors
U=(Ue(2:end,:)+Ue(1:end-1,:))/2; %Average of adjacent
U=[Ue(1,:); U];
U=vecnormalize(U);

%Define 3rd direction
e3=vecnormalize(Vg(2,:)-Vg(1,:));

%Set initial 2nd dir. as most orthogonal 
I=eye(3,3);
d=dot(e3(ones(3,1),:),I,2);
[~,indBest]=min(d);
e2=I(indBest,:);
e1=vecnormalize(cross(e3,e2));
e2=vecnormalize(cross(e3,e1));

%Define rotation matrices for all points allong curve
R_curve=zeros(size(Vg,1),9);
R1=[e1; e2; e3]';
R_curve(1,:)=R1([1 4 7 2 5 8 3 6 9]);

for q=2:size(Vg,1)
    a=U(q-1,:); b=U(q,:);    
    theta=real(acos(dot(a,b))); %Complex if dot product is out of range [-1.1] due to precission issues
    w=vecnormalize(cross(b,a));    
    if norm(w)>0.5
        [Rn]=vecAngle2Rot(theta,w);
        R1=Rn'*R1;
    end      
    R_curve(q,:)=R1([1 4 7 2 5 8 3 6 9]);
end

R_mat=R_curve(indRepPath,:);

%% Rotate circles

[VC]=vectorTensorProductArray(VC,R_mat);

%% Translate circles

VS=VC+Vg(indRepPath,:);

%% Create coordinate grids for mesh formation

X=reshape(VS(:,1),[nr,numSteps])';
Y=reshape(VS(:,2),[nr,numSteps])';
Z=reshape(VS(:,3),[nr,numSteps])';

%% Fix self intersections

if fixOpt==1
    
    D=distND(VS,Vg);
    logicSelf=repmat(indRepPath,1,size(Vg,1))==repmat(1:1:size(Vg,1),size(VS,1),1);
    D(logicSelf)=NaN;
    
    logicInValid = any(D<RC,2);
    logicInValidGrid=reshape(logicInValid,[nr,numSteps])';
    
%     %Force start and end to zero    
%     if any(logicInValidGrid(1,:))
%         v=[X(1,:)' Y(1,:)' Z(1,:)']; %Current curve
%         v=evenlySampleCurve(v(~logicInValidGrid(1,:),:),nr,'spline',0);        
%         X(1,:)=v(:,1);
%         Y(1,:)=v(:,2);
%         Z(1,:)=v(:,3);
%     end
%     
%     if any(logicInValidGrid(end,:))
%         v=[X(end,:)' Y(end,:)' Z(end,:)']; %Current curve
%         v=evenlySampleCurve(v(~logicInValidGrid(end,:),:),nr,'spline',0);
%         X(end,:)=v(:,1);
%         Y(end,:)=v(:,2);
%         Z(end,:)=v(:,3);
%     end
%     
%     logicInValidGrid(1,:)=0;
%     logicInValidGrid(end,:)=0;
    
    indColumnFix=find(any(logicInValidGrid,1));
    
    % for q=indColumnFix
    %     logicInvalidNow=logicInValidGrid(:,q);
    %
    %     dd=diff(logicInvalidNow);
    %     indStart=find([0; dd==1]);
    %     indEnd=find([dd==-1; 0]);
    %
    %     v=[X(:,q) Y(:,q) Z(:,q)]; %Current curve
    %
    %     for qs=1:1:numel(indStart)
    %         indStartNow=indStart(qs);
    %         indEndNow=indEnd(qs);
    %         ind=indStartNow:indEndNow;
    %         v(ind,:)=linspacen(v(indStartNow,:),v(indEndNow,:),numel(ind))';
    %
    %         d=pathLength(v);
    %         di=d(logicInvalidNow,:);
    %         for qq=1:1:3
    %             v(logicInvalidNow,qq)=interp1(d(~logicInvalidNow),v(~logicInvalidNow,qq),di,'spline');
    %         end
    %
    %         d=pathLength(v);
    %         di=linspace(d(indStartNow),d(indEndNow,:),numel(ind))';
    %         numel(di)
    %         nnz(logicInvalidNow)
    %         for qq=1:1:3
    %             v(ind,qq)=interp1(d,v(:,qq),di,'linear');
    %         end
    %     end
    %     X(:,q)=v(:,1);
    %     Y(:,q)=v(:,2);
    %     Z(:,q)=v(:,3);
    % end
    
    for q=1:1:size(X,2)
        v=[X(:,q) Y(:,q) Z(:,q)]; %Current curve
        v=evenlySampleCurve(v(~logicInValidGrid(:,q),:),size(v,1),'spline',0);
        X(:,q)=v(:,1);
        Y(:,q)=v(:,2);
        Z(:,q)=v(:,3);
    end
    
end

[Fs,Vs,~,Cs_ind]=grid2patch(X,Y,Z,colorCirc,[0 1]);

Fs=fliplr(Fs); 

%% Create colors
Cs_d=repmat(dPath,[1,nr]);
Cs_d=Cs_d(:);

% Cs_rgb=vecnormalize(abs(U(Cs_ind,:)));
Cs_rgb=vecnormalize(faceToVertexMeasure(Fs,Vs,abs(Vs(Fs(:,3),:)-Vs(Fs(:,2),:))));

% Cs_C=repmat(C,[1,nr]);
% Cs_C=Cs_C(:);

Cs_C=C(Cs_ind);

%%

switch patchType
    case 'quad'
        if closeOpt
            error('Capping cylinder is currently only supported for the tri patch type');
        end
    case 'tri'        
        [Fs,Vs,Cs_C]=quad2tri(Fs,Vs,'a',Cs_C);

        if closeOpt
            indTop=numSteps:numSteps:size(Vs,1);
            [Ft,Vt]=regionTriMesh3D({Vs(indTop,:)},[],0);
            
            indBottom=1:numSteps:size(Vs,1);
            [Fb,Vb]=regionTriMesh3D({Vs(indBottom,:)},[],0);
            
            [Fs,Vs,Cs]=joinElementSets({Fs,Ft,Fb},{Vs,Vt,Vb});
            [Fs,Vs,ind1]=mergeVertices(Fs,Vs);
            
            Cs_C=[Cs_C; C(size(C,1).*ones(size(Vt,1),1),:); C(ones(size(Vb,1),1),:)];
            Cs_C=Cs_C(ind1,:);
            
            %Fix normal directions for caps if needed
            E1=patchEdges(Fs(Cs==1,:),0);
            for q=2:1:3
                Eb=patchBoundary(Fs(Cs==q,:),Vs);
                %If cap boundary edges are a member of the tube mesh edges the cap
                %faces are pointing the wrong way
                if any(isrowmember(Eb,E1,0))
                    Fs(Cs==q,:)=fliplr(Fs(Cs==q,:));
                end
            end
        end
end

%% Collect output

varargout{1}=Fs;
varargout{2}=Vs;
varargout{3}=Cs_C;
varargout{4}=Cs_rgb;
varargout{5}=Cs_d;
 
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
