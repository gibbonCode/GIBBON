function [F,V]=polyExtrude(Vc,cPar)

%% COMPUTER CURVE LENGTH
D=max(pathLength(Vc)); %Compute curve length for point sampling
numPoints=size(Vc,1);

%% PARSE cPar

%Check what mode
if isfield(cPar,'pointSpacing') && isfield(cPar,'numSteps')
    error('Either specify pointSpacing or numSteps, not both.');
end

%Check extrudeMode, numSteps and pointSpacing settings
if isfield(cPar,'numSteps')
    extrudeMode=1;
else
    extrudeMode=0;
    if ~isfield(cPar,'pointSpacing')
        cPar.pointSpacing=[]; %Default is average point spacing
    end
    if isempty(cPar.pointSpacing)
        cPar.pointSpacing=D/numPoints; %Default is average point spacing
    end
end

%Check depth
if ~isfield(cPar,'depth');
    error('cPar.depth was not specified.');
end

%Check dir
if ~isfield(cPar,'dir');
    cPar.dir=0; %Default symmetric
end

%Check patchType
if ~isfield(cPar,'patchType');
    cPar.patchType='tri'; %Default triangles
end

%Check direction vector
if ~isfield(cPar,'n');
    [R_fit]=pointSetPrincipalDir(Vc);
    cPar.n=R_fit(:,3); %Default normal direction to polygon
end
cPar.n=vecnormalize(cPar.n);
cPar.n=cPar.n(:)';

%Check closeLoopOpt
if ~isfield(cPar,'closeLoopOpt');
    cPar.closeLoopOpt=0; %Default off
end

%% Extruding the skethed profile

switch extrudeMode
    case 1 %Use number of points
        
    case 0 %Resample curve and use pointSpacing
        %Set point spacing
        cPar.numSteps=round(cPar.depth./cPar.pointSpacing);
        
        %Resampling sketch with desired spacing
        D=max(pathLength(Vc)); %Computer curve length for point sampling
        n=round(D./cPar.pointSpacing); %Determine number of points based on curve length
        interpMethod='linear';
        [Vc]=evenlySampleCurve(Vc,n,interpMethod,cPar.closeLoopOpt); %Resampling curve
end

%Create coordinates in extrusion direction
switch cPar.dir
    case 0
        V_add=(cPar.depth/2).*cPar.n;
        Vc_start=Vc-V_add(ones(size(Vc,1),1),:);
        Vc_end=Vc+V_add(ones(size(Vc,1),1),:);
    case 1
        Vc_start=Vc;
        V_add=cPar.depth.*cPar.n;
        Vc_end=Vc+V_add(ones(size(Vc,1),1),:);
    case -1
        Vc_start=Vc;
        V_add=cPar.depth.*cPar.n;
        Vc_end=Vc-V_add(ones(size(Vc,1),1),:);
end

[F,V]=polyLoftLinear(Vc_start,Vc_end,cPar);

% X=linspacen(Vc_start(:,1),Vc_end(:,1),cPar.numSteps)';
% Y=linspacen(Vc_start(:,2),Vc_end(:,2),cPar.numSteps)';
% Z=linspacen(Vc_start(:,3),Vc_end(:,3),cPar.numSteps)';
% c=(1:1:size(Z,1))';
% C=c(:,ones(1,size(Z,2)));
% 
% %Create quad patch data
% [F,V,C] = surf2patch(X,Y,Z,C);
% 
% %Close patch if required
% if cPar.closeLoopOpt
%     I=[(2:size(Z,1))' (2:size(Z,1))' (1:size(Z,1)-1)' (1:size(Z,1)-1)'];
%     J=[ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1)];
%     F_sub=sub2ind(size(Z),I,J);
%     F=[F;F_sub];
%     [C]=vertexToFaceMeasure(F,C);
%     C(end-size(F_sub,1):end,:)=C(end-size(F_sub,1):end,:)+0.5; 
% else
%     [C]=vertexToFaceMeasure(F,C);
% end
% C=round(C);
% 
% switch cPar.patchType
%     case 'quad' 
%     case 'tri_slash' %Convert quads to triangles by slashing
%         F=[F(:,1) F(:,3) F(:,2); F(:,1) F(:,4) F(:,3)];
%     case 'tri' %Convert quads to approximate equilateral triangles
% 
%         logicSlashType=repmat(iseven(C),2,1);
%         
%         Xi=X;
%         x=X(1,:);
%         dx=diff(x)/2;
%         dx(end+1)=(x(1)-x(end))/2;
%         for q=2:2:size(X,1)
%             X(q,:)=X(q,:)+dx;
%         end
%         if ~cPar.closeLoopOpt
%             X(:,1)=Xi(:,1);
%             X(:,end)=Xi(:,end);
%         end
%         
%         Yi=Y;
%         y=Y(1,:);
%         dy=diff(y)/2;
%         dy(end+1)=(y(1)-y(end))/2;
%         for q=2:2:size(Y,1)
%             Y(q,:)=Y(q,:)+dy;
%         end
%         if ~cPar.closeLoopOpt
%             Y(:,1)=Yi(:,1);
%             Y(:,end)=Yi(:,end);
%         end
%         
%         Zi=Z;
%         z=Z(1,:);
%         dz=diff(z)/2;
%         dz(end+1)=(z(1)-z(end))/2;
%         for q=2:2:size(Z,1)
%             Z(q,:)=Z(q,:)+dz;
%         end
%         if ~cPar.closeLoopOpt
%             Z(:,1)=Zi(:,1);
%             Z(:,end)=Zi(:,end);
%         end
%         
%         V=[X(:) Y(:) Z(:)];
%         
%         F1=[F(:,1) F(:,3) F(:,2); F(:,1) F(:,4) F(:,3)]; 
%         F2=[F(:,1) F(:,4) F(:,2); F(:,2) F(:,4) F(:,3)];
%         F=[F1(~logicSlashType,:);F2(logicSlashType,:)];
%         
% %         C=repmat(C,2,1);
% %         C=[C(~logicSlashType,:);C(logicSlashType,:)];
%         
% end
% 
