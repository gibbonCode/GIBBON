function [F_implant,V_implant,C_implant,curveSet]=parHipImplant(inputStruct)

% function [F_implant,V_implant,C_implant,curveSet]=parHipImplant(inputStruct)
% ------------------------------------------------------------------------
%
% 
% ------------------------------------------------------------------------

%% Parse input structure

defaultInputStruct.ballRadius=20;
defaultInputStruct.stickRadius=7;
defaultInputStruct.stickLength=21;
defaultInputStruct.stickLengthStraight=defaultInputStruct.stickLength-6;
defaultInputStruct.neckRadius=15;
defaultInputStruct.neckEllipseScale=2;
defaultInputStruct.collarThickness=3; 
defaultInputStruct.loftOffset=20;
defaultInputStruct.loftLenght=40;
defaultInputStruct.stemRadius=8;
defaultInputStruct.stemLength=50;
defaultInputStruct.stemAngle=0.25*pi;
defaultInputStruct.pointSpacing=2;

[inputStruct]=structComplete(inputStruct,defaultInputStruct,1); %Complement provided with default if missing or empty

ballRadius=inputStruct.ballRadius;
stickRadius=inputStruct.stickRadius;
stickLength=inputStruct.stickLength;
stickLengthStraight=inputStruct.stickLengthStraight;
neckRadius=inputStruct.neckRadius;
neckEllipseScale=inputStruct.neckEllipseScale;
collarThickness=inputStruct.collarThickness;
loftOffset=inputStruct.loftOffset;
loftLenght=inputStruct.loftLenght;
stemRadius=inputStruct.stemRadius;
stemLength=inputStruct.stemLength;
stemAngle=inputStruct.stemAngle;
pointSpacing=inputStruct.pointSpacing;

%%
xc=cos(asin(stickRadius/ballRadius))*ballRadius;
dx=ballRadius-xc;
Qz=euler2DCM([0 0 -stemAngle]);
n2=[1 0 0];
n2=n2*Qz;
pTip=[ballRadius+stickLength+loftOffset+collarThickness 0 0];
pTip=pTip+(loftLenght+stemLength).*n2;

p1=[ballRadius+stickLength+collarThickness 0 0];
p2=[ballRadius+stickLength+collarThickness+loftOffset 0 0]+(loftLenght.*n2);
n1=[1 0 0];

%%
nRef=1;
while 1
    [F_ball,V_ball]=geoSphere(nRef,ballRadius);
    pointSpacingNow=mean(patchEdgeLengths(F_ball,V_ball));
    if pointSpacingNow<=pointSpacing
        break
    else
        nRef=nRef+1;
    end
end

logicRight=V_ball(:,1)>xc;
logicCut=any(logicRight(F_ball),2);
logicCut=triSurfLogicSharpFix(F_ball,logicCut,3);
F_ball=F_ball(~logicCut,:);
[F_ball,V_ball]=patchCleanUnused(F_ball,V_ball);
Eb_ball=patchBoundary(F_ball); 
indB=unique(Eb_ball(:));
[T,P,R] = cart2sph(V_ball(:,2),V_ball(:,3),V_ball(:,1));
P(indB)=atan2(xc,ballRadius*sin(acos(xc./ballRadius)));
[V_ball(:,2),V_ball(:,3),V_ball(:,1)] = sph2cart(T,P,R);

indBallHole=edgeListToCurve(Eb_ball);
indBallHole=indBallHole(1:end-1);

pointSpacingBall=mean(sqrt(sum(diff(V_ball(indBallHole,:),1,1).^2,2)));

%%

numStepsExtrude=ceil(stickLength./pointSpacingBall);
numStepsExtrude=numStepsExtrude+double(iseven(numStepsExtrude));

clear cParExtrude;
cParExtrude.depth=stickLength+dx; 
cParExtrude.patchType='tri'; 
cParExtrude.dir=1;
cParExtrude.n=[1 0 0];
cParExtrude.closeLoopOpt=1; 
cParExtrude.numSteps=numStepsExtrude;
[F_stick,V_stick]=polyExtrude(V_ball(indBallHole,:),cParExtrude);
F_stick=fliplr(F_stick);
logicLow=V_stick(:,2)<0;
logicRight=V_stick(:,1)>(ballRadius+stickLengthStraight);

w=V_stick(:,1)-(ballRadius+stickLengthStraight);
w=(w./max(w(:)));
[T,R,Z] = cart2pol(V_stick(:,2),V_stick(:,3),V_stick(:,1));
R(logicRight)=R(logicRight)+(neckRadius-stickRadius).*w(logicRight);
[V_stick(:,2),V_stick(:,3),V_stick(:,1)] = pol2cart(T,R,Z);

w=V_stick(:,1)-(ballRadius+stickLengthStraight);
w=w./max(w(:));
w=w.*(neckEllipseScale-1);
w=w+1;
V_stick(logicRight & logicLow,2)=V_stick(logicRight & logicLow,2).*w(logicRight & logicLow);

indEndStick=numStepsExtrude:numStepsExtrude:size(V_stick,1);
pointSpacingStick=mean(patchEdgeLengths(F_stick,V_stick));

%%

vn2=V_stick(indEndStick,:);
ve=[vn2(2:end,:); vn2(1,:)]-vn2;
n=[1 0 0];
vd=vecnormalize(cross(n(ones(size(ve,1),1),:),ve));
vn2(:,1)=vn2(:,1)+collarThickness;

nc=ceil((pi*collarThickness/2)/pointSpacingStick);
if nc<4
    nc=4;
end
vc=linspacen(V_stick(indEndStick,:),vn2,nc);
X=squeeze(vc(:,1,:));
Y=squeeze(vc(:,2,:));
Z=squeeze(vc(:,3,:));
t=repmat(linspace(-1,1,nc),size(Z,1),1);
a=acos(t);
X=X+collarThickness/2.*sin(a).*repmat(vd(:,1),1,nc);
Y=Y+collarThickness/2.*sin(a).*repmat(vd(:,2),1,nc);
Z=Z+collarThickness/2.*sin(a).*repmat(vd(:,3),1,nc);

[F_collar,V_collar]=grid2patch(X,Y,Z,[],[1 0 0]);
[F_collar,V_collar]=quad2tri(F_collar,V_collar,'a');

%%

[F_head,V_head,C_head]=joinElementSets({F_ball,F_stick,F_collar},{V_ball,V_stick,V_collar});
[F_head,V_head,~,indFix]=mergeVertices(F_head,V_head);
Eb_head=patchBoundary(F_head);
indBallHole=indFix(indBallHole);

%%

nRef=1;
while 1
    [F_tip,V_tip]=hemiSphereMesh(nRef,stemRadius,0);
    pointSpacingNow=mean(patchEdgeLengths(F_tip,V_tip));
    if pointSpacingNow<=pointSpacing
        break
    else
        nRef=nRef+1;
    end
end

Qy=euler2DCM([0 -0.5*pi 0]);
V_tip=V_tip*Qy;
V_tip=V_tip*Qz;
V_tip=V_tip+pTip(ones(size(V_tip,1),1),:);
Eb_tip=patchBoundary(F_tip);
indBoundaryCurveTip=edgeListToCurve(Eb_tip);
indBoundaryCurveTip=indBoundaryCurveTip(1:end-1);
pointSpacingTip=mean(patchEdgeLengths(F_tip,V_tip));

%%
clear cParExtrude;
numStepsExtrude=ceil(stemLength/pointSpacingTip);
numStepsExtrude=numStepsExtrude+double(iseven(numStepsExtrude));
cParExtrude.depth=stemLength; 
cParExtrude.patchType='tri'; 
cParExtrude.dir=1;
cParExtrude.n=-n2;
cParExtrude.closeLoopOpt=1; 
cParExtrude.numSteps=numStepsExtrude;
[F_stem_straight,V_stem_straight]=polyExtrude(V_tip(indBoundaryCurveTip,:),cParExtrude);
F_stem_straight=fliplr(F_stem_straight);
[F_stem,V_stem,C_stem]=joinElementSets({F_stem_straight,F_tip},{V_stem_straight,V_tip});
[F_stem,V_stem,~,indFix]=mergeVertices(F_stem,V_stem);
Eb_stem=patchBoundary(F_stem);
indBoundaryCurveTip=indFix(indBoundaryCurveTip+size(V_stem_straight,1));

%%
if size(Eb_head,1)<size(Eb_stem,1)
    numEdgesNeeded=size(Eb_stem,1);
    [F_head,V_head,Eb_head,C_head]=triSurfSplitBoundary(F_head,V_head,Eb_head,numEdgesNeeded,C_head);
else
    numEdgesNeeded=size(Eb_head,1);
    [F_stem,V_stem,Eb_stem,C_stem]=triSurfSplitBoundary(F_stem,V_stem,Eb_stem,numEdgesNeeded,C_stem);
end

indBoundaryCurve_head=edgeListToCurve(Eb_head);
indBoundaryCurve_head=indBoundaryCurve_head(1:end-1);
indBoundaryCurve_stem=edgeListToCurve(Eb_stem);
indBoundaryCurve_stem=indBoundaryCurve_stem(1:end-1);

%%

V_loft1=V_head(indBoundaryCurve_head,:);
[~,indMax]=max(V_loft1(:,3));
if indMax>1
    V_loft1=V_loft1([indMax:size(V_loft1,1) 1:indMax-1],:);
end
V_loft2=V_stem(indBoundaryCurve_stem,:);
[~,indMax]=max(V_loft2(:,3));
if indMax>1
    V_loft2=V_loft2([indMax:size(V_loft2,1) 1:indMax-1],:);
end
V_loft2=flipud(V_loft2);

% f=0.0005;
% pp=1;
% pointSpacing=mean(patchEdgeLengths(F_head,V_head));
% numStepsCurve=ceil(sqrt(sum((p1-p2).^2))/pointSpacing);
% numStepsCurve=numStepsCurve+double(iseven(numStepsCurve));
% [Vg]=sweepCurveSmooth(p1,p2,n1,n2,numStepsCurve,pp,f);

pointSpacingHead=mean(patchEdgeLengths(F_head,V_head));
d=sqrt(sum((p1-p2).^2));
numStepsCurve=ceil(d/pointSpacingHead);
numStepsCurve=numStepsCurve+double(iseven(numStepsCurve));
f=d/3;
p=[p1;p1+f*n1; p2-f*n2;p2];

Vg=bezierCurve(p,numStepsCurve);

[F_loft,V_loft,C_loft]=sweepLoft(V_loft1,V_loft2,n1,n2,Vg);
 F_loft=fliplr(F_loft);
E=F_loft(~iseven(C_loft),[1 2]);
VE=patchCentre(E,V_loft);
V_loft(E(:,1),:)=VE;
indBoundaryCurve_head=1:numStepsCurve:size(V_loft,1);
indBoundaryCurve_stem=numStepsCurve:numStepsCurve:size(V_loft,1);
[F_loft,V_loft,~]=quad2tri(F_loft,V_loft,'a',C_loft);

%%

C_loft=(max(C_head)+1)*ones(size(F_loft,1),1);
C_stem=C_stem+max(C_loft);
[F_implant,V_implant,C_implant]=joinElementSets({F_head,F_loft,F_stem},{V_head,V_loft,V_stem},{C_head,C_loft,C_stem});
[F_implant,V_implant,~,indFix]=mergeVertices(F_implant,V_implant);
indBallHole=indFix(indBallHole);
indBoundaryCurve_head=indFix(indBoundaryCurve_head+size(V_head,1));
indBoundaryCurve_stem=indFix(indBoundaryCurve_stem+size(V_head,1));
indBoundaryCurveTip=indFix(indBoundaryCurveTip+size(V_head,1)++size(V_loft,1));

%%

curveSet{1}=indBallHole;
curveSet{2}=indBoundaryCurve_head;
curveSet{3}=indBoundaryCurve_stem;
curveSet{4}=indBoundaryCurveTip;

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
