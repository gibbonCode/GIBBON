function [V]=filletCurve(V,r,np,closedLoopOpt)

% [V]=filletCurve(V,r,np,closedLoopOpt)
% ------------------------------------------------------------------------
% This function fillets a curve based on the input radius r using np points
% per fillet arc. If closedLoopOpt==1 then closed end conditions are used
% such that the end and start regions are also filleted.
%
% %% EXAMPLE:
% Vt=[0 0 0; 10 0 0; 5 10 0; 10 0 10; 0 10 10; ];
% r=2; %Fillet radius
% np=25; %Number of points used to construct each fillet edge
% closedLoopOption=0; %Use 1 if curve represents a closed loop but containes unique points
% [VN]=filletCurve(Vt,r,np,closedLoopOption);
%
% figure; hold on;
% plotV(Vt,'k.-.','lineWidth',2,'MarkerSize',25);
% plotV(VN,'r.-','lineWidth',3);
% axis equal; view(3); axis tight;  grid on;
% drawnow;
% %%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/03/19
% 2017/06/03 Fixed bug in relation to angles being 0 or 180 degrees
%------------------------------------------------------------------------

%%

angleTolerance=1e-3;
%%

%Cope with 2D input
if size(V,2)==2
    V(:,3)=0;
    nDim=2;
else
    nDim=3;
end

numPoints=size(V,1);
if numPoints>2
    numSteps=numPoints-2;
    indStart=1;
    for q=1:1:numSteps
        [Vr]=filletEdgeSet(V(indStart:indStart+2,:),r,np,angleTolerance);
        V=[V(1:indStart,:);Vr;V(indStart+2:end,:)];
        indStart=indStart+size(Vr,1);
    end
else
    error('Input curve V has too few points!');
end

if closedLoopOpt==1
    V_closed=[V(end-1:end,:); V(1,:);];
    [Vr]=filletEdgeSet(V_closed,r,np,angleTolerance);
    V=[V(1:end-1,:); Vr];
    
    V_closed=[V(end,:); V(1:2,:);];
    [Vr]=filletEdgeSet(V_closed,r,np,angleTolerance);
    V=[Vr; V(2:end,:); ];
end

if nDim==2
    V=V(:,1:2);
end
end

function [Vr]=filletEdgeSet(V,r,np,angleTolerance)
Vc=V-V(2*ones(size(V,1),1),:);

P1=Vc(1,:);
M1=norm(P1);
N1=P1./M1;

P2=Vc(3,:);
M2=norm(P2);
N2=P2./M2;

M_min=min([M1 M2]);

rMax=sqrt(sum((M_min.*N1-M_min.*N2).^2))/2;

if r>rMax
    error(['Radius is too big! Current max: ',num2str(rMax)]);
end

Nm=(N1+N2)/2;
Nm=Nm./norm(Nm);

a=real(acos(dot(P1,Nm)./norm(P1)));

if mod(abs(a),pi)<angleTolerance
    Vr=[];
else
    d=r./sin(a);
    d1=r./tan(a);
    d2=r./tan(a);
    
    Pc=Nm.*d;
    
    P1c=N1.*d1;
    P2c=N2.*d2;
    
    P1cc=P1c-Pc;
    P2cc=P2c-Pc;
    
    Vn=linspacen(P1cc,P2cc,np)';
    [theta_Vn,phi_Vn,~] = cart2sph(Vn(:,1),Vn(:,2),Vn(:,3));
    [Vn(:,1),Vn(:,2),Vn(:,3)]=sph2cart(theta_Vn,phi_Vn,r.*ones(size(phi_Vn)));
    Vn=Vn+Pc(ones(1,size(Vn,1)),:);
    
    Vr=Vn+V(2*ones(size(Vn,1),1),:);
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
