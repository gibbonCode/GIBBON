function [A] = ellipseFit(varargin)

% function [A] = ellipseFit(V,optMethod,numSample)
% ------------------------------------------------------------------------
% This function aims to fit an ellips to the points specified in the vertex
% array V in a least squares sense. The returned vector A contains the
% center, radii, and orientation of the ellipse:
% A = [uCentre, vCentre, Ru, Rv, thetarad];
%
% If the number of points is 1 then A=[]
% If the number of points is >1 and <4 then a circle is fitted instead
% where the centre of the cirlce is simply the mean of the points and its
% radius is the mean radius with respect to that centre.
%
% If >3 points are provided two methods can be chosen the two are largely
% equivalent but handle "rounding errors" differently. Both methods are
% based on the publication: %
% Authors: Andrew Fitzgibbon, Maurizio Pilu, Bob Fisher
% Reference: "Direct Least Squares Fitting of Ellipses", IEEE T-PAMI, 1999
%
%  @Article{Fitzgibbon99,
%   author = "Fitzgibbon, A.~W.and Pilu, M. and Fisher, R.~B.",
%   title = "Direct least-squares fitting of ellipses",
%   journal = pami,
%   year = 1999,
%   volume = 21,
%   number = 5,
%   month = may,
%   pages = "476--480"
%  }
%
% If optMethod==1 then the fitting procedure is a minor adjustment of the
% code provided here:
%http://research.microsoft.com/en-us/um/people/awf/ellipse/fitellipse.html
%
% If optMethod==2 then the fitting procedure follows the code given here:
% http://www.mathworks.nl/matlabcentral/fileexchange/22684-ellipse-fit-direct-method
%
% The above holds if the final input numSample=[]. If numSample is not
% empty cubic interpolation is used to up/down sample the polygon to
% numSample points. This can be used to avoid circle fitting polygons
% sample with only few points (e.g. <=4).
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2013/24/09
% 2016/10/14 Added varargin handling
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        V=varargin{1};
        optMethod=2;
        numSample=[];
    case 2
        V=varargin{1};
        optMethod=varargin{2};
        numSample=[];
    case 3
        V=varargin{1};
        optMethod=varargin{2};
        numSample=varargin{3};
end

%%
numPoints=size(V,1);
if ~isempty(numSample) %Up/Down sample polygon to numSample points
    if numPoints>2
        [V]=evenlySampleCurve(V,numSample,'pchip',1);
        numPoints=size(V,1);
    else
        warning(['Skipped upsampling of points since only ', num2str(numPoints), ' points were provided!']);
    end
end

switch numPoints
    case 1
        warning(['Only ', num2str(numPoints), ' points provided!'])
        warning('Returning empty output');
        A=[];
    case {2,3,4}
        warning(['Only ', num2str(numPoints), ' points provided!'])
        warning('Fitted circle instead centered on mean with average radius');
        
        x=V(:,1); y=V(:,2);
        
        %Get centroid
        uCentre = mean(x);
        vCentre = mean(y);
        
        %Get mean radius
        x = (x-uCentre);
        y = (y-vCentre);
        R=mean(sqrt(x.^2+y.^2));
        Ru=R;
        Rv=R;
        
        %Arbitrary angle since its a circle
        thetarad=0;
        
        A = [uCentre, vCentre, Ru, Rv, thetarad];
        
        %     case 3 Fit circumcirlce to triangle
        %         warning('Only 3 points provided!')
        %         warning('Output defines the circumcircle of the triangel!');
        %
        %         x=V(:,1); y=V(:,2);
        %
        %         TR = triangulation([1 2 3],x,y);
        %         CC = circumcenter(TR);
        %
        %         %Get centroid
        %         uCentre = CC(1);
        %         vCentre = CC(2);
        %
        %         %Get mean radius
        %         x = (x-uCentre);
        %         y = (y-vCentre);
        %         R=mean(sqrt(x.^2+y.^2));
        %         Ru=R;
        %         Rv=R;
        %
        %         %Arbitrary angle since its a circle
        %         thetarad=0;
        %
        %         A = [uCentre, vCentre, Ru, Rv, thetarad];
        
    otherwise
        switch optMethod
            case 1
                [par]=ellipseFitMethod1(V);
            case 2
                [par]=ellipseFitMethod2(V);
        end
        
        % Convert to geometric radii, and centers
        thetarad = 0.5*atan2(par(2),par(1) - par(3));
        cost = cos(thetarad);
        sint = sin(thetarad);
        sin_squared = sint.*sint;
        cos_squared = cost.*cost;
        cos_sin = sint .* cost;
        
        Ao = par(6);
        Au =   par(4) .* cost + par(5) .* sint;
        Av = - par(4) .* sint + par(5) .* cost;
        Auu = par(1) .* cos_squared + par(3) .* sin_squared + par(2) .* cos_sin;
        Avv = par(1) .* sin_squared + par(3) .* cos_squared - par(2) .* cos_sin;
        
        tuCentre = - Au./(2.*Auu);
        tvCentre = - Av./(2.*Avv);
        wCentre = Ao - Auu.*tuCentre.*tuCentre - Avv.*tvCentre.*tvCentre;
        
        uCentre = tuCentre .* cost - tvCentre .* sint;
        vCentre = tuCentre .* sint + tvCentre .* cost;
        
        Ru = -wCentre./Auu;
        Rv = -wCentre./Avv;
        
        Ru = sqrt(abs(Ru)).*sign(Ru);
        Rv = sqrt(abs(Rv)).*sign(Rv);
        
        %Fixing angle to represent major axis
        if Ru<Rv %If the first radius is smaller than the second than add offset to angle
            
            Rminor=Ru;
            Rmajor=Rv;
            
            thetarad=thetarad+(pi/2);
            
            Ru=Rmajor;
            Rv=Rminor;
        end
        thetarad=rem(thetarad,pi);
        A = [uCentre, vCentre, Ru, Rv, thetarad];
        
end
end

function [par]=ellipseFitMethod1(V)

x=V(:,1); y=V(:,2);

% normalize data
mx = mean(x);
my = mean(y);
sx = (max(x)-min(x))/2;
sy = (max(y)-min(y))/2;

x = (x-mx)/sx;
y = (y-my)/sy;

% Force to column vectors
x = x(:);
y = y(:);

% Build design matrix
D = [ x.*x  x.*y  y.*y  x  y  ones(size(x)) ];

% Build scatter matrix
S = D'*D;

% Build 6x6 constraint matrix
C(6,6) = 0; C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;

% Solve eigensystem

% Old way, numerically unstable if not implemented in matlab
[gevec, geval] = eig(S,C);

% Find the negative eigenvalue
checkFactor=1e-8; %N.B. CHECK THIS
L = real(diag(geval)) < checkFactor & ~isinf(diag(geval));

% Extract eigenvector corresponding to negative eigenvalue
A = real(gevec(:,L));

% unnormalize
par = [
    A(1)*sy*sy,   ...
    A(2)*sx*sy,   ...
    A(3)*sx*sx,   ...
    -2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy,   ...
    -A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy,   ...
    A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my   ...
    - A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my   ...
    + A(6)*sx*sx*sy*sy   ...
    ]';
end

function [par]=ellipseFitMethod2(V)
centroid = mean(V);   % the centroid of the data set

D1 = [(V(:,1)-centroid(1)).^2, (V(:,1)-centroid(1)).*(V(:,2)-centroid(2)),...
    (V(:,2)-centroid(2)).^2];
D2 = [V(:,1)-centroid(1), V(:,2)-centroid(2), ones(size(V,1),1)];
S1 = D1'*D1;
S2 = D1'*D2;
S3 = D2'*D2;
T = -inv(S3)*S2';
M = S1 + S2*T;
M = [M(3,:)./2; -M(2,:); M(1,:)./2];
[evec,~] = eig(M);
cond = 4*evec(1,:).*evec(3,:)-evec(2,:).^2;
A1 = evec(:,cond>0);
A = [A1; T*A1];
A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
    A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
A(4) = A4;  A(5) = A5;  A(6) = A6;
A = A/norm(A);
par=A;
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
