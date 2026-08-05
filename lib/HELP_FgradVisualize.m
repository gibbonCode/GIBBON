%% FgradVisualize
% Below is a demonstration of the features of the |FgradVisualize| function

%%
clear; close all; clc;

%% Syntax
% FgradVisualize(Fbig,boxDim,F,V,plt,ang)

%% Description
    %Creates an interactive visualization of a deformation or a series of
    %deformations, on a single element geometry. Final output plot has an
    %animation slider, and can be exported as a gif for external use. 
    
    %----MANDATORY ARGUMENTS----%
    %Fbig: 1xn Cell array containing n deformation gradients to be plotted (3x3 matrices)
    
    %----OPTIONAL ARGUMENTS (leave as [] if not needed)----:
    %- Fibre angle (ang): A 1x3 or 3x1 matrix of directional cosines corresponding
    %to the element fibre direction. Include to have that visualized in the
    %plot
    
    %- Plot Handle (plt): A plot handle that the figure can optionally be
    %plotted on. Use to generate the interactive figure on an existing
    %figure/subplot. If left as [], interactive figure will be generated on a fresh
    %plot.
    
    %Geometry: Two options available:
    %  - If you have custom geometry, provide patch faces (F) and nodes
    %  (V).
    %  - If not, a generic rectangular element can be used. Provide
    %  dimensions of the box (boxDim) as a 1x3 vector [L, W, T].     
    %  
    
    %  If none of the above is provided, a single-element 1x1 cube will be
    %  used.
    
    %NOTE: If using F,V, leave boxDim = [], and if using boxDim leave F & V
    %as []. If both are provided, F,V will be used.

    
    %Made by Darshan Senthil, 2026
%% Examples
%% Example 1: Visualize uniaxial deformations on a cube
% Visualize uniaxial strain on a 1x1x1 cube in each direction

Fx = [2,0,0;0,1,0;0,0,1]; %Varying strain for each
Fy = [1,0,0;0,5.5,0;0,0,1]; 
Fz = [1,0,0;0,1,0;0,0,0.4]; 

%Plotting results
plt=cFigure;
subplot(1,3,1)
FgradVisualize({Fx},[],[],[],plt,[]) %Since we only require a unit cube, specifying geometry inputs is not required.
subplot(1,3,2)
FgradVisualize({Fy},[],[],[],plt,[])
subplot(1,3,3)
FgradVisualize({Fz},[],[],[],plt,[])

%% Example 2: Shearing a Tetrahedron
% Visualize a shear experiment on a single 4-noded tetrahedron element

%Setting up geometry
V = [0.0  0.0  0.0; 5.0  0.0  0.0;0.0  5.0  0.0;0.0  0.0  5.0];
E = [1  2  3  4];
[Ff,~]=element2patch(E,elementType='tet4');

ang = [0.65, 0.35, 0.45]; %Defining a fibre direction
beta = 0:0.01:0.2; %20 data points
plt=cFigure;

for x = 1:length(beta)
    F = eye(3);
    F(1,2)=beta(x); %Shear in x-y plane
    Fbig{x}=F;
end

FgradVisualize(Fbig,[],Ff,V,plt,ang) %Note that you can also track the principal stretch directions and volume change (J)

%% Example 3: Uniaxial Test on a rectangular aorta sample
% Visualize a uniaxial stress test on a sample of aorta


%Initialize Parameters of interest to calculate stress
k1 = 0.0368; %[MPa]
k2 = 0.25; %[-]
angle1 = 41; %[deg]
angle2 = -angle1; %[deg]
D1 = 0.0001; %[1/MPa]
C10 = 0.001; %[MPa]
params = [k1, k2, D1, C10, angle1];
xT=[1;1];

boxDim = [2,1,5]; %Sample Dimensions
lam = (1:0.001:1.8)'; %Go up to 50% strain, or 1.5x stretch
ang = [cosd(angle1), 0, sind(angle1)]; %Defining a fibre direction
plt=cFigure;


for x = 1:length(lam)
    lam3=lam(x);
    %Uniaxial stress - numerically solve for values of lam1, lam2 that make
    %S11 = S22 = 0.
    xT=fminsearch(@(xT)transStretchSolver(xT,lam3,params),[lam3^(-1/2.5); lam3^(-1/2.5)]);
    %Get full deformation gradient; solve for final stress
    F=diag([xT(1),xT(2),lam3]);
    P = stressCalc(F,params);
    %Save for plotting 
    Fbig{x}=F;
end

%Making plot of deformed shape
FgradVisualize(Fbig,boxDim,[],[],plt,ang); %Note that you can also track the principal stretch directions and volume change (J)

grid on; grid minor
%%
%======HELPER FUNCTIONS FOR EXAMPLE 3======%

%Helper function to calculate stress
function P = stressCalc(F,params)

% load circumferentialdata
% e = circumData(:,2);

%Initialize Parameters of interest
    k1 = params(1); %[MPa]
    k2 = params(2); %[-]
    angle1 = params(5); %[rad]
    angle2 = angle1; %[rad]
    D1 = params(3); %[MPa]
    C10 = params(4); %[1/MPa]
%Right & Left Cauchy Stress Tensors
    C = F'*F;
    B = F*F';
    J = det(F);
%First Invariant
    I1 = B(1,1)+B(2,2)+B(3,3);
%Isochoric Cauchy Stress
    Siso = (2*C10/J)*(B*(J^(-2/3))-I1*eye(3)*(3*J^(2/3))^(-1));
%Volumetric Cauchy Stress
    Svol = (2/D1)*(J-1)*eye(3);
%Angle vector
    m4 = [cosd(angle1);0;sind(angle1)];
    m6 = [cosd(angle2);0;sind(angle2)];
%Isochoric Aniso Invariants
    I4 = m4'*C*m4;
    I6 = m6'*C*m6;
%Anisotropic Cauchy Stress Tensor
    %I4 contributioon
    a4 = F*m4;
    a6 = F*m6;
    Si4 = 2*k1*(I4-1)*exp(k2*(I4-1)^2)*(a4*a4');
    %I6 contribution
    Si6 = 2*k1*(I6-1)*exp(k2*(I6-1)^2)*(a6*a6');
    %Sum them up
    Saniso = Si4 + Si6;
    %Final Cauchy Stress
    Sc = Siso+Svol+Saniso;
%Piola-Kirchoff Stress Tensor
    P = J*Sc*(inv(F))';
end

function residual = transStretchSolver(xT,x3,params)
%Calculates 11 and 22 stresses (residual) from 11 and 22 stretches.
    F = diag([xT(1),xT(2),x3]);
    P = stressCalc(F,params);
    residual = sqrt(P(1,1)^2 + P(2,2)^2);
end

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2022 Kevin Mattheus Moerman and the GIBBON contributors
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
