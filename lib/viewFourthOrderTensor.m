function viewFourthOrderTensor(varargin)

% function viewFourthOrderTensor(C,numDigits,fontSizeIm,fontSize)


%% Parse input

numDigits=[];
fontSizeIm=[];
fontSize=[];
toleranceLevel=eps(1);
switch nargin
    case 1
        C=varargin{1};        
        toleranceLevel=eps(1);
    case 2
        C=varargin{1};
        numDigits=varargin{2};                
    case 3
        C=varargin{1};
        numDigits=varargin{2};
        fontSizeIm=varargin{3};        
    case 4
        C=varargin{1};
        numDigits=varargin{2};
        fontSizeIm=varargin{3};
        fontSize=varargin{4};        
    case 5        
        C=varargin{1};
        numDigits=varargin{2};
        fontSizeIm=varargin{3};
        fontSize=varargin{4};
        toleranceLevel=varargin{5};
end

if isempty(C)
    I=eye(3,3); %The 2nd order identity tensor
    C=dyadicProduct(I,I,1); %A type of 4th order identity tensor
end

if isempty(numDigits)
    numDigits=3; 
end

if isempty(fontSizeIm)
    fontSizeIm=15; 
end

if isempty(fontSize)
    fontSize=15; 
end

%%
% Plot settings
faceAlpha=0.5;

%%

ind_C_all=1:numel(C);
[I,J,K,L]=ind2sub(size(C),ind_C_all);

kronD_IJ=(I==J);
kronD_KL=(K==L);

switch class(C)
    case 'sym'
        CV=sym(zeros(6,6));
        logicSym=1;
    otherwise
        CV=zeros(6,6);
        logicSym=0;
end


%Create mapping indices
P=I.*kronD_IJ+(1-kronD_IJ).*(9-I-J);
Q=K.*kronD_KL+(1-kronD_KL).*(9-K-L);

%Convert to linear indices
[ind_CV]=sub2ind(size(CV),P,Q);

%Get the unique 36 entries
[~,indUni,~]=unique(ind_CV,'first');
ind_CV_uni=ind_CV(indUni);

ind_C_uni=ind_C_all(indUni);

CV(ind_CV_uni)=C(ind_C_uni);

%%

%Derive 9x9 matrix representation
[CM]=fourthOrderMat(C);

if logicSym
    symbolicVariableSet=symvar(CM);
    C_dummy=double(subs(CM,symbolicVariableSet,1:numel(symbolicVariableSet)));
    logicEntry=abs(C_dummy)>toleranceLevel;
    C_dummy_V=double(subs(CV,symbolicVariableSet,1:numel(symbolicVariableSet)));
    logicEntry_V=abs(C_dummy_V)>toleranceLevel;
else
    logicEntry=abs(CM)>toleranceLevel;
    logicEntry_V=abs(CV)>toleranceLevel;
end

%%

C_V=zeros(size(C));
C_V(ind_C_uni)=ind_CV_uni;
[CM_V]=fourthOrderMat(C_V);
logicEntry_CV=CM_V>eps(CM_V);

%%
hf1=cFigure;
title('9x9 array mapping of $4^{th}$ order tensor','fontSize',fontSize,'Interpreter','LATEX');
xlabel('$q=3(j-1)+l$','fontSize',fontSize,'Interpreter','LATEX');
ylabel('$p=3(i-1)+k$','fontSize',fontSize,'Interpreter','LATEX');

if logicSym
    [Fp,Vp,Cp]=ind2patch(logicEntry,C_dummy,'sk');    
    gpatch(Fp,Vp,Cp,0.25.*ones(1,3),faceAlpha);
else
    [Fp,Vp,Cp]=ind2patch(logicEntry,CM,'sk');
    gpatch(Fp,Vp,Cp,0.25.*ones(1,3),faceAlpha);
end
colormap gjet;
if ~logicSym
    colorbar;
end

m=zeros(3,3);
[Fp,Vp,~]=ind2patch(true(size(m)),m,'sk');
Vp(:,[1 2])=3*Vp(:,[1 2])-1;
patch('Faces',Fp,'Vertices',Vp,'FaceColor','none','EdgeColor','k','lineWidth',5);

[Fp,Vp,~]=ind2patch(logicEntry_CV,logicEntry_CV,'sk');
patch('Faces',Fp,'Vertices',Vp,'FaceColor',0.25.*ones(1,3),'EdgeColor','k','lineWidth',1,'FaceAlpha',0.25);

image_numeric(CM,hf1,numDigits,fontSizeIm,'k','tex');

view(2); axis tight; axis ij; axis square;
set(gca,'FontSize',fontSize);
drawnow;

%%

hf2=cFigure;
title('Voigt array mapping of $4^{th}$ order tensor','fontSize',fontSize,'Interpreter','LATEX');
xlabel('$q=k\delta_{kl}+(1-\delta_{kl})(9-k-l)$','fontSize',fontSize,'Interpreter','LATEX');
ylabel('$p=i\delta_{ij}+(1-\delta_{ij})(9-i-j)$','fontSize',fontSize,'Interpreter','LATEX');

if logicSym
    [Fp,Vp,Cp]=ind2patch(logicEntry_V,C_dummy_V,'sk');
    gpatch(Fp,Vp,Cp,0.25.*ones(1,3),faceAlpha);
else
    [Fp,Vp,Cp]=ind2patch(logicEntry_V,CV,'sk');
    gpatch(Fp,Vp,Cp,0.25.*ones(1,3),faceAlpha);
end
colormap gjet;
if ~logicSym
    colorbar;
end

m=zeros(2,2);
[Fp,Vp,~]=ind2patch(true(size(m)),m,'sk');
Vp(:,[1 2])=3*Vp(:,[1 2])-1;
patch('Faces',Fp,'Vertices',Vp,'FaceColor','none','EdgeColor','k','lineWidth',5);

image_numeric(CV,hf2,numDigits,fontSizeIm,'k','tex');
view(2); axis tight; axis ij; axis square;
set(gca,'FontSize',fontSize);
drawnow;

 
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
