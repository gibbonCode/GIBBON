function viewFourthOrderTensor2D(C,fontSizeIm)

%%
% Plot settings
fig_color='w'; 
fig_colordef='white'; 
faceAlpha=0.5;
fontSize=15; 

%%
Cv=voigtMap(C);

ind_C_all=1:numel(C);
[I,J,K,L]=ind2sub(size(C),ind_C_all);

kronD_IJ=(I==J);
kronD_KL=(K==L);

switch class(C)
    case 'sym'
        CV=sym(zeros(6,6));
    otherwise
        CV=zeros(6,6);
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

switch class(C)
    case 'sym'
        logicSym=1;   
        
        symbolicVariableSet=symvar(CM);
        C_dummy=double(subs(CM,symbolicVariableSet,1:numel(symbolicVariableSet)));
        logicEntry=C_dummy>eps(C_dummy);    
        
        C_dummy_V=double(subs(CV,symbolicVariableSet,1:numel(symbolicVariableSet)));
        logicEntry_V=C_dummy_V>eps(C_dummy_V);    
        
    otherwise
        logicSym=0;
        logicEntry=CM>eps(CM);        
        logicEntry_V=CV>eps(CV);        
end

%%

C_V=zeros(size(C));
C_V(ind_C_uni)=ind_CV_uni;
[CM_V]=fourthOrderMat(C_V);
logicEntry_CV=CM_V>eps(CM_V);
        
%%
figuremax(fig_color,fig_colordef);
title('9x9 array mapping of 4^{th} order tensor','fontSize',fontSize);
xlabel('q=3*(j-1)+l','fontSize',fontSize);ylabel('p=3*(i-1)+k','fontSize',fontSize); 

if logicSym
    [Fp,Vp,Cp]=ind2patch(logicEntry,C_dummy,'sk');    
    patch('Faces',Fp,'Vertices',Vp,'FaceColor','flat','CData',Cp,'EdgeColor','none','FaceAlpha',faceAlpha);
    colormap jet; 
else
    [Fp,Vp,Cp]=ind2patch(logicEntry,CM,'sk');
    patch('Faces',Fp,'Vertices',Vp,'FaceColor','flat','CData',Cp,'EdgeColor','none','FaceAlpha',faceAlpha);
    colormap jet; colorbar;
end

m=zeros(3,3);
[Fp,Vp,~]=ind2patch(true(size(m)),m,'sk');
Vp(:,[1 2])=3*Vp(:,[1 2])-1;
patch('Faces',Fp,'Vertices',Vp,'FaceColor','none','EdgeColor','k','lineWidth',5);

[Fp,Vp,~]=ind2patch(logicEntry_CV,logicEntry_CV,'sk');
patch('Faces',Fp,'Vertices',Vp,'FaceColor',0.25.*ones(1,3),'EdgeColor','k','lineWidth',1,'FaceAlpha',0.25);

image_numeric(CM,gcf,2,fontSizeIm);
% view(2); axis equal; axis tight; axis ij;
view(2); axis tight; axis ij;
set(gca,'FontSize',fontSize);
drawnow;
        
%%
figuremax(fig_color,fig_colordef);
title('Voigt array mapping of 4^{th} order tensor','fontSize',fontSize);
xlabel('q=k*\delta_{kl}+(1-\delta_{kl})*(9-k-l)','fontSize',fontSize);ylabel('p=i*\delta_{ij}+(1-\delta_{ij})*(9-i-j)','fontSize',fontSize); 

if logicSym
    [Fp,Vp,Cp]=ind2patch(logicEntry_V,C_dummy_V,'sk');
    patch('Faces',Fp,'Vertices',Vp,'FaceColor','flat','CData',Cp,'EdgeColor','none','FaceAlpha',faceAlpha);
    colormap jet; 
else
    [Fp,Vp,Cp]=ind2patch(logicEntry_V,CV,'sk');
    patch('Faces',Fp,'Vertices',Vp,'FaceColor','flat','CData',Cp,'EdgeColor','none','FaceAlpha',faceAlpha);
    colormap jet; colorbar;
end

m=zeros(2,2);
[Fp,Vp,~]=ind2patch(true(size(m)),m,'sk');
Vp(:,[1 2])=3*Vp(:,[1 2])-1;
patch('Faces',Fp,'Vertices',Vp,'FaceColor','none','EdgeColor','k','lineWidth',5);

image_numeric(CV,gcf,2,fontSizeIm);
% view(2); axis equal; axis tight; axis ij;
view(2); axis tight; axis ij;
set(gca,'FontSize',fontSize);
drawnow;

