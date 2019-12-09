function [Fn,Vn,Cn]=element2lattice(varargin)

% function [Fn,Vn,Cn]=element2lattice(E,V,cPar)

% 2017/04/26 fixed bug for older MATLAB versions in relation to subtracting
% columns from matrices e.g. subtracting an nx1 column from all columns in
% an nxm matrix.
% 2018/05/07 Updated handling of default input (no change in functionality)
% 2018/05/07 Added default boundary face use if not provided at all. 

%% Parse input

switch nargin
    case 2
        E=varargin{1};
        V=varargin{2};
        cPar=[];
    case 3
        E=varargin{1};
        V=varargin{2};
        cPar=varargin{3};
end

%The default control parameters
cParDefault.shrinkFactor=0.25;
cParDefault.latticeSide=1;
cParDefault.numDigitKeep=5;
cParDefault.meshType='tri';
cParDefault.indBoundary=[];
cParDefault.hexSplit=0;
cParDefault.hexMethod=2;
cParDefault.elementType=[];

%Complement input structure with default
[cPar]=structComplete(cPar,cParDefault,0);

if isempty(cPar.elementType)
    switch size(E,2)
        case 4
            cPar.elementType='tet4';        
        case 8
            cPar.elementType='hex8';
        case 14 
            cPar.elementType='rhomdo14';
        otherwise
            error('Element type missing');
    end
end

elementType=cPar.elementType;
switch elementType
    case 'tet4'
    case 'hex8'
    case 'rhomdo14'        
    otherwise
        error('Element type missing or not supported');
end
[F]=element2patch(E,elementType); %Element faces

%If no boundary entry is provided 
if isempty(cPar.indBoundary)     
    cPar.indBoundary=tesBoundary(F,V);
end

%%

if cPar.latticeSide==1
    cPar.shrinkFactor=1-cPar.shrinkFactor;
end

if size(cPar.shrinkFactor,1)==size(V,1)
   spatVarMode=true(1);
else
    spatVarMode=false(1);
end
    
%Get shrunk face coordinates
Vc=zeros(size(F,1)*size(F,2),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q);
    if size(F,1)==1
        FX=X(F)';
    else
        FX=X(F);
    end
    
    if spatVarMode
        S=mean(cPar.shrinkFactor(F),2);
    else
        S=cPar.shrinkFactor;
    end
    
    FX_mean=mean(FX,2);
    FX_mean=FX_mean(:,ones(size(FX,2),1));
    FX=((FX-FX_mean).*S)+FX_mean;
    Vc(:,q)=FX(:);
end
Fc=reshape(1:size(Vc,1),size(F,1),size(F,2));

%Get shrunk element coordinates
Vcc=zeros(size(E,1)*size(E,2),size(V,2));
Vcc_Ecc_mean=zeros(size(E,1),3);
for q=1:1:size(V,2)
    X=V(:,q);
    if size(E,1)==1
        EX=X(E)';
    else
        EX=X(E);
    end
    if spatVarMode
        S=mean(cPar.shrinkFactor(E),2);
    else
        S=cPar.shrinkFactor;
    end
    EX_mean=mean(EX,2);
    Vcc_Ecc_mean(:,q)=EX_mean;
    EX_mean=EX_mean(:,ones(size(EX,2),1));
    EX=((EX-EX_mean).*S)+EX_mean;
    Vcc(:,q)=EX(:);
end
Ecc=reshape(1:size(Vcc,1),size(E,1),size(E,2));
[Fcc,~]=element2patch(Ecc,[]);

%Get shrunk edge coordinates
[Fe]=patchEdges(F,0);

Ve=zeros(size(Fe,1)*2,size(V,2));
for q=1:1:size(V,2)
    X=V(:,q);
    if size(Fe,1)==1
        EX=X(Fe)';
    else
        EX=X(Fe);
    end
    
    if spatVarMode
        S=mean(cPar.shrinkFactor(Fe),2);
    else
        S=cPar.shrinkFactor;
    end
    
    EX_mean=mean(EX,2);
    EX_mean=EX_mean(:,ones(size(EX,2),1));
    EX=((EX-EX_mean).*S)+EX_mean;
    Ve(:,q)=EX(:);
end
Fec=reshape(1:size(Ve,1),size(Fe,1),size(Fe,2));
Fec1=reshape(Fec(:,1),[size(F,2),size(Fec,1)/size(F,2)])';
Fec2=reshape(Fec(:,2),[size(F,2),size(Fec,1)/size(F,2)])';

switch cPar.latticeSide
    case 1
        if strcmp(cPar.meshType,'hex')
            switch elementType
                case 'tet4'
                    ind1=(1:size(E,1));
                    ind2=ind1+size(E,1);
                    ind3=ind2+size(E,1);
                    ind4=ind3+size(E,1);
                    
                    Fn=[...
                        F(ind1,2)+size(Vc,1)+size(Vcc,1) F(ind1,1)+size(Vc,1)+size(Vcc,1) Fc(ind1,1) Fc(ind1,2)    Fc(ind2,1) Fc(ind2,2)  Fcc(ind1,1)+size(Vc,1) Fcc(ind1,2)+size(Vc,1);...
                        Fc(ind3,2) Fc(ind3,1) Fcc(ind3,1)+size(Vc,1) Fcc(ind3,2)+size(Vc,1)   F(ind3,2)+size(Vc,1)+size(Vcc,1) F(ind3,1)+size(Vc,1)+size(Vcc,1) Fc(ind1,1) Fc(ind1,3);...
                        Fc(ind4,2) Fc(ind4,1) Fcc(ind4,1)+size(Vc,1) Fcc(ind4,2)+size(Vc,1)   F(ind4,2)+size(Vc,1)+size(Vcc,1) F(ind4,1)+size(Vc,1)+size(Vcc,1) Fc(ind1,3) Fc(ind1,2);...
                        F(ind2,1)+size(Vc,1)+size(Vcc,1) F(ind2,3)+size(Vc,1)+size(Vcc,1) Fc(ind2,3) Fc(ind2,1)    Fc(ind4,2) Fc(ind4,3)  Fcc(ind4,3)+size(Vc,1) Fcc(ind4,2)+size(Vc,1);...
                        F(ind2,2)+size(Vc,1)+size(Vcc,1) F(ind2,3)+size(Vc,1)+size(Vcc,1) Fc(ind3,3) Fc(ind3,1) Fc(ind2,2) Fc(ind2,3) Fcc(ind2,3)+size(Vc,1) Fcc(ind2,2)+size(Vc,1);...
                        Fc(ind4,1) Fc(ind4,3) Fcc(ind4,3)+size(Vc,1) Fcc(ind4,1)+size(Vc,1)   F(ind4,1)+size(Vc,1)+size(Vcc,1) F(ind4,3)+size(Vc,1)+size(Vcc,1) Fc(ind3,3) Fc(ind3,2);...
                        ];
                    
                    Vn=[Vc;Vcc;V];
                    
                    if cPar.hexSplit>0
                        [Fn,Vn]=subHex(Fn,Vn,cPar.hexSplit,3);
                    end
                    Cn=zeros(size(Fn,1),1);
                case 'hex8'
                    
                    ind1=1:size(E,1);
                    ind2=ind1+size(E,1);
                    ind3=ind2+size(E,1);
                    ind4=ind3+size(E,1);
                    ind5=ind4+size(E,1);
                    ind6=ind5+size(E,1);
                    
                    offSet_e=size(Vc,1)+size(Vcc,1)+size(V,1);
                    offSet_cc=size(Vc,1);
                    offSet_F=size(Vc,1)+size(Vcc,1);
                    Fcc=Fcc+offSet_cc;
                    F=F+offSet_F;
                    Fec1=Fec1+offSet_e;
                    Fec2=Fec2+offSet_e;
               
                    switch cPar.hexMethod
                        case 1
                            Fn=[...
                                Fc(ind4,2) Fc(ind4,1) Fcc(ind1,2) Fcc(ind1,1)   F(ind1,1) F(ind1,2) Fc(ind1,2) Fc(ind1,1) ;... %1
                                F(ind1,1) F(ind1,4) Fc(ind1,4) Fc(ind1,1)       Fc(ind6,3) Fc(ind6,4) Fcc(ind1,4) Fcc(ind1,1);... %2
                                F(ind1,4) F(ind1,3) Fc(ind1,3) Fc(ind1,4)       Fc(ind3,1) Fc(ind3,2) Fcc(ind1,3) Fcc(ind1,4) ;... %3
                                Fc(ind5,2) Fc(ind5,1) Fcc(ind1,3) Fcc(ind1,2)   F(ind1,2) F(ind1,3) Fc(ind1,3) Fc(ind1,2);...%4
                                ...
                                Fc(ind4,4) Fc(ind4,3) Fcc(ind2,4) Fcc(ind2,3)   F(ind2,3) F(ind2,4) Fc(ind2,4) Fc(ind2,3);...%5
                                F(ind2,1) F(ind2,4) Fc(ind2,4) Fc(ind2,1)       Fc(ind6,1) Fc(ind6,2) Fcc(ind2,4) Fcc(ind2,1);...%6
                                Fc(ind3,4) Fc(ind3,3) Fcc(ind2,2) Fcc(ind2,1)   F(ind2,1) F(ind2,2) Fc(ind2,2) Fc(ind2,1);...%7
                                Fc(ind5,4) Fc(ind5,3) Fcc(ind2,3) Fcc(ind2,2)   F(ind2,2) F(ind2,3) Fc(ind2,3) Fc(ind2,2);...%8
                                ...
                                Fc(ind6,1) Fc(ind6,4) Fcc(ind3,1) Fcc(ind3,4)   F(ind3,4) F(ind3,1) Fc(ind3,1) Fc(ind3,4);...
                                F(ind3,3) F(ind3,2) Fc(ind3,2) Fc(ind3,3)       Fc(ind5,4) Fc(ind5,1) Fcc(ind3,2) Fcc(ind3,3);...
                                F(ind5,3) F(ind5,2) Fc(ind5,2) Fc(ind5,3)       Fc(ind4,4) Fc(ind4,1) Fcc(ind5,2) Fcc(ind5,3); ...
                                Fc(ind4,3) Fc(ind4,2) Fcc(ind6,3) Fcc(ind6,2)   F(ind6,2) F(ind6,3) Fc(ind6,3) Fc(ind6,2); ...
                                ];
                            Vn=[Vc;Vcc;V;];
                        
                            if cPar.hexSplit>0                                
                                [Fn,Vn]=subHex(Fn,Vn,cPar.hexSplit,3);                                
                            end
                            Cn=zeros(size(Fn,1),1);
                        case 2             
                            Vn=[Vc;Vcc;V;Ve];
                            El=[Fc(ind4,2) Fc(ind4,1) Fcc(ind1,2) Fcc(ind1,1)    Fec1(ind1,1) Fec2(ind1,1) Fc(ind1,2) Fc(ind1,1) ;... %1
                                Fec2(ind1,4) Fec1(ind1,4) Fc(ind1,4) Fc(ind1,1)  Fc(ind6,3) Fc(ind6,4) Fcc(ind1,4) Fcc(ind1,1);... %2
                                Fec2(ind1,3) Fec1(ind1,3) Fc(ind1,3) Fc(ind1,4)  Fc(ind3,1) Fc(ind3,2) Fcc(ind1,3) Fcc(ind1,4) ;... %3
                                Fc(ind5,2) Fc(ind5,1) Fcc(ind1,3) Fcc(ind1,2)    Fec1(ind1,2) Fec2(ind1,2) Fc(ind1,3) Fc(ind1,2);...%4
                                ...
                                Fc(ind4,4) Fc(ind4,3) Fcc(ind2,4) Fcc(ind2,3)    Fec1(ind2,3) Fec2(ind2,3) Fc(ind2,4) Fc(ind2,3);...%5
                                Fec2(ind2,4) Fec1(ind2,4) Fc(ind2,4) Fc(ind2,1)  Fc(ind6,1) Fc(ind6,2) Fcc(ind2,4) Fcc(ind2,1);...%6
                                Fc(ind3,4) Fc(ind3,3) Fcc(ind2,2) Fcc(ind2,1)    Fec1(ind2,1) Fec2(ind2,1) Fc(ind2,2) Fc(ind2,1);...%7
                                Fc(ind5,4) Fc(ind5,3) Fcc(ind2,3) Fcc(ind2,2)    Fec1(ind2,2) Fec2(ind2,2) Fc(ind2,3) Fc(ind2,2);...%8
                                ...
                                Fc(ind6,1) Fc(ind6,4) Fcc(ind3,1) Fcc(ind3,4)    Fec1(ind3,4) Fec2(ind3,4) Fc(ind3,1) Fc(ind3,4);...%9
                                Fec2(ind3,2) Fec1(ind3,2) Fc(ind3,2) Fc(ind3,3)  Fc(ind5,4) Fc(ind5,1) Fcc(ind3,2) Fcc(ind3,3);...%10
                                Fec2(ind5,2) Fec1(ind5,2) Fc(ind5,2) Fc(ind5,3)  Fc(ind4,4) Fc(ind4,1) Fcc(ind5,2) Fcc(ind5,3); ...%11
                                Fc(ind4,3) Fc(ind4,2) Fcc(ind6,3) Fcc(ind6,2)    Fec1(ind6,2) Fec2(ind6,2) Fc(ind6,3) Fc(ind6,2); ...%12
                                ];
                            
                            Es=[Fec2(ind6,2) Fc(ind4,2) Fcc(ind1,1) Fc(ind6,3)   F(ind1,1) Fec1(ind1,1) Fc(ind1,1) Fec2(ind1,4);...
                                Fec1(ind6,4) Fc(ind6,4) Fcc(ind1,4) Fc(ind3,1)   F(ind1,4) Fec1(ind1,4) Fc(ind1,4) Fec2(ind1,3);...
                                Fc(ind3,2)  Fcc(ind1,3)  Fc(ind5,1)  Fec2(ind5,4) Fec1(ind1,3) Fc(ind1,3) Fec2(ind1,2) F(ind1,3);...
                                Fc(ind4,1) Fec1(ind5,2) Fc(ind5,2)  Fcc(ind1,2)  Fec1(ind4,1)  F(ind1,2) Fec1(ind1,2) Fc(ind1,2)
                                ...
                                Fec1(ind3,4) Fc(ind3,4) Fcc(ind2,1) Fc(ind6,1) F(ind2,1)  Fec1(ind2,1) Fc(ind2,1) Fec2(ind2,4);...
                                Fec2(ind3,2) Fc(ind5,4) Fcc(ind2,2) Fc(ind3,3) F(ind2,2)  Fec1(ind2,2) Fc(ind2,2) Fec2(ind2,1);...
                                Fec2(ind5,2) Fc(ind4,4) Fcc(ind2,3) Fc(ind5,3) F(ind2,3)  Fec1(ind2,3) Fc(ind2,3) Fec2(ind2,2);...
                                Fec2(ind4,2) Fc(ind6,2) Fcc(ind2,4) Fc(ind4,3) F(ind2,4)  Fec1(ind2,4) Fc(ind2,4) Fec2(ind2,3);...
                                ];
                            
                            if cPar.hexSplit>0                                                               
                                [El,Vl]=patchCleanUnused(El,Vn);
                                [Es,Vs]=patchCleanUnused(Es,Vn);
                                [El,Vl]=subHex(El,Vl,cPar.hexSplit,3);
                                Vn=[Vl;Vs];
                                Es=Es+size(Vl,1);
                            end
                            Fn=[El; Es];                            
                            Cn=[zeros(size(El,1),1);ones(size(Es,1),1)];
                    end  
                case 'rhomdo14'
                                        
                    Vn=[Vc;Vcc;V;];
                    Fcc=Fcc+size(Vc,1);
                    F=F+size(Vc,1)+size(Vcc,1);                    
                    
                    ind1 = 1:size(E,1);
                    ind2 = ind1+size(E,1);
                    ind3 = ind2+size(E,1);
                    ind4 = ind3+size(E,1);
                    ind5 = ind4+size(E,1);
                    ind6 = ind5+size(E,1);
                    ind7 = ind6+size(E,1);
                    ind8 = ind7+size(E,1);
                    ind9 = ind8+size(E,1);
                    ind10= ind9+size(E,1);
                    ind11= ind10+size(E,1);
                    ind12= ind11+size(E,1);
                    
                    Fn=[... 
                        F(ind5,3)  Fc(ind5,3)  Fcc(ind5,3)  Fc(ind2,3)  F(ind5,2)  Fc(ind5,2)  Fcc(ind5,2)  Fc(ind2,4) ;... %1
                        F(ind1,3)  Fc(ind1,3)  Fcc(ind1,3)  Fc(ind5,1)  F(ind1,2)  Fc(ind1,2)  Fcc(ind1,2)  Fc(ind5,2);... %2
                        F(ind5,1)  Fc(ind5,1)  Fcc(ind5,1)  Fc(ind8,3)  F(ind5,4)  Fc(ind5,4)  Fcc(ind5,4)  Fc(ind8,4);... %3
                        F(ind5,4)  Fc(ind5,4)  Fcc(ind5,4)  Fc(ind6,4)  F(ind5,3)  Fc(ind5,3)  Fcc(ind5,3)  Fc(ind6,1);... %4
                        F(ind1,4)  Fc(ind1,4)  Fcc(ind1,4)  Fc(ind8,2)  F(ind1,3)  Fc(ind1,3)  Fcc(ind1,3)  Fc(ind8,3);... %5
                        F(ind1,1)  Fc(ind1,1)  Fcc(ind1,1)  Fc(ind12,3) F(ind1,4)  Fc(ind1,4)  Fcc(ind1,4)  Fc(ind12,4);... %6
                        F(ind1,2)  Fc(ind1,2)  Fcc(ind1,2)  Fc(ind9,4)  F(ind1,1)  Fc(ind1,1)  Fcc(ind1,1)  Fc(ind9,1);...%7
                        F(ind8,2)  Fc(ind8,2)  Fcc(ind8,2)  Fc(ind4,2)  F(ind8,1)  Fc(ind8,1)  Fcc(ind8,1)  Fc(ind4,3) ;... %8
                        F(ind8,1)  Fc(ind8,1)  Fcc(ind8,1)  Fc(ind7,3)  F(ind8,4)  Fc(ind8,4)  Fcc(ind8,4)  Fc(ind7,4);... %9
                        F(ind9,4)  Fc(ind9,4)  Fcc(ind9,4)  Fc(ind2,4)  F(ind9,3)  Fc(ind9,3)  Fcc(ind9,3)  Fc(ind2,1)  ;... %10
                        F(ind2,2)  Fc(ind2,2)  Fcc(ind2,2)  Fc(ind10,4) F(ind2,1)  Fc(ind2,1)  Fcc(ind2,1)  Fc(ind10,1);... %11
                        F(ind2,3)  Fc(ind2,3)  Fcc(ind2,3)  Fc(ind6,1)  F(ind2,2)  Fc(ind2,2)  Fcc(ind2,2)  Fc(ind6,2) ;... %12
                        F(ind6,3)  Fc(ind6,3)  Fcc(ind6,3)  Fc(ind3,3)  F(ind6,2)  Fc(ind6,2)  Fcc(ind6,2)  Fc(ind3,4);... %13
                        F(ind6,4)  Fc(ind6,4)  Fcc(ind6,4)  Fc(ind7,4)  F(ind6,3)  Fc(ind6,3)  Fcc(ind6,3)  Fc(ind7,1);... %14
                        F(ind3,3)  Fc(ind3,3)  Fcc(ind3,3)  Fc(ind7,1)  F(ind3,2)  Fc(ind3,2)  Fcc(ind3,2)  Fc(ind7,2) ;... %15
                        F(ind4,4)  Fc(ind4,4)  Fcc(ind4,4)  Fc(ind7,2)  F(ind4,3)  Fc(ind4,3)  Fcc(ind4,3)  Fc(ind7,3);... %16
                        F(ind3,1)  Fc(ind3,1)  Fcc(ind3,1)  Fc(ind10,3) F(ind3,4)  Fc(ind3,4)  Fcc(ind3,4)  Fc(ind10,4) ;... %17
                        F(ind3,2)  Fc(ind3,2)  Fcc(ind3,2)  Fc(ind11,4) F(ind3,1)  Fc(ind3,1)  Fcc(ind3,1)  Fc(ind11,1) ;... %18
                        F(ind9,3)  Fc(ind9,3)  Fcc(ind9,3)  Fc(ind10,1) F(ind9,2)  Fc(ind9,2)  Fcc(ind9,2)  Fc(ind10,2);... %19
                        F(ind9,2)  Fc(ind9,2)  Fcc(ind9,2)  Fc(ind12,2) F(ind9,1)  Fc(ind9,1)  Fcc(ind9,1)  Fc(ind12,3);... %20
                        F(ind10,3) Fc(ind10,3) Fcc(ind10,3) Fc(ind11,1) F(ind10,2) Fc(ind10,2) Fcc(ind10,2) Fc(ind11,2);... %21
                        F(ind11,3) Fc(ind11,3) Fcc(ind11,3) Fc(ind12,1) F(ind11,2) Fc(ind11,2) Fcc(ind11,2) Fc(ind12,2);... %22
                        F(ind4,1)  Fc(ind4,1)  Fcc(ind4,1)  Fc(ind11,3) F(ind4,4)  Fc(ind4,4)  Fcc(ind4,4)  Fc(ind11,4);... %23
                        F(ind4,2)  Fc(ind4,2)  Fcc(ind4,2)  Fc(ind12,4) F(ind4,1)  Fc(ind4,1)  Fcc(ind4,1)  Fc(ind12,1) ;... %24
                        ];
                    
                    if cPar.hexSplit>0
                        [Fn,Vn]=subHex(Fn,Vn,cPar.hexSplit,2);
                    end                    
                    Cn=zeros(size(Fn,1),1);
                    
            end
        else
            switch elementType
                case 'tet4'
                    Fn=[Fc(:,1) Fcc(:,1)+size(Vc,1) Fcc(:,2)+size(Vc,1) Fc(:,2);...
                        Fc(:,2) Fcc(:,2)+size(Vc,1) Fcc(:,3)+size(Vc,1) Fc(:,3);...
                        Fc(:,3) Fcc(:,3)+size(Vc,1) Fcc(:,1)+size(Vc,1) Fc(:,1)];
                    
                    Cn=zeros(size(Fn,1),1);
                    
                    %Boundary parts
                    if ~isempty(cPar.indBoundary)
                        Fb=[F(cPar.indBoundary,1)+size(Vc,1)+size(Vcc,1) Fc(cPar.indBoundary,1) Fc(cPar.indBoundary,2) F(cPar.indBoundary,2)+size(Vc,1)+size(Vcc,1);...
                            F(cPar.indBoundary,2)+size(Vc,1)+size(Vcc,1) Fc(cPar.indBoundary,2) Fc(cPar.indBoundary,3) F(cPar.indBoundary,3)+size(Vc,1)+size(Vcc,1);...
                            F(cPar.indBoundary,3)+size(Vc,1)+size(Vcc,1) Fc(cPar.indBoundary,3) Fc(cPar.indBoundary,1) F(cPar.indBoundary,1)+size(Vc,1)+size(Vcc,1)];
                        Fn=[Fn;Fb];
                        Cn=[Cn;ones(size(Fb,1),1)];
                    end
                    Fn=fliplr(Fn);
                case 'hex8'
                    Fn=[Fc(:,1) Fcc(:,1)+size(Vc,1) Fcc(:,2)+size(Vc,1) Fc(:,2);...
                        Fc(:,2) Fcc(:,2)+size(Vc,1) Fcc(:,3)+size(Vc,1) Fc(:,3);...
                        Fc(:,3) Fcc(:,3)+size(Vc,1) Fcc(:,4)+size(Vc,1) Fc(:,4);...
                        Fc(:,4) Fcc(:,4)+size(Vc,1) Fcc(:,1)+size(Vc,1) Fc(:,1)];
                    
                    Cn=zeros(size(Fn,1),1);
                    
                    %Boundary parts
                    if ~isempty(cPar.indBoundary)
                        Fb=[F(cPar.indBoundary,1)+size(Vc,1)+size(Vcc,1) Fc(cPar.indBoundary,1) Fc(cPar.indBoundary,2) F(cPar.indBoundary,2)+size(Vc,1)+size(Vcc,1);...
                            F(cPar.indBoundary,2)+size(Vc,1)+size(Vcc,1) Fc(cPar.indBoundary,2) Fc(cPar.indBoundary,3) F(cPar.indBoundary,3)+size(Vc,1)+size(Vcc,1);...
                            F(cPar.indBoundary,3)+size(Vc,1)+size(Vcc,1) Fc(cPar.indBoundary,3) Fc(cPar.indBoundary,4) F(cPar.indBoundary,4)+size(Vc,1)+size(Vcc,1);...
                            F(cPar.indBoundary,4)+size(Vc,1)+size(Vcc,1) Fc(cPar.indBoundary,4) Fc(cPar.indBoundary,1) F(cPar.indBoundary,1)+size(Vc,1)+size(Vcc,1);...
                            ];
                        Fn=[Fn;Fb];
                        Cn=[Cn;ones(size(Fb,1),1)];
                    end
                    Fn=fliplr(Fn);
            end
            Vn=[Vc;Vcc;V];
            switch cPar.meshType
                case 'tri'
                    Fn=[Fn(:,[1 2 3]); Fn(:,[3 4 1])];
                    Cn=[Cn;Cn];
            end
        end
        
    case 2
        %%
        if strcmp(cPar.meshType,'hex')
            switch elementType
                case 'tet4'                    
                    switch cPar.hexMethod
                        case 1
                            [Fcq,Vcq]=tri2quad(Fc,Vc);
                            [Fccq,Vccq]=tri2quad(Fcc,Vcc);
                            Vccq(end-size(Fcc)+1:end,:)=repmat(Vcc_Ecc_mean,[4 1]);
                            Fn=[Fcq Fccq+size(Vcq,1)];
                            Fn=Fn(:,[1 5 6 2 4 8 7 3]);    
                            Fn=Fn(:,[5:8 1:4]);
                            Vn=[Vcq;Vccq];
                            if cPar.hexSplit>0
                                [Fn,Vn]=subHex(Fn,Vn,cPar.hexSplit,3);
                            end
                            Cn=zeros(size(Fn,1),1);
                        case 2
                            [Fcq,Vcq]=tri2quad(Fc,Vc);
                            [Fccq,Vccq]=tri2quad(Fcc,Vcc);                            
                            [Es,Vs]=tet2hex(Ecc,Vcc);
                            Vn=[Vcq;Vccq;Vs];
                            Es=Es+size(Vcq,1)+size(Vccq,1);
                            El=[Fccq+size(Vcq,1) Fcq];
                            El=El(:,[1 5 6 2 4 8 7 3]);
                            if cPar.hexSplit>0
                                [El,Vl]=patchCleanUnused(El,Vn);
                                [Es,Vs]=patchCleanUnused(Es,Vn);
                                [El,Vl]=subHex(El,Vl,cPar.hexSplit,3);
                                Vn=[Vl;Vs];
                                Es=Es+size(Vl,1);
                            end                            
                            Fn=[El; Es];                            
                            Cn=[zeros(size(El,1),1);ones(size(Es,1),1)];
                    end
                case 'hex8'                    
                    switch cPar.hexMethod
                        case 1
                            [Fcq,Vcq]=subQuad(Fc,Vc,1);
                            [Fccq,Vccq]=subQuad(Fcc,Vcc,1);
                            Vccq(end-size(Fcc)+1:end,:)=repmat(Vcc_Ecc_mean,[6 1]);
                            Fn=[Fcq Fccq+size(Vcq,1)];
                            Fn=Fn(:,[1 5 6 2 4 8 7 3]);
                            Fn=Fn(:,[5:8 1:4]);
                            Vn=[Vcq;Vccq];                            
                            if cPar.hexSplit>0
                                [Fn,Vn]=subHex(Fn,Vn,cPar.hexSplit,3);
                            end                            
                            Cn=zeros(size(Fn,1),1);                            
                        case 2
                            Vn=[Vc;Vcc;];
                            El=[Fcc+size(Vc,1) Fc];
                            El=El(:,[1 5 6 2 4 8 7 3]);                            
                            Es=Ecc+size(Vc,1);                            
                            if cPar.hexSplit>0
                                [El,Vl]=patchCleanUnused(El,Vn);
                                [Es,Vs]=patchCleanUnused(Es,Vn);
                                [El,Vl]=subHex(El,Vl,cPar.hexSplit,3);
                                Vn=[Vl;Vs];
                                Es=Es+size(Vl,1);
                            end
                            Fn=[El; Es];
                            Cn=[zeros(size(El,1),1);ones(size(Es,1),1)];
                    end
                case 'rhomdo14'
                    V_ce=patchCentre(E,V);
                    Vn=[Vc;Vcc;V_ce];
                    Fcc=Fcc+size(Vc,1);                    
                    ind_ce=(1:1:size(V_ce,1))'+size(Vc,1)+size(Vcc,1);
                    
                    ind1 = 1:size(E,1);
                    ind2 = ind1+size(E,1);
                    ind3 = ind2+size(E,1);
                    ind4 = ind3+size(E,1);
                    ind5 = ind4+size(E,1);
                    ind6 = ind5+size(E,1);
                    ind7 = ind6+size(E,1);
                    ind8 = ind7+size(E,1);
                    ind9 = ind8+size(E,1);
                    ind10= ind9+size(E,1);
                    ind11= ind10+size(E,1);
                    ind12= ind11+size(E,1);
                    
                    
                    E1=[fliplr(Fcc(ind1,:)) fliplr([ind_ce Fcc(ind5,3) Fcc(ind5,4) Fcc(ind8,1)]) ];
                    E2=[fliplr(Fcc(ind2,:))  fliplr([Fcc(ind9,2) Fcc(ind10,3) ind_ce Fcc(ind9,1)]) ];
                    E3=[fliplr(Fcc(ind3,:)) fliplr([ind_ce  Fcc(ind7,3) Fcc(ind6,4) Fcc(ind6,1)]) ];
                    E4=[fliplr((Fcc(ind4,:))) fliplr([Fcc(ind12,2) Fcc(ind12,3) ind_ce Fcc(ind11,1)])  ];
                    
                    Fn=[...
                        Fcc(ind1,:)  Fc(ind1,:);... %1
                        Fcc(ind2,:)  Fc(ind2,:);... %2
                        Fcc(ind3,:)  Fc(ind3,:);... %3
                        Fcc(ind4,:)  Fc(ind4,:);... %4
                        Fcc(ind5,:)  Fc(ind5,:);... %5
                        Fcc(ind6,:)  Fc(ind6,:);... %6
                        Fcc(ind7,:)  Fc(ind7,:);... %7
                        Fcc(ind8,:)  Fc(ind8,:);... %8
                        Fcc(ind9,:)  Fc(ind9,:);... %9
                        Fcc(ind10,:) Fc(ind10,:);... %10
                        Fcc(ind11,:) Fc(ind11,:);... %11
                        Fcc(ind12,:) Fc(ind12,:);... %12
                        ];
                    
                    if cPar.hexSplit>0
                        [Fn,Vn]=subHex(Fn,Vn,cPar.hexSplit,2);
                    end
                    Cn=[zeros(size(Fn,1),1);ones(4*size(E1,1),1)];                    
                    Fn=[Fn;E1;E2;E3;E4];
                    
            end
        else
            switch elementType
                case 'tet4'
                    switch cPar.meshType
                        case 'tri'
                            Fn=[Fcc(:,1)+size(Vc,1) Fcc(:,1)+size(Vc,1) Fcc(:,2)+size(Vc,1) Fcc(:,2)+size(Vc,1);...
                                Fcc(:,2)+size(Vc,1) Fcc(:,2)+size(Vc,1) Fcc(:,3)+size(Vc,1) Fcc(:,3)+size(Vc,1);...
                                Fcc(:,3)+size(Vc,1) Fcc(:,3)+size(Vc,1) Fcc(:,1)+size(Vc,1) Fcc(:,1)+size(Vc,1)];
                            
                            Cn=zeros(size(Fn,1),1);
                            Vn=[Vc;Vcc;];
                            
                            %Tri
                            Fn=[Fn(:,[1 2 3]); Fn(:,[3 4 1])];
                            Cn=[Cn;Cn];
                            
                            %Boundary parts
                            if ~isempty(cPar.indBoundary)
                                Fb=Fc(cPar.indBoundary,:);
                                Fn=[Fn;Fb];
                                Cn=[Cn;ones(size(Fb,1),1)];
                            end
                            
                        case 'quad'
                            Fn=[Fc(:,1) Fcc(:,1)+size(Vc,1) Fcc(:,2)+size(Vc,1) Fc(:,2);...
                                Fc(:,2) Fcc(:,2)+size(Vc,1) Fcc(:,3)+size(Vc,1) Fc(:,3);...
                                Fc(:,3) Fcc(:,3)+size(Vc,1) Fcc(:,1)+size(Vc,1) Fc(:,1)];
                            Cn=zeros(size(Fn,1),1);
                            Vn=[Vc;Vcc;];
                            
                            %Boundary parts
                            if ~isempty(cPar.indBoundary)
                                %Boundary parts are triangles and require conversion to
                                %quad
                                [Fq,Vq]=tri2quad(Fc(cPar.indBoundary,:),Vc);
                                
                                %subdevide existing quads to match up with quads on
                                %boundary
                                [Fn,Vn]=subQuad(Fn,Vn,1);
                                Cn=zeros(size(Fn,1),1);
                                Cn=[Cn;ones(size(Fq,1),1)];
                                Fn=[Fn;Fq+size(Vn,1)];
                                Vn=[Vn;Vq];
                            end
                    end
                    %                 Fn=fliplr(Fn);
                case 'hex8'
                    
                    %Quad
                    Fn=[Fc(:,1) Fcc(:,1)+size(Vc,1) Fcc(:,2)+size(Vc,1) Fc(:,2);...
                        Fc(:,2) Fcc(:,2)+size(Vc,1) Fcc(:,3)+size(Vc,1) Fc(:,3);...
                        Fc(:,3) Fcc(:,3)+size(Vc,1) Fcc(:,4)+size(Vc,1) Fc(:,4);...
                        Fc(:,4) Fcc(:,4)+size(Vc,1) Fcc(:,1)+size(Vc,1) Fc(:,1);];
                    Cn=zeros(size(Fn,1),1);
                    
                    %Boundary parts
                    if ~isempty(cPar.indBoundary)
                        Fb=Fc(cPar.indBoundary,:);
                        Fn=[Fn;Fb];
                        Cn=[Cn;ones(size(Fb,1),1)];
                    end
                    Vn=[Vc;Vcc;];
                    
                    switch cPar.meshType
                        case 'tri'
                            Fn=[Fn(:,[1 2 3]); Fn(:,[3 4 1])];
                            Cn=[Cn;Cn];
                    end
            end
        end
end

%Removing double vertices
[Fn,Vn]=mergeVertices(Fn,Vn);
 
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
