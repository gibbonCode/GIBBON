function [Fn,Vn,Cn]=element2lattice(varargin)

% function [Fn,Vn,Cn]=element2lattice(E,V,cPar)

%% Parse input

switch nargin
    case 2
        E=varargin{1};
        V=varargin{2};
        cpar=[];
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

if isempty(cPar)
    cPar=cParDefault;
else
    if ~isfield(cPar,'shrinkFactor')
        cPar.shrinkFactor=cParDefault.shrinkFactor;
    end
    if ~isfield(cPar,'latticeSide')
        cPar.latticeSide=cParDefault.latticeSide;
    end
    if ~isfield(cPar,'numDigitKeep')
        cPar.numDigitKeep=cParDefault.numDigitKeep;
    end
    if ~isfield(cPar,'meshType')
        cPar.meshType=cParDefault.meshType;
    end
end

%%

switch size(E,2)
    case 4
        elementType='tet4';
    case 8
        elementType='hex8';
    otherwise
        error('Element type not supported');
end

[F]=element2patch(E); %Element faces

if cPar.latticeSide==1
    cPar.shrinkFactor=1-cPar.shrinkFactor;
end

Vc=zeros(size(F,1)*size(F,2),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q);
    if size(F,1)==1
        FX=X(F)';
    else
        FX=X(F);
    end
    FX_mean=mean(FX,2);
    FX=((FX-FX_mean)*cPar.shrinkFactor)+FX_mean;
    Vc(:,q)=FX(:);
end
Fc=reshape(1:size(Vc,1),size(F,1),size(F,2));

Vcc=zeros(size(E,1)*size(E,2),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q);
    if size(E,1)==1
        EX=X(E)';
    else
        EX=X(E);
    end
    EX_mean=mean(EX,2);
    EX=((EX-EX_mean)*cPar.shrinkFactor)+EX_mean;
    Vcc(:,q)=EX(:);
end

Ec=reshape(1:size(Vcc,1),size(E,1),size(E,2));

[Fcc,~]=element2patch(Ec,[]);

switch cPar.latticeSide
    case 1  
        
        if strcmp(cPar.meshType,'hex')
            
            switch elementType
                case 'tet4'
                    ind1=1:size(E,1);
                    ind2=ind1+size(E,1);
                    ind3=ind2+size(E,1);
                    ind4=ind3+size(E,1);
             
                    Fn=[...
                        F(ind1,1)+size(Vc,1)+size(Vcc,1) Fc(ind1,1) Fc(ind1,2) F(ind1,2)+size(Vc,1)+size(Vcc,1) Fc(ind4,3) Fcc(ind1,1)+size(Vc,1) Fcc(ind1,2)+size(Vc,1) Fc(ind4,2);...                        
                        F(ind2,1)+size(Vc,1)+size(Vcc,1) Fc(ind2,1) Fc(ind2,2) F(ind2,2)+size(Vc,1)+size(Vcc,1) Fc(ind3,1) Fcc(ind2,1)+size(Vc,1) Fcc(ind2,2)+size(Vc,1) Fc(ind3,3);...  
                        F(ind3,1)+size(Vc,1)+size(Vcc,1) Fc(ind3,1) Fc(ind3,2) F(ind3,2)+size(Vc,1)+size(Vcc,1) Fc(ind4,1) Fcc(ind3,1)+size(Vc,1) Fcc(ind3,2)+size(Vc,1) Fc(ind4,3);...   
                        F(ind1,2)+size(Vc,1)+size(Vcc,1) Fc(ind1,2) Fc(ind1,3) F(ind1,3)+size(Vc,1)+size(Vcc,1) Fc(ind2,3) Fcc(ind1,2)+size(Vc,1) Fcc(ind1,3)+size(Vc,1) Fc(ind2,2);...  
                        F(ind1,3)+size(Vc,1)+size(Vcc,1) Fc(ind1,3) Fc(ind1,1) F(ind1,1)+size(Vc,1)+size(Vcc,1) Fc(ind3,3) Fcc(ind1,3)+size(Vc,1) Fcc(ind1,1)+size(Vc,1) Fc(ind3,2);...
                        F(ind2,3)+size(Vc,1)+size(Vcc,1) Fc(ind2,3) Fc(ind2,1) F(ind2,1)+size(Vc,1)+size(Vcc,1) Fc(ind4,2) Fcc(ind4,2)+size(Vc,1) Fcc(ind4,1)+size(Vc,1) Fc(ind4,1);...
                        ];

                    Cn=repmat((1:1:size(E,1))',[6 1]);
       
                    %Boundary parts
                    if ~isempty(cPar.indBoundary)
                        Fb=[F(cPar.indBoundary,1)+size(Vc,1)+size(Vcc,1) Fc(cPar.indBoundary,1) Fc(cPar.indBoundary,2) F(cPar.indBoundary,2)+size(Vc,1)+size(Vcc,1);...
                            F(cPar.indBoundary,2)+size(Vc,1)+size(Vcc,1) Fc(cPar.indBoundary,2) Fc(cPar.indBoundary,3) F(cPar.indBoundary,3)+size(Vc,1)+size(Vcc,1);...
                            F(cPar.indBoundary,3)+size(Vc,1)+size(Vcc,1) Fc(cPar.indBoundary,3) Fc(cPar.indBoundary,1) F(cPar.indBoundary,1)+size(Vc,1)+size(Vcc,1)];
                    end
                                        
                    Vn=[Vc;Vcc;V];
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
                     size(Cn)
                     
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
        switch elementType
            case 'tet4'
                switch cPar.meshType
                    case 'tri'
                        %Quad
                        Fn=[Fc(:,1) Fcc(:,1)+size(Vc,1) Fcc(:,2)+size(Vc,1) Fc(:,2);...
                            Fc(:,2) Fcc(:,2)+size(Vc,1) Fcc(:,3)+size(Vc,1) Fc(:,3);...
                            Fc(:,3) Fcc(:,3)+size(Vc,1) Fcc(:,1)+size(Vc,1) Fc(:,1)];
                        %Tri
                        Fn=[Fn(:,[1 2 3]); Fn(:,[3 4 1])];
                        
                        Cn=zeros(size(Fn,1),1);
                        
                        %Boundary parts
                        if ~isempty(cPar.indBoundary)
                            Fb=Fc(cPar.indBoundary,:);
                            Fn=[Fn;Fb];
                            Cn=[Cn;ones(size(Fb,1),1)];
                        end
                        
                        Vn=[Vc;Vcc;];
                        
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

%%

%Removing double vertices
[~,IND_V,IND_IND]=unique(pround(Vn,cPar.numDigitKeep),'rows');
Vn=Vn(IND_V,:);
if size(Fn,1)==1
    Fn=IND_IND(Fn)'; %Fix indices in F
else
    Fn=IND_IND(Fn); %Fix indices in F
end
    

