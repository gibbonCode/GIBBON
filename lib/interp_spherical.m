function [Ri]=interp_spherical(T,P,R,Ti,Pi,interpMethod,numberInterpSteps)

%%
demoMode=(nargin==0);
if demoMode
    disp('DEMO MODE');
    [hp2,Fp,Tp,Pp]=demoFunction1;
end

%%
if numberInterpSteps==1
    [Ri]=interp_spherical_part(T,P,R,Ti,Pi,interpMethod);
else
    %Getting source vertices
    V=zeros(size(R,1),3);
    [V(:,1),V(:,2),V(:,3)] = sph2cart(T,P,R);
    
    %Getting target vertices
    Vi=zeros(size(Pi,1),3);
    [Vi(:,1),Vi(:,2),Vi(:,3)] = sph2cart(Ti,Pi,ones(size(Pi)));
    
    %Determine the rotation settings
    rotationRange=linspace(0,pi,numberInterpSteps+1);
    rotationRange=rotationRange(1:end-1);
    phiThreshold=diff(rotationRange(1:2))/2;
    
    %Step wise spherical interpolation around equator
    indDone=[];
    Ri=ones(size(Vi,1),1);
    for q=1:1:numel(rotationRange)
        %Rotate so that points of interest are at equator
        [DCM,~]=euler2DCM([rotationRange(q) 0 0]); %Direction cosines matrix
        Vr=V*DCM; %Rotated source vertices
        Vri=Vi*DCM; %Rotated target vertices
        
        %Convert to spherical coordinates
        [Tr,Pr,Rr]=cart2sph(Vr(:,1),Vr(:,2),Vr(:,3)); %Source
        [Ti,Pi,~]=cart2sph(Vri(:,1),Vri(:,2),Vri(:,3)); %Target
        
        %Get indices for the current target vertices
        indDoNow=find(Pi>=-phiThreshold & Pi<=phiThreshold);
        indDoNow=indDoNow(~ismember(indDoNow,indDone)); %remove ones done already
        
        %Interpolate region
        [Ri_step]=interp_spherical_part(Tr,Pr,Rr,Ti(indDoNow),Pi(indDoNow),interpMethod);
        
        Ri(indDoNow)=Ri_step; %Setting new radii
        
        indDone=unique([indDone; indDoNow]); %Adding vertix indices to done list
        
        if demoMode
            p_color=hsv(numel(rotationRange));
            [hp2]=demoFunction2(hp2,Fp,Tp,Pp,p_color);
        end
    end
end


%%

    function [Ri]=interp_spherical_part(T,P,R,Ti,Pi,interpMethod)
        %Tesselate data above, below, left and right to aid interpolation (and
        %avoid some polar- and edge artifacts)
        T=[T+2*pi; T+2*pi; T+2*pi; T; T; T; T-2*pi; T-2*pi; T-2*pi];
        P=[P-pi; P; P+pi; P-pi; P; P+pi; P-pi; P; P+pi];
        R=repmat(R,[9,1]);
        
        %Removing double points
        fRound=1e5; %Rounding factor for unique test
        [~,indUni,~]=unique(round([T P R]*fRound)/fRound,'rows');
        P=P(indUni);
        T=T(indUni);
        R=R(indUni);
        
        if strcmp(interpMethod,'natural') || strcmp(interpMethod,'linear') || strcmp(interpMethod,'nearest') %TriScatterdInterp function
            F_delaunay=scatteredInterpolant([T P],R,interpMethod); %interpolator
            Ri=F_delaunay([Ti Pi]);
        elseif strcmp(interpMethod,'cubic') %Griddata function
            Ri = griddata(T,P,R,Ti,Pi,interpMethod);
        else
            %Add warning here
        end
    end

%% DEMO functions
    function [hp2,Fp,Tp,Pp]=demoFunction1
        %% Plot settings
        figColor='w'; 
        figColorDef='white';
        font_size=10;
        
        %% Creating source vertices
        %Start with sphere triangulation (Buckminster-Fuller type)
        [F,V,Vs]=geoSphere(3,1);
        T=Vs(:,1); P=Vs(:,2); %Spherical coordinates
        
        %Modifying to obtain complex shape
        [DCM,~]=euler2DCM([0.25.*pi -0.25*pi 0]); %Rotate
        Vr=(V*DCM);
        [Tr,Pr,Rr]=cart2sph(Vr(:,1),Vr(:,2),Vr(:,3)); %convert to spherical coordinates
        [V(:,1),V(:,2),V(:,3)]=sph2cart(T,P,1+0.25*cos(2*(Tr))+0.25*cos(8*(P))); %Convert back with unity radii
        V(:,3)=V(:,3)*2; %Scale
        [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3)); %Convert to spherical coordinates        
        
        %% Creating target vertices
        
        %Defining denser point set for interpolation
        [Fi,Vi,Vsi]=geoSphere(4,1);
        Ti=Vsi(:,1); Pi=Vsi(:,2); %Spherical coordinates

        %%
        
        figuremax(figColor,figColorDef);        
        for q1=1:3
            subplot(1,3,q1);
            hold on;
            xlabel('X','FontWeight','bold','FontSize',font_size); ylabel('Y','FontWeight','bold','FontSize',font_size); zlabel('Z','FontWeight','bold','FontSize',font_size);
            view(3); grid on; 
            set(gca,'FontSize',font_size);
        end
        
        subplot(1,3,1);
        title('Coarse original form');
        hp1=patch('Faces',F,'Vertices',V);
        set(hp1,'FaceColor','g','EdgeColor','w');
        axis equal tight vis3d;
        camlight('headlight'); lighting phong;
        drawnow;
        ha=axis;
        
        subplot(1,3,2);
        hp1=patch('Faces',Fi,'Vertices',Vi);
        set(hp1,'FaceColor',0.5*ones(1,3),'EdgeColor','w');
        axis equal tight vis3d;
        camlight('headlight'); lighting phong;
        drawnow;
        
        subplot(1,3,3);
        title('Interpolated form');
        hp2=patch('Faces',Fi,'Vertices',Vi);
        set(hp2,'FaceColor','r','EdgeColor','none'); lighting phong;
        axis equal tight vis3d;
        axis(ha);        
        camlight('headlight'); lighting phong;
        drawnow;      
                
        numberInterpSteps=8; %Number of interpolation steps
        interpMethod='cubic';
        
        Fp=Fi; Tp=Ti; Pp=Pi;
        
    end

    function [hp2]=demoFunction2(hp2,Fp,Tp,Pp,p_color)
                
        subplot(1,3,2); colormap(p_color); colorbar;caxis([1-0.5,size(p_color,1)+0.5]);
        
        [V2i(:,1),V2i(:,2),V2i(:,3)]=sph2cart(Tp,Pp,Ri);
        
        %Plotting
        subplot(1,3,3);
        delete(hp2);
        hp2=patch('Faces',Fp,'Vertices',V2i);
        set(hp2,'FaceColor','r','EdgeColor','none'); lighting phong;
        
        subplot(1,3,2);
        title({'Visualizing interpolation steps'; [num2str(q),' of ',num2str(numel(rotationRange))]});
        
        hp3=plot3(Vi(indDoNow,1),Vi(indDoNow,2),Vi(indDoNow,3),'k.');
        set(hp3,'markersize',20,'color',p_color(q,:));
        axis tight;
        drawnow;
    end
end

 
%% <-- GIBBON footer text --> 
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
