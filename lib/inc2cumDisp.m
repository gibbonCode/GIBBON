function [V_path]=inc2cumDisp(V1_C,D_C,V1i,Method,ExtrapolationMethod,siz,plotOn)

if plotOn==1    
    font_size=25;

    hf1=cFigure;;
    subplot(1,2,1);
    xlabel('X (mm)','FontSize',font_size);ylabel('Y (mm)','FontSize',font_size); %zlabel('Z','FontSize',font_size);
    hold on;
    plotV(V1i,'k.');
    view(2); grid on; axis equal; axis tight;
    drawnow;
    
    subplot(1,2,2);
    xlabel('X (mm)','FontSize',font_size);ylabel('Y (mm)','FontSize',font_size); %zlabel('Z','FontSize',font_size);
    hold on;
    plotV(V1i,'k.');
    view(2); grid on; axis equal; axis tight;
    drawnow;
    
    hp1=[]; hp2=[]; hp3=[]; hp4=[];
end


V_path=repmat(V1i,[1 1 numel(V1_C)+1]);
for q=1:numel(V1_C)

    %% KNOWN INITIALS AND MATCHED DISPLACEMENTS
    v1=V1_C{q}; %Initials where we know displacement
    d=D_C{q}; %Displacement for these initials
    
    %Remove invalid points (NaNs)
    L=any(isnan(v1),2) | any(isnan(d),2);
    v1=v1(~L,:);
    d=d(~L,:);
    
    %% POINTS AT WHICH DISPLACEMENT ASSESSMENT IS DESIRED
    v1i=V_path(:,:,q); %The current path ending is where we wish to interpolate
    
    %Remove invalid points (NaNs)
    L=any(isnan(v1i),2);
    v1ii=v1i(~L,:);
    
    %Plotting
    if plotOn==1
        v2=v1+d;
        c=sqrt(sum(d.^2,2));
        [Fs,Vs,Cs]=quiver3Dpatch(v1(:,1),v1(:,2),v1(:,3),d(:,1),d(:,2),d(:,3),c,[min(c) max(c)]);
        subplot(1,2,1);
        delete(hp1);
        hp1=patch('Faces',Fs,'Vertices',Vs,'EdgeColor','none', 'CData',Cs,'FaceColor','flat','FaceAlpha',1);
        delete(hp2);
        hp2=plotV(v2,'r.');
        view(2); grid on; axis equal; axis tight;
        colorbar;
        drawnow;
    end
    
    %% INTERPOLATE FOR DESIRED POINTS
    
    %Construct interpolants
    F_DX = scatteredInterpolant(v1,d(:,1),Method,ExtrapolationMethod);
    F_DY = scatteredInterpolant(v1,d(:,2),Method,ExtrapolationMethod);
    F_DZ = scatteredInterpolant(v1,d(:,3),Method,ExtrapolationMethod);
    
    %Interpolate displacement
    dd=[F_DX(v1ii) F_DY(v1ii) F_DZ(v1ii)];    
    
    %Store results in path
    v1=v1i;
    d=nan(size(v1));
    d(~L,:)=dd;    
    
    if ~isempty(siz)
%         f=ones(3,3);f=f./sum(f(:)); %Averaging filter
        f=gauss_kernel(3,2,2,'width'); 
        m=sqrt(sum(d.^2,2));
        %Reshaping results to mesh form
        dx=reshape(d(:,1),siz);
        dy=reshape(d(:,2),siz);
        dz=reshape(d(:,3),siz);
        
        %Filtering
        Xc=convn(dx,f,'same');
        Yc=convn(dy,f,'same');
        Zc=convn(dz,f,'same');
        
        d=[Xc(:) Yc(:) Zc(:)];
        
        %Push back magnitude
        m2=sqrt(sum(d.^2,2));
        d=m(:,ones(1,3)).*(d./m2(:,ones(1,3)));
    end
    v2=v1+d;
    V_path(:,:,q+1)=v2;
    
    %Plotting
    if plotOn==1
        c=sqrt(sum(d.^2,2));
        [Fs,Vs,Cs]=quiver3Dpatch(v1(:,1),v1(:,2),v1(:,3),d(:,1),d(:,2),d(:,3),c,[min(c) max(c)]);
        subplot(1,2,2);
%         delete(hp3);
        hp3=patch('Faces',Fs,'Vertices',Vs,'EdgeColor','none', 'CData',Cs,'FaceColor','flat','FaceAlpha',1);
        delete(hp4);
        hp4=plotV(v2,'k.');
        view(2); grid on; axis equal; axis tight;
        colorbar;
        drawnow;
    end
end

if plotOn==1; 
    close(hf1); 
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
