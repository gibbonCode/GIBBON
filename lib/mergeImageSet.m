function [N,G,hf]=mergeImageSet(Nc,Gc,nFac,plotOn)

%Initialise N and G using the first image set
N=Nc{1}; 
G=Gc{1};
for q=2:1:numel(Nc)
    N2=Nc{q}; 
    G2=Gc{q};     
    [N,G]=mergeImageData(N,G,N2,G2,nFac,0);
end

if plotOn==1; 
    %Plot settings

    fontSize=20;    
    tP=1; %Pause time

    %Open figure
    hf=cFigure;  set(gcf,'renderer','OpenGL');
    xlabel('X (mm)','FontSize',fontSize); ylabel('Y (mm)','FontSize',fontSize); zlabel('Z (mm)','FontSize',fontSize);
    hold on;
    
    %Plotting results
    hp1=[]; 
    for q=1:size(N,4);
        title(['Dynamic ',num2str(q),' of ',num2str(size(N,4))],'FontSize',fontSize);
        delete(hp1); %Delete previous plot handle       
        M=N(:,:,:,q); %Define image to plot
        
        %Define slice selection logic
        L=false(size(M)); 
        L(round(size(M,1)/2),:,:)=1; 
        L(:,round(size(M,2)/2),:)=1; 
        L(:,:,round(size(M,3)/2))=1;        
        
        %Define patch data
        [F,V,C]=ind2patch(L,M,'vb'); %Creating patch data for selection of low voxels
        
        %Convert coordinates
        [V(:,1),V(:,2),V(:,3)]=im2MRcart(V(:,2),V(:,1),V(:,3),G.v,G.OR,G.r,G.c); %Convert image to cartesian coordinates
        
        %Patch the image data
        hp1= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','none');
        
        %Axis settings
        axis equal; view(3); axis tight;  grid on;
        colormap('gray'); colorbar;
        set(gca,'FontSize',fontSize);
        drawnow;
        pause(tP); 
    end 
%     pause(tP);
%     close(hf);
else
    hf=NaN;
end

%TO DO: Proper memory allocation for N, currently it grows in the loop
 
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
