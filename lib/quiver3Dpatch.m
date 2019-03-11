function [varargout]=quiver3Dpatch(x,y,z,ux,uy,uz,c,a)

% function [F,V,C]=quiver3Dpatch(x,y,z,ux,uy,uz,c,a)
% ------------------------------------------------------------------------
%
% This function allows plotting of colored 3D arrows by generating patch
% data (faces �F�, vertices �V� and color data �C�). The patch data which
% allows plotting of 3D quiver arrows with specified (e.g. colormap driven)
% color. To save memory n arrows are created using only n*6 faces and n*7
% vertices. The vector "a" defines arrow length scaling where a(1) is the
% smallest arrow length and a(2) the largest. Use the PATCH command to plot
% the arrows:  
%
% [F,V,C]=quiver3Dpatch(x,y,z,ux,uy,uz,a)
% patch('Faces',F,'Vertices',V,'CData',C,'FaceColor','flat'); 
%
% Below is an example illustrating color specifications for (combined)
% patch data.  
%
%%% EXAMPLE
% % Simulating 3D volume and vector data
% n=25;
% [X,Y,Z]=meshgrid(linspace(-4.77,4.77,n));
% phi=(1+sqrt(5))/2;
% M=2 - (cos(X + phi*Y) + cos(X - phi*Y) + cos(Y + phi*Z) + cos(Y - phi*Z) + cos(Z - phi*X) + cos(Z + phi*X));
% 
% % Simulating vector data 
% %Vector data here based on the gradient of the image
% [u,v,w] = gradient(M); 
% G=hypot(hypot(u,v),w); %Vector lenghts
% 
% a=[min(G(:)) max(G(:))]; %Arrow length scaling
% L=G>0.7; %Logic indices for arrows to make a selection if required
% [Fv,Vv,Cv]=quiver3Dpatch(X(L),Y(L),Z(L),u(L),v(L),w(L),G(L),a);
% 
% figure; 
% xlabel('X','FontSize',10);ylabel('Y','FontSize',10);zlabel('Z','FontSize',10);
% title('Colormap driven vectors','FontSize',15);
% patch('Faces',Fv,'Vertices',Vv,'EdgeColor','none', 'CData',Cv,'FaceColor','flat','FaceAlpha',1); 
% 
% colormap jet; colorbar; 
% % caxis([min(Cv(:)) max(Cv(:))]);
% view(3); grid on; axis equal; axis vis3d; 
% set(gca,'FontSize',10);
% camlight headlight; lighting flat;
% drawnow;
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2014/01/13 Updated example in help
% 2014/11/11 Updated to allow for RGB color input data (e.g. C is an nx3 array)
% 2016/06/30 Added handling of empty "a" parameter. Reverts to default
% normal length range. 
%------------------------------------------------------------------------

%% 

%Convert to columns
x=x(:); y=y(:); z=z(:); 
ux=ux(:); uy=uy(:); uz=uz(:); 

%Spherical coordinates
[THETA_vec,PHI_vec,R_vec] = cart2sph(ux,uy,uz);
  
%% Setting arrow size properties

if isempty(a)
    a=[min(R_vec(:)) max(R_vec(:))];
end

%Arrow length
if  min(R_vec(:))==max(R_vec(:)) %If all radii are equal, or if just 1 vector is used
    arrow_length=a(2)*ones(size(R_vec)); %All arrow lengths become a(2)
else
    %Scale arrow lengths between a(1) and a(2)
    arrow_length=R_vec-min(R_vec(:));
    arrow_length=a(1)+((arrow_length./max(arrow_length(:))).*(a(2)-a(1)));
end

%Other arrow dimensions as functions of arrow length and phi (golden ratio)
phi=(1+sqrt(5))/2;
head_size=arrow_length./(phi);
head_width=head_size./(2.*phi);
stick_width=head_width./(2.*phi);
if numel(a)==3
    head_width=a(3)*head_width;
    stick_width=a(3)*stick_width;
end
h=sin((2/3).*pi).*stick_width;
ha=sin((2/3).*pi).*head_width;

%% Creating arrow triangle vertices coordinates

X_tri=[zeros(size(x))  zeros(size(x)) zeros(size(x))...
    head_size.*ones(size(x)) head_size.*ones(size(x)) head_size.*ones(size(x))...
    arrow_length];
Y_tri=[-(0.5.*stick_width).*ones(size(x)) (0.5.*stick_width).*ones(size(x))  zeros(size(x))...
    -(0.5.*head_width).*ones(size(x))  (0.5.*head_width).*ones(size(x))  zeros(size(x))...
    zeros(size(x))];
Z_tri=[-(0.5.*stick_width.*tan(pi/6)).*ones(size(x))...
    -(0.5.*stick_width.*tan(pi/6)).*ones(size(x))...
    (h-(0.5.*stick_width.*tan(pi/6))).*ones(size(x))...
    -(0.5.*head_width.*tan(pi/6)).*ones(size(x))...
    -(0.5.*head_width.*tan(pi/6)).*ones(size(x))...
    (ha-(0.5.*head_width.*tan(pi/6))).*ones(size(x))...
    zeros(size(x))];

% Rotating vertices
[THETA_ar,PHI_ar,R_vec_ar] = cart2sph(X_tri,zeros(size(Y_tri)),Z_tri);
PHI_ar=PHI_ar+PHI_vec*ones(1,size(THETA_ar,2));
[X_arg,Y_arg,Z_arg] = sph2cart(THETA_ar,PHI_ar,R_vec_ar);
Y_arg=Y_arg+Y_tri;
[THETA_ar,PHI_ar,R_vec_ar] = cart2sph(X_arg,Y_arg,Z_arg);
THETA_ar=THETA_ar+THETA_vec*ones(1,size(THETA_ar,2));
[X_arg,Y_arg,Z_arg] = sph2cart(THETA_ar,PHI_ar,R_vec_ar);

X_arg=X_arg+x*ones(1,size(THETA_ar,2)); X_arg=X_arg';
Y_arg=Y_arg+y*ones(1,size(THETA_ar,2)); Y_arg=Y_arg';
Z_arg=Z_arg+z*ones(1,size(THETA_ar,2)); Z_arg=Z_arg';

V=[X_arg(:) Y_arg(:) Z_arg(:)];

%% Creating faces matrix

%Standard vertex order for 6 face arrow style
F_order=[1 2 7; 2 3 7; 3 1 7; 4 5 7; 5 6 7; 6 4 7;];
no_nodes=size(X_tri,2);
b=(no_nodes.*((1:1:numel(x))'-1))*ones(1,3);

%Loops size(F_order,1) times
F=zeros(numel(x)*size(F_order,1),3); %Allocating faces matrix
for f=1:1:size(F_order,1)
    Fi=ones(size(x))*F_order(f,:)+b;
    F(1+(f-1)*numel(x):f*numel(x),:)=Fi;
end

%     %Alternative without for loop, not faster for some tested problems
%     F_order=(ones(numel(x),1)*[1 2 7 2 3 7 3 1 7 4 5 7 5 6 7 6 4 7])';
%     F_order1=ones(numel(x),1)*(1:6);
%     F_order2=ones(numel(x),1)*[2 3 1 5 6 4];
%     F_order=[F_order1(:) F_order2(:)];
%     F_order(:,3)=7;
%     b=repmat(((no_nodes.*(0:1:numel(x)-1)')*ones(1,3)),[6,1]);
%     F=F_order+b;

%% Collect face and vertex output
varargout{1}=F;
varargout{2}=V;

%% Add color specification if requested

if nargout==3
    if isempty(c) %If empty specify vector magnitudes as colormap driven color
        C=repmat(R_vec,[size(F_order,1),1]);
    else %If user specified color replicate to match # of added faces for arrow
        %c my be an nx3 array to allow for RGB type color data
        C=repmat(c,[size(F_order,1),1]);
    end
    varargout{3}=C;
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
