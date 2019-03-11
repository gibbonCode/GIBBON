function H=errorbarC_XY(x,y,em,ep,w,C,d)
hold on;
switch d
    case 1
        y_u=y+ep;
        y_l=y-em;
        y_v=[y_u(:) y_l(:)];
        x_v=[x(:) x(:) ];

        y_ul=y_u;
        y_ur=y_u;
        y_ll=y_l;
        y_lr=y_l;
        x_ul=x-0.5*w;
        x_ur=x+0.5*w;
        x_ll=x-0.5*w;
        x_lr=x+0.5*w;
        y_hu=[y_ul(:) y_ur(:)];
        x_hu=[x_ul(:) x_ur(:)];
        y_hl=[y_ll(:) y_lr(:)];
        x_hl=[x_ll(:) x_lr(:)];

        H=[];
        for i=1:1:numel(x);
            %Points
            h1=plot(x(i),y(i),'b.');
            %Vertical lines
            h2=plot(x_v(i,:),y_v(i,:),'b-');
            %Horizontal lines upper
            h3=plot(x_hu(i,:),y_hu(i,:),'b-');
            %Horizontal lines lower
            h4=plot(x_hl(i,:),y_hl(i,:),'b-');

            %Set color
            set([h1;h2;h3;h4],'Color',C(i,:));

            H=[H;h1;h2;h3;h4];
        end
    case 2
        x_r=x+ep;
        x_l=x-em;
        
        y_h=[y(:) y(:)];
        x_h=[x_l(:) x_r(:) ];

        x_ul=x_l;
        x_ur=x_r;
        x_ll=x_l;
        x_lr=x_r;
        y_ul=y+0.5*w;
        y_ur=y+0.5*w;
        y_ll=y-0.5*w;
        y_lr=y-0.5*w;
        
        y_vl=[y_ul(:) y_ll(:)];
        x_vl=[x_ul(:) x_ll(:)];
        y_vr=[y_ur(:) y_lr(:)];
        x_vr=[x_ur(:) x_lr(:)];

        H=[];
        for i=1:1:numel(x);
            %Points
            h1=plot(x(i),y(i),'b.');
            %Horintal lines
            h2=plot(x_h(i,:),y_h(i,:),'b-');
            %Vertical lines left
            h3=plot(x_vl(i,:),y_vl(i,:),'b-');
            %Vertical lines left
            h4=plot(x_vr(i,:),y_vr(i,:),'b-');
            %Set color
            set([h1;h2;h3;h4],'Color',C(i,:));

            H=[H;h1;h2;h3;h4];
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
