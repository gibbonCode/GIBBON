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
 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
