function [varargout]=image_numeric(varargin)

% function image_numeric(M)
% ------------------------------------------------------------------------
% Plots the intensities of an image (rounded to fit 4 digits) in the image
% at the pixel coordinates.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2008/08/28: Created
% 
% 
% ------------------------------------------------------------------------

%%

switch nargin
    case 1
        M=varargin{1};
        hf=gcf;
        numDigits=2;
        fontSize=10;
        textColor='k';
        interpreter='tex';
    case 2
        M=varargin{1};
        hf=varargin{2};
        numDigits=2;
        fontSize=10;
        textColor='k';
        interpreter='tex';
    case 3
        M=varargin{1};
        hf=varargin{2};
        numDigits=varargin{3};
        fontSize=10;
        textColor='k';
        interpreter='tex';
    case 4        
        M=varargin{1};
        hf=varargin{2};
        numDigits=varargin{3};
        fontSize=varargin{4};
        textColor='k';
        interpreter='tex';
    case 5        
        M=varargin{1};
        hf=varargin{2};
        numDigits=varargin{3};
        fontSize=varargin{4};
        textColor=varargin{5};
        interpreter='tex';
    case 6
        M=varargin{1};
        hf=varargin{2};
        numDigits=varargin{3};
        fontSize=varargin{4};
        textColor=varargin{5};
        interpreter=varargin{6};
end

if isempty(hf)
    hf=gcf;
end

% textFormat=['%6.',num2str(numDigits),'e'];
textFormat=['%0.',num2str(numDigits),'f';];                    
 
%%
figure(hf);
H=zeros(1,numel(M));
qh=1;
for qi=1:size(M,1)
    for qj=1:size(M,2)
        m=M(qi,qj);        
        switch class(m)
            case'sym'
                try %to convert to double
                    m=double(m);                    
                    image_text=sprintf(textFormat,m);
                catch %Contains symbolic expressions
                    switch interpreter
                        case 'none'
                            image_text=char(m);
                        case 'latex'
                            image_text=latex(m);
                        case 'tex'
                            image_text=texlabel(m);
                    end
                    %image_text = strsplit(image_text,' '); %Places space
                    %separated elements on new line
                end
            otherwise                
                image_text=sprintf(textFormat,m);
        end
        H(qh)=text(qj,qi, image_text,'horizontalAlignment','center','color',textColor,'FontWeight','demi','FontSize',fontSize,'interpreter',interpreter);
        qh=qh+1;
    end
end
if nargout==1
    varargout{1}=H;
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
