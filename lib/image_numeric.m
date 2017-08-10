function [varargout]=image_numeric(varargin)

% function image_numeric(M)
% ------------------------------------------------------------------------
% Plots the intensities of an image (rounded to fit 4 digits) in the image
% at the pixel coordinates.
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 28/08/2008
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
