function [glyphImage]=textImage(varargin)

%%

switch nargin
    case 0
        textInput='GIBBON';
        FontName=vision.internal.getDefaultFont();
        FontSize=50;
        padAmount=0;
    case 1
        textInput=varargin{1};
        FontName=vision.internal.getDefaultFont();
        FontSize=50;
        padAmount=0;
    case 2
        textInput=varargin{1};
        FontName=varargin{2};
        FontSize=50;
        padAmount=0;
    case 3        
        textInput=varargin{1};
        FontName=varargin{2};
        FontSize=varargin{3};        
        padAmount=0;        
    case 4
        textInput=varargin{1};
        FontName=varargin{2};
        FontSize=varargin{3};
        padAmount=varargin{4};
end

position=[0 0];

if isa(textInput,'cell')
    textString=char(textInput);
else
    textString=textInput;
end

%%
FONT=listTrueTypeFonts(FontName);

try
[glyphBitmapArray, ...
    glyphIdxFromCharcode, ...
    glyphBitmapStartIdx, ...
    glyphWidths, ...
    glyphHeights, ...
    glyphXAdvances, ...
    glyphLeftBearings, ...
    glyphTopBearings, ...
    fontAscend, ...
    fontDescend, ...
    fontLinespace, ...
    maxBitmapSize]=visionPopulateGlyphBuffer(FONT.fileName,FONT.faceIndex, ...
    FontSize,false);
catch ME
    warning(['Check if ',FontName,' is a valid font use listTrueTypeFonts to see available fonts']);
    rethrow(ME);
end

glyphStruct.glyphBitmapArray     = glyphBitmapArray;
glyphStruct.glyphIdxFromCharcode = glyphIdxFromCharcode;
glyphStruct.glyphBitmapStartIdx  = glyphBitmapStartIdx;
glyphStruct.glyphWidths          = glyphWidths;
glyphStruct.glyphHeights         = glyphHeights;
glyphStruct.glyphXAdvances       = glyphXAdvances;
glyphStruct.glyphLeftBearings    = glyphLeftBearings;
glyphStruct.glyphTopBearings     = glyphTopBearings;
%
glyphStruct.fontAscend            = fontAscend;
glyphStruct.fontDescend           = fontDescend;
glyphStruct.fontLinespace         = fontLinespace;

fontStruct.fontAscend            = fontAscend;
fontStruct.fontDescend           = fontDescend;
fontStruct.fontLinespace         = fontLinespace;

glyphStruct.maxBitmapSize         = maxBitmapSize;

%%

glyphIdxFromCharcode=glyphStruct.glyphIdxFromCharcode;
glyphXAdvances=glyphStruct.glyphXAdvances;
fontHeightWLinegap  = fontStruct.fontLinespace;
fontHeightWOLinegap = fontStruct.fontAscend - fontStruct.fontDescend;


b1 = 1; % for converting 0-based to 1-based
spaceGlyphIdx = glyphIdxFromCharcode(32+b1);
if (spaceGlyphIdx==0)
    spaceCharWidth = int32(fontHeightWOLinegap/4);
else
    spaceCharWidth = int32(glyphXAdvances(spaceGlyphIdx+b1));
end


MARGIN_LeftOrRight = spaceCharWidth; % used 3 before
MARGIN_TopOrBottom = spaceCharWidth; % used 3 before

%%

position = int32(position);

tbLocationXY = int32(position);

tbLocationX = tbLocationXY(1);
tbLocationY = tbLocationXY(2);

tbTopLeftY = tbLocationY;
tbTopLeftX = tbLocationX;

textLocationXY.x = tbTopLeftX + int32(MARGIN_LeftOrRight);
textLocationXY.y = tbTopLeftY + int32(MARGIN_TopOrBottom);


%%

textStringU16 = uint16(textString);

penX = int32(textLocationXY.x);

% go to reference baseline (near the middle of the glyph)
penY = int32(textLocationXY.y) + fontStruct.fontAscend;

oneI32 = int32(1);

isNewLineChar = (textStringU16 == uint16(10));

numChars = numel(textStringU16);

glyphImages=cell(1,numChars);
glyph_I=zeros(numChars,2);
glyph_J=zeros(numChars,2);
q=1;
for ii=1:size(textStringU16,1)    
    %Set position
    if ii>1
        penY = penY + fontHeightWLinegap; % go to next line
    end
    penX = textLocationXY.x;
    for jj=1:size(textStringU16,2)
        
        %see logic in mdlOutputs of sviptextrender.cpp
        
        % reset x position to the beginning on a line
        
        
        thisCharcode = textStringU16(ii,jj);
        thisCharcode_1b = thisCharcode+b1;
        thisGlyphIdx = glyphStruct.glyphIdxFromCharcode(thisCharcode_1b);
        thisGlyphIdx_1b = thisGlyphIdx + b1;
        glyphExists = (thisGlyphIdx ~= 0);
        if ~glyphExists
            penX = penX + int32(spaceCharWidth);
        else
            thisGlyphW = glyphStruct.glyphWidths(thisGlyphIdx_1b);
            thisGlyphH = glyphStruct.glyphHeights(thisGlyphIdx_1b);
            
            xx=penX+int32(glyphStruct.glyphLeftBearings(thisGlyphIdx_1b));
            yy=penY-int32(glyphStruct.glyphTopBearings(thisGlyphIdx_1b));
            
            start_I = yy;
            end_I = yy+int32(thisGlyphH)- oneI32;
            start_J = xx;
            end_J = xx+int32(thisGlyphW)- oneI32;
            
            bitmapStartIdx_1b = glyphStruct.glyphBitmapStartIdx(thisGlyphIdx_1b) + b1;
            bitmapEndIdx_1b =  bitmapStartIdx_1b + uint32(thisGlyphW*thisGlyphH) - uint32(1);
            thisGlyphBitmap = glyphStruct.glyphBitmapArray(bitmapStartIdx_1b:bitmapEndIdx_1b);
            thisGlyphBitmap = reshape(thisGlyphBitmap,[thisGlyphW thisGlyphH])';
            
            glyphImages{q}=thisGlyphBitmap;
            glyph_I(q,:)=[start_I end_I];
            glyph_J(q,:)=[start_J end_J];
            
            % update X position for next character
            penX=penX+int32(glyphStruct.glyphXAdvances(thisGlyphIdx_1b));
            
        end
        q=q+1;
    end
end

%Crop around text
glyph_I=glyph_I-min(glyph_I(:))+1;
glyph_J=glyph_J-min(glyph_J(:))+1;

%% Create image

glyphImage=zeros(max(glyph_I(:))+min(glyph_I(:))-1,max(glyph_J(:))+min(glyph_J(:))-1);
for q=1:1:numChars
    glyphImage(glyph_I(q,1):glyph_I(q,2),glyph_J(q,1):glyph_J(q,2))=glyphImages{q};
end
glyphImage=double(glyphImage);

if padAmount>0    
    siz=size(glyphImage)+2*padAmount;
    glyphImagePad=zeros(siz);
    glyphImagePad(padAmount+1:padAmount+size(glyphImage,1),padAmount+1:padAmount+size(glyphImage,2))=glyphImage;
    glyphImage=glyphImagePad;
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
