function [mousePointerType]=specialPointerShape(pointerType)

switch pointerType
    case 'smallHand'
        PointerShapeCData=[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,2,2,NaN,NaN,NaN,NaN,NaN,1,2,2,NaN,NaN,NaN,NaN,NaN;NaN,2,2,2,NaN,NaN,1,2,NaN,1,2,2,NaN,NaN,NaN,NaN;NaN,1,2,2,2,2,NaN,1,2,2,2,2,NaN,NaN,NaN,NaN;NaN,NaN,1,2,2,2,2,2,2,2,1,2,2,NaN,NaN,NaN;NaN,NaN,NaN,1,2,2,2,2,1,2,2,1,2,NaN,NaN,NaN;NaN,NaN,NaN,NaN,1,2,2,2,2,1,2,2,2,NaN,NaN,NaN;NaN,NaN,NaN,1,1,1,2,2,2,2,2,2,2,NaN,NaN,NaN;NaN,NaN,NaN,1,1,1,2,2,2,2,2,2,NaN,2,2,NaN;NaN,NaN,NaN,NaN,1,1,2,2,2,2,2,2,2,2,2,NaN;NaN,NaN,NaN,NaN,NaN,1,1,1,1,1,NaN,2,2,1,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,2,1,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,1,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
        PointerShapeHotSpot=[3,3];
    case 'ulc' %Upper left corner
        PointerShapeCData=makeUpperLeftCornerPointer;
        PointerShapeHotSpot=[3,3];
    case 'llc' %Lower left corner
        PointerShapeCData=rot90(makeUpperLeftCornerPointer,1);
        PointerShapeHotSpot=[16-2,3];
    case 'lrc' %Lower right corner
        PointerShapeCData=rot90(makeUpperLeftCornerPointer,2);
        PointerShapeHotSpot=[16-2,16-2];
    case 'urc' %Upper right corner
        PointerShapeCData=rot90(makeUpperLeftCornerPointer,3);
        PointerShapeHotSpot=[3,16-2];
    case 'x'
        PointerShapeCData=2*eye(16,16); %2's on diagonal
        PointerShapeCData=max(PointerShapeCData,2.*diag(ones(1,16-1),-1)); %2's off diagonal
        PointerShapeCData=max(PointerShapeCData,2.*diag(ones(1,16-1),1)); %2's off diagonal
        PointerShapeCData=max(PointerShapeCData,diag(ones(1,16-2),2)); %1's on 2nd diagonal
        PointerShapeCData=max(PointerShapeCData,diag(ones(1,16-2),-2)); %1's on 2nd diagonal
        PointerShapeCData=max(PointerShapeCData,fliplr(PointerShapeCData));
        PointerShapeCData(PointerShapeCData==0)=nan;
        PointerShapeHotSpot=[8,8];
    case 'cut'
        PointerShapeCData=[2,2,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2;2,2,2,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2,2;1,2,2,2,1,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2,2,1;NaN,1,2,2,2,1,NaN,NaN,NaN,NaN,1,2,2,2,1,NaN;NaN,NaN,1,2,2,2,1,NaN,NaN,1,2,2,2,1,NaN,NaN;NaN,NaN,NaN,1,2,2,2,1,1,2,2,2,1,NaN,NaN,NaN;NaN,NaN,NaN,NaN,1,2,2,2,2,2,2,1,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,1,2,2,2,2,1,NaN,NaN,NaN,NaN,NaN;1,2,2,2,NaN,1,2,2,2,2,1,NaN,2,2,2,1;2,1,2,2,1,2,2,2,2,2,2,1,2,2,1,2;2,1,1,2,2,2,2,1,1,2,2,2,2,1,1,2;2,1,1,2,2,2,1,NaN,NaN,1,2,2,2,1,1,2;2,1,2,2,2,1,NaN,NaN,NaN,NaN,1,2,2,2,1,2;2,2,2,2,1,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2,2,2;2,2,2,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2,2;2,2,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2];
        PointerShapeHotSpot=[8,8];
    case 'pen'
        PointerShapeCData=[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2,2,1;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2,2,1,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2,2,2,1,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2,2,1,2,1,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,2,2,2,1,NaN,2,1,NaN;NaN,NaN,NaN,NaN,NaN,NaN,1,2,2,2,1,NaN,2,1,1,NaN;NaN,NaN,NaN,NaN,NaN,1,2,2,2,1,NaN,2,1,1,NaN,NaN;NaN,NaN,NaN,NaN,1,2,2,2,1,NaN,2,1,1,NaN,NaN,NaN;NaN,NaN,NaN,1,2,2,2,1,NaN,2,1,1,NaN,NaN,NaN,NaN;NaN,NaN,1,2,2,2,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,1,2,2,2,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,1,2,2,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;1,2,1,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;2,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
        PointerShapeHotSpot=[16,1];
end

mousePointerType.PointerShapeCData=PointerShapeCData;
mousePointerType.PointerShapeHotSpot=PointerShapeHotSpot;

end

function PointerShapeCData=makeUpperLeftCornerPointer
PointerShapeCData=ones(16,16);
PointerShapeCData(2:end,2:end)=2;
PointerShapeCData(3:end,3:end)=2;
PointerShapeCData(4:end,4:end)=2;
PointerShapeCData(5:end,5:end)=1;
PointerShapeCData(6:end,6:end)=NaN;
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
