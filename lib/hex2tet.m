function [TET,Vtet,C]=hex2tet(HEX,V,C,tetOpt)

C=C(:);
switch tetOpt
    case 1 %Add central node and cross side faces
        
        [F,~]=element2patch(HEX,C);
        
        numV=size(V,1);
        numE=size(HEX,1);
        
        %The original vertices
        X=V(:,1); Y=V(:,2); Z=V(:,3);
        
        %The mid-element points
        if numE==1;
            Vm=[mean(X(HEX),1) mean(Y(HEX),1) mean(Z(HEX),1)];
        else
            Vm=[mean(X(HEX),2) mean(Y(HEX),2) mean(Z(HEX),2)];
        end
        
        %The mid-face points
        Vf=[mean(X(F),2) mean(Y(F),2) mean(Z(F),2)];
        
        %TET point collection
        Vtet=[V; Vm; Vf];
        
        %Defining tetrahedral node list per hex element
        numV2=numV+numE+1;
        indAdd=reshape(numV2:(numV2+6*numE)-1,numE,6);        
        TET_set=[sort(HEX,2) (numV+1:numV+numE)' indAdd];
        
        %TET format based on Delaunay tesselation of a single set
        TET_format=[3,1,11,12;2,6,15,10;13,5,15,9;10,8,13,14;10,6,15,9;...
            6,5,15,13;1,2,15,12;2,10,12,4;6,13,15,9;10,8,14,4;13,8,7,14;...
            11,1,15,12;14,3,12,4;5,11,15,9;14,3,9,12;10,13,6,9;10,6,13,8;...
            13,11,5,9;2,9,15,12;2,10,9,12;10,9,12,4;9,14,12,4;10,14,9,4;...
            10,13,9,14;9,11,15,12;5,1,15,11;5,11,13,7;2,10,15,9;3,11,9,12;...
            9,11,3,14;11,7,3,14;13,7,11,14;13,11,9,14];        
        
        %Reform TET_set as an nx4
        TET_set_reform1=TET_set(:,TET_format(:))';        
        TET_set_reform2=reshape(TET_set_reform1,size(TET_format,1),numel(TET_set_reform1)/size(TET_format,1))';
        TET=reshape(TET_set_reform2,4,numel(TET_set_reform2)/4)';
        
        %Fix color information
        C=repmat(C,size(TET_format,1),1);        
        
    case 2 %Delaunay of single hex applied to all
        HEX=HEX(:,[1 2 4 3 5 6 8 7]);
        tetInd =[5     1     2     3;...
            6     5     2     3;...
            6     7     5     3;...
            6     4     7     3;...
            6     2     4     3;...
            6     8     7     4];
        a=tetInd';
        a=a(:)';
        A=HEX(:,a);
        TET=reshape(A',4,6.*size(HEX,1))';
        if ~isempty(C)
            C=(ones(6,1)*C');
            C=C(:);
        end
        Vtet=V;
    case 3 %Same as 2 flip top to bottom
        
        %Switch top and bottom
        HEX=HEX(:,[8:-1:5 4:-1:1]);
        
        HEX=HEX(:,[1 2 4 3 5 6 8 7]);
        tetInd =[5     1     2     3;...
            6     5     2     3;...
            6     7     5     3;...
            6     4     7     3;...
            6     2     4     3;...
            6     8     7     4];
        a=tetInd';
        a=a(:)';
        A=HEX(:,a);
        TET=reshape(A',4,6.*size(HEX,1))';
        if ~isempty(C)
            C=(ones(6,1)*C');
            C=C(:);
        end
        Vtet=V;
end

