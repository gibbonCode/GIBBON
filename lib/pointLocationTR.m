function [TI,BC]=pointLocationTR(TR,QP,distCropOpt,chullCropOpt,waitbarOpt)

% function [TI,BC]=pointLocationTR(TR,QP,distCropOpt,chullCropOpt,waitbarOpt)
% ------------------------------------------------------------------------
%
% This function finds the triangles or tetrahedrons in which the query
% points QP are contained. In addition it outputs the barycentric
% coordinates.
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/04/01
%------------------------------------------------------------------------

E=TR.ConnectivityList;
numPerE=size(E,2);
numElements=size(E,1);

if distCropOpt==1
    [TR_centre_all,TR_rad_all] = circumcenter(TR,(1:numElements)');   
end

if chullCropOpt==1
    FBtri = freeBoundary(TR);
    DT=delaunayTriangulation(TR.Points(unique(FBtri(:)),:));
    ti_DT= pointLocation(DT,QP);
    
    logicDT=~isnan(ti_DT);
    QP=QP(logicDT,:);    
else
    logicDT=true(size(QP,1),1);
end

%% START LOOP

if waitbarOpt==1
    waitBarString='Finding points inside tesselation elements...';
    hw = waitbar(0,waitBarString);
    nUpdateSteps=25;
    waitbarSet=unique(round(linspace(1,numElements,nUpdateSteps)));
end

ti=nan(size(QP,1),1);
bc=nan(size(QP,1),numPerE);
for q=1:1:numElements
    
    if distCropOpt==1
        %Get current circum centre and radius
        TR_centre=TR_centre_all(q,:);        
        TR_rad_max=TR_rad_all(q);
        
        %Check distance with current points        
        QP_min=QP-TR_centre(ones(size(QP,1),1),:); %Query points with current element centre subtracted
        QP_dist=sqrt(sum(QP_min.^2,2)); %Distances of query points to current element centre        
        logicClose=QP_dist<=(TR_rad_max+eps(TR_rad_max)); 
        bcClose=nan(size(logicClose,1),numPerE);
        QP_test=QP(logicClose,:);
        if nnz(logicClose)>0
            [logicClose_sub,bcClose_sub]=isInsideTR(TR,QP_test,q*ones(size(QP_test,1),1));
            bcClose(logicClose,:)=bcClose_sub;            
            logicClose(logicClose)=logicClose_sub;        
        end             
    else
        [logicClose,bcClose]=isInsideTR(TR,QP,q*ones(size(QP,1),1));        
    end
    
    ti(logicClose)=q;
    
    bc(logicClose,:)=bcClose(logicClose,:);
    
    if waitbarOpt==1 && ismember(q,waitbarSet)
        waitbar(q/numElements,hw,waitBarString);
    end
end

if waitbarOpt==1
    close(hw);
end

TI=nan(size(logicDT));
TI(logicDT)=ti;

BC=nan(size(logicDT,1),numPerE);
BC(logicDT,:)=bc;

end