function [P]=triIncenter(F,V)

% function [P]=triIncenter(F,V)
% ------------------------------------------------------------------------
% Computes the triangle incenter for the triangles defined by the face
% array F and the vertex array V. 
% 
% Change log: 
% 2022/04/28 Created
% ------------------------------------------------------------------------

%%

try %Try to use triangulation method
    warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
    TR=triangulation(F,V);
    P=incenter(TR);
    warning('on','MATLAB:triangulation:PtsNotInTriWarnId');
catch %Use custom approach instead
    a=sqrt( sum( (V(F(:,2),:)-V(F(:,3),:)).^2 ,2) );
    b=sqrt( sum( (V(F(:,1),:)-V(F(:,3),:)).^2 ,2) );
    c=sqrt( sum( (V(F(:,1),:)-V(F(:,2),:)).^2 ,2) );

    P=zeros(size(F,1),size(V,2));
    ABC=[a b c];
    sABC=sum(ABC,2);
    for q=1:1:size(V,2)
        X=V(:,q);
        P(:,q)=sum(X(F).*ABC,2)./sABC;
    end
end