function [R,Ri]=euler2DCM(E)

% ----------------------------------------------------------------------
% function [R,Ri]=euler2DCM(E)
% 
% This function generates a rotation matrices based on the Euler angles in
% E.
%
% Kevin Mattheus Moerman, 01/06/2011
% kevinmoerman@hotmail.com
% ----------------------------------------------------------------------


switch class(E)
    case 'double'
        R=zeros(3,3,size(E,1));
    case 'sym'
        R=sym(zeros(3,3,size(E,1)));
end
Ri=R;

for i=1:1:size(E,1);
    
    Rx=[1        0        0;...
        0        cos(E(i,1))  -sin(E(i,1));...
        0        sin(E(i,1))  cos(E(i,1))];
    
    Ry=[cos(E(i,2))  0        sin(E(i,2));...
        0        1        0;...
        -sin(E(i,2)) 0        cos(E(i,2))];
    
    Rz=[cos(E(i,3))  -sin(E(i,3)) 0;...
        sin(E(i,3))  cos(E(i,3))  0;...
        0        0        1];
    
    Rxyz=Rx*Ry*Rz;
    R(:,:,i)=Rxyz;
    Ri(:,:,i)=inv(Rxyz);
    
end

end


