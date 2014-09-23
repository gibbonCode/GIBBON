function [F,C]=element2faces(E,C)

%Converts elements to faces enabling PATCH based visualisation colordata is
%copied for each face. Double faces for bordering elements may occur. 
%
%
%19/01/2012, Kevin Mattheus Moerman

switch size(E,2)
    case 4 %4 node tetrahedral elements
        F=[E(:,[1 2 3]);... 
           E(:,[1 2 4]);... 
           E(:,[2 3 4]);... 
           E(:,[3 1 4])]; 
       C=repmat(C(:),4,1);
    case 8 %8 node hexahedral elements
        F=[E(:,[1 2 3 4]);... %top
           E(:,[5 6 7 8]);... %bottom
           E(:,[1 2 6 5]);... %side 1
           E(:,[3 4 8 7]);... %side 2
           E(:,[2 3 7 6]);... %front
           E(:,[1 4 8 5]);]; %back       
       C=repmat(C(:),6,1);
    otherwise
        error('MATLAB:ELEMENT2FACES: size(E,1) not consistent with tetrahedral or hexahedral element');
end

end