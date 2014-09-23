function [subIndOut,L_valid]=snapSubInd(subIndIn,siz)

%Initialize outputs
L_valid=false(size(subIndIn));
subIndOut=subIndIn;

%Fixing out of range indices
for q=1:1:size(subIndIn,2)
    subInd=subIndIn(:,q); %The current subscript index set
    L_valid(:,q)=(subInd>0) & (subInd<=siz(q)); %Logic for valid indices not out of range
    
    subIndToFix=subInd(~L_valid(:,q)); %Indices to fix
    
    %Snapping out of range indices to boundary
    subIndToFix=(subIndToFix.*(subIndToFix>1))+(subIndToFix<1); %Fix smaller than 1
    subIndToFix=(subIndToFix.*(subIndToFix<=siz(q)))+siz(q).*(subIndToFix>siz(q)); %Fix larger than siz(q)

    %Storing fixed indices in output
    subIndOut_current=subInd;
    subIndOut_current(~L_valid(:,q))=subIndToFix;
    subIndOut(:,q)=subIndOut_current;
end

