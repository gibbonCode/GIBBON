function [CM]=fourthOrderCell(C)

CM=cell(3,3);

for i=1:1:3
    for j=1:1:3
        switch class(C)
            case 'sym'
                C_sub=sym(zeros(3,3));
            otherwise
                C_sub=zeros(3,3);                
        end
        
        for k=1:1:3
            for l=1:1:3                
                C_sub(k,l)=C(i,j,k,l);
            end
        end
        CM{i,j}=C_sub;
    end
end