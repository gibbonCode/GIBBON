function cVoigt=voigtMap(c)

siz2=3*ones(1,2);
siz4=3*ones(1,4);

[linearIndexVoigt,linearIndexFourthOrder]=tensor2voigtMap(c);
secondOrder=0;
if ndims(c)==4 %4d array e.g. 4th order tensor
    if all(size(c)==siz4) %4th order tensor
        switch class(c)
            case 'double'
                cVoigt=zeros(6,6);
            case 'sym'
                cVoigt=sym(zeros(6,6));
        end
    else
        %TO DO add warning here
    end
elseif all(size(c)==siz2) %2nd order tensor
    secondOrder=1;
        switch class(c)
            case 'double'
                cVoigt=zeros(6,1);
            case 'sym'
                cVoigt=sym(zeros(6,1));
        end
else
    %TO DO add warning here
end

cVoigt(linearIndexVoigt)=c(linearIndexFourthOrder);
if secondOrder==1
    cVoigt(4:end)=2.*cVoigt(4:end);
end

end

