function c=voigtUnMap(cVoigt)


if isvector(cVoigt) %assume that c is a 4th order tensor
    siz_c=[3 3];
    cVoigt(4:end)=(1/2).*cVoigt(4:end); %Undo doubling
else
    siz_c=[3 3 3 3];
end

[linearIndexVoigt,linearIndexFourthOrder]=tensor2voigtMap(zeros(siz_c));
    
switch class(cVoigt)
    case 'double'
        c=zeros(siz_c);
    case 'sym'
        c=sym(zeros(siz_c));
end

c(linearIndexFourthOrder)=cVoigt(linearIndexVoigt);


end

