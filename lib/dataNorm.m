function Mn=dataNorm(M,normOpt)

switch normOpt
    case 1 %Using max and min
        Mn=double(M);
        Mn=Mn-min(Mn(:));
        Mn=Mn./max(Mn(:));
    case 2 %Using max only and letting M(M<0)= 0
        Mn=double(M);
        Mn(Mn<0)=0;        
        Mn=Mn./max(Mn(:));
end