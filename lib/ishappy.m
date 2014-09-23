function L=ishappy(n)

N=n;
L=0;
NHistory=n;
while L==0
    N=sum(str2num(num2str(N).').^2); %Sum digits 
    if N==1 %n is a happy number
        L=1;
    elseif ismember(N,NHistory); %Kill we are in a loop
        break
    end
    NHistory=[NHistory N]; 
end
NHistory=[NHistory N]; 
