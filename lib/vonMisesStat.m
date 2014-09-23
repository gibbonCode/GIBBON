function [MU,k]=vonMisesStat(T)

MU=angle(mean(exp(1i.*T)));
Rsq=mean(cos(T)).^2+mean(sin(T)).^2;
R=sqrt(Rsq);

p=2; 
ki=R.*(p-Rsq)./(1-Rsq);
diffTol=ki/1000;
qIter=1;
while 1;
    kn=ki;
    A=besseli(p/2,ki)./ besseli((p/2)-1,ki);
    ki=ki-((A-R)./(1-A^2-((p-1)/ki)*A));                
    if abs(kn-ki)<=diffTol
        break
    end
    qIter=qIter+1
end
k=ki;

