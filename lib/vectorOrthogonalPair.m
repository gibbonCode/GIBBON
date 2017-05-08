function [a,d]=vectorOrthogonalPair(f)

% function [a,d]=vectorOrthogonalPair(f)
% ------------------------------------------------------------------------
%
% Based on the input vector f this function generates the normalized output
% vectors a and d which are orthogonal to eachother and to f. 
% 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2017/05/08: Updated for GIBBON
% 2017/05/08: Fixed bug in relation to co-linear output.

%%

f=vecnormalize(f);

e1=[1 0 0];
e2=[0 1 0];
e3=[0 0 1];

E1=e1(ones(size(f,1),1),:);
E2=e2(ones(size(f,1),1),:);
E3=e3(ones(size(f,1),1),:);
 
% T=abs(real([acos(dot(f,E1,2)) acos(dot(f,E2,2)) acos(dot(f,E3,2))]));
% T(T>pi)=(T(T>pi)-pi);
% T=abs(T-0.5*pi);
% [~,J_min]=min(T,[],2);

[~,J_min]=min(abs(f),[],2);

d=zeros(size(f));
for q=1:1:3
   logicNow=J_min==q;   
   if any(logicNow)
       switch q
           case 1
               E=E1(logicNow,:);
           case 2
               E=E2(logicNow,:);
           case 3
               E=E3(logicNow,:);
       end
       d(logicNow,:)=cross(E,f(logicNow,:),2);
   end
end
[d]=vecnormalize(d);
a=cross(f,d,2); [a]=vecnormalize(a); %d is orthogonal to f and a
d=cross(f,a,2); [d]=vecnormalize(d); %a is reset to be orthogonal to both f and d
