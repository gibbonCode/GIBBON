function Fmin=fmin_symbolic_form(Pn,OPT_struct,F)

% %Getting function variables 
% for i=1:1:numel(P)    
%     eval(['P',num2str(i),'=P(i);']);
% end

X=OPT_struct.FX;

%Evaluate function
P=sym('P',size(Pn)); %Create P as symbolic
Ff=subs(subs(F,P,Pn)); %Substitution of P and other

%Derive optimisation function value
switch OPT_struct.method
    case 'fminsearch'
        Fmin=sum((Ff(:)-OPT_struct.FY(:)).^2);
    case 'lsqnonlin' 
        Fmin=(Ff-OPT_struct.FY);
end
Fmin=double(Fmin);

end