function Fmin=fmin_symbolic_form(P,OPT_struct,F)

%Getting function variables 
for i=1:1:numel(P)    
    eval(['P',num2str(i),'=P(i);']);
end
X=OPT_struct.FX;

%Evaluate function
Ff=subs(F);

%Derive optimisation function value
switch OPT_struct.method
    case 'fminsearch'
        Fmin=sum((Ff(:)-OPT_struct.FY(:)).^2);
    case 'lsqnonlin' 
        Fmin=(Ff-OPT_struct.FY);
end


end