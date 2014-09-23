function d=divisors(x,arg)

% Determine divisors of input x (which needs to be natural number)
% By default the output provides a vector containing all positive divisors
% but this can be altered depending on arg. 


%Default if no option is given
if exist('arg')==0
    arg='pos';
end

x=abs(x);
X=0:1:x;

d=X(ismember(X,abs(x)./X));
switch arg
    case 'pos'
    case 'neg'
        d=-fliplr(d);
    case 'all'
        d=[-fliplr(d) d];
    otherwise
        error('False input!');
end

end
