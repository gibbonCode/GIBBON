function P=bezierCurve(p,n)

%%

% if size(p,2)==2
%     p(:,3)=0;
% end

%%

t=linspace(0,1,n)';
N=size(p,1);
nn=0:1:N-1;
f=factorial(nn); 
s=factorial(N-1)./( f.*flip(f) ); %Sigma
P=( s.* ((1-t).^flip(nn)) .* (t.^nn) )*p; %Output
