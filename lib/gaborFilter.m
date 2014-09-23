function GB=gaborFilter(Vr,S,f,tagShift,d,siz,p)

S=S(:)';
Sm=ones(size(Vr,1),1)*S;

%Gaussian component
G=exp( -0.5.*(sum((Vr./Sm).^2,2)) );
G=G./sum(G(:)); %Normalising

%Harmonic component
H=real(exp(1i.*(2*pi*f*(Vr(:,d)-tagShift))));

%Defining Gabor filter form
GB=G.*H;

%Normalising
% normFactor=abs(2.*sqrt(2).*exp(-2.*f.^2.*S(d).^2.*pi.^(3/2)).*prod(S).*cos(phaseOffset));
% GB=GB./normFactor;
% GB=GB./abs(sum(GB(:))); 

% Defining Gabor filter
GB=reshape(GB,siz); 

%Flip sign if required
GB=p.*GB;




