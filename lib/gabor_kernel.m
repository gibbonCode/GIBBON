function GB=gabor_kernel(Vr,S,f,d,siz,p)

%Gaussian component
G=exp( -0.5.*((Vr(:,1)./S(1)).^2+(Vr(:,2)./S(2)).^2+(Vr(:,3)./S(3)).^2) );

%Harmonic component
H=p.*cos(2*pi*f*Vr(:,d));

% Defining Gabor filter
GB=reshape(G.*H,siz); 


% % Gaussian component
% G=exp( -0.5.*((Vr(:,1)./S(1)).^2+(Vr(:,2)./S(2)).^2+(Vr(:,3)./S(3)).^2) );
% G=G./sum(G(:));  %Normalise to sum to 1
% 
% %Harmonic component
% H=p.*cos(2*pi*f*Vr(:,d));
% 
% % Defining Gabor filter
% GB=reshape(G.*H,siz); 
% GB=GB./sum(abs(GB(:))); %Normalise such that sum(abs(GB))==1
% % GB=GB./abs(sum(GB(:))); %Normalise such that sum(GB)==1



