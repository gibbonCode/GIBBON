function cmap=ccmap(num,n)

switch num
    case 1 %jet
        cmap=jet(n);
    case 2 %rgb
        red_part=linspace(-1,1,n)';
        red_part(red_part<0)=0;
        green_part=flipud(red_part);
        x=linspace(0,1,n)';
        blue_part=x;
        blue_part(x<=0.5)=2.*x(x<=0.5);
        blue_part(x>=0.5)=2.*(1-x(x>=0.5));
        cmap=[red_part(:) green_part(:) blue_part(:)];
    case 3 %rwb
        red_part=linspace(-1,1,n)';
        red_part(red_part<0)=0;
        blue_part=flipud(red_part);
        x=linspace(0,1,n)';
        white_part=x;
        white_part(x<=0.5)=2.*x(x<=0.5);
        white_part(x>=0.5)=2.*(1-x(x>=0.5));
        white_part=white_part*ones(1,3);
        
        cmap=zeros(n,3);
        cmap(:,1)=red_part(:);
        cmap(:,3)=blue_part(:);
        cmap=cmap+white_part;
    case 4 %rkb
        red_part=linspace(-1,1,n)';
        red_part(red_part<0)=0;
        blue_part=flipud(red_part);
        cmap=0.1.*ones(n,3);
        cmap(:,1)=cmap(:,1)+0.9.*red_part(:);
        cmap(:,3)=cmap(:,3)+0.9.*blue_part(:);
    case 5 %ryg
        red_part=linspace(-1,1,n)';
        red_part(red_part<0)=0;
        green_part=flipud(red_part);
        x=linspace(0,1,n)';
        yellow_part=x;
        yellow_part(x<=0.5)=2.*x(x<=0.5);
        yellow_part(x>=0.5)=2.*(1-x(x>=0.5));
        cmap=zeros(n,3);
        cmap(:,1)=red_part(:)+yellow_part;
        cmap(:,2)=green_part(:)+yellow_part;
    case 6 %ryg classic
        red_part=linspace(-1,1,n)';
        red_part(red_part<0)=0;
        green_part=flipud(red_part);
        x=linspace(0,1,n)';
        yellow_part=x;
        yellow_part(x<=0.5)=2.*x(x<=0.5);
        yellow_part(x>=0.5)=2.*(1-x(x>=0.5));
        cmap=zeros(n,3);
        cmap(:,1)=0.7.*red_part(:)+yellow_part;
        cmap(:,2)=0.2.*green_part(:)+yellow_part;
    case 7 %rkg classic
        red_part=linspace(-1,1,n)';
        red_part(red_part<0)=0;
        green_part=flipud(red_part);
        cmap=0.1.*ones(n,3);
        cmap(:,1)=cmap(:,1)+0.7.*red_part(:);
        cmap(:,2)=cmap(:,2)+0.3.*green_part(:);
end

end