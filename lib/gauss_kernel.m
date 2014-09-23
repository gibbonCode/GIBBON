function hg=gauss_kernel(k,nd,f,m)

switch m
    case 'sigma'
        %NOTE: This method is equivalent to using hg = fspecial('gaussian',[k k], S)
        S=f;
        switch nd
            case 1
                x = linspace(-((k-1)/2),((k-1)/2),k);
                hg=exp(-(x.^2)./(2*S^2));
            case 2
                [x,y] = meshgrid(linspace(-((k-1)/2),((k-1)/2),k));
                hg=exp(-(x.^2 + y.^2)./(2*S^2));
            case 3
                [x,y,z] = meshgrid(linspace(-((k-1)/2),((k-1)/2),k));
                hg=exp(-(x.^2 + y.^2 + z.^2)./(2*S^2));
        end
    case 'width'
        % Here "f" defines "where the bell curve is" at the edges of the
        % kernel e.g. f=2 results in twice the standard deviation.
        switch nd
            case 1
                x = linspace(-f,f,k);
                hg=exp(-(x.^2)./2);
            case 2
                [x,y] = meshgrid(linspace(-f,f,k));
                hg=exp(-(x.^2 + y.^2)./2);
            case 3
                [x,y,z] = meshgrid(linspace(-f,f,k));
                hg=exp(-(x.^2 + y.^2 + z.^2)./2);
        end
    otherwise
        warning('False input for method! Valid options are "sigma" and "width"');
end

hg(hg<eps*max(hg(:)))=0;

%Making sure mask sums to 1
hg=hg./sum(hg(:));


end