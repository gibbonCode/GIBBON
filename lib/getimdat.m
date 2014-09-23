function [D]=getimdat(file_name_IMDAT,opt,d)

D=load(file_name_IMDAT);
IMDAT=D.IMDAT;
siz=IMDAT.dat;
siz_par=IMDAT.par;

switch opt
    case 'dat'
        D=zeros([siz(1) siz(2) siz(3) numel(d)],'uint16');
    case 'par'
        D=repmat(IMDAT.M_info,[siz_par(1) numel(d)]);
end

switch IMDAT.type
    case 'full'
        disp('WARNING! Upload not done. See IMDAT.dat for data');
    case 'split'
        switch opt
            case 'dat'
                for i=1:1:numel(d)
                    IMDAT.load_names{d(i)}
                    load(IMDAT.load_names{d(i)});
                    D(:,:,:,i)=m;
                end
            case 'par'
                for i=1:1:numel(d)
                    load(IMDAT.load_names_par{d(i)});
                    D(:,i)=p;
                end
        end
end

