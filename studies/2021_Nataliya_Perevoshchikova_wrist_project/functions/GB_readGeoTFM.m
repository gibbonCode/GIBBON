function [RT] = GB_readGeoTFM(transform_fname)
% function [RT] = loadGeomagicTransform(transform_fname)
% tries to read the transform file, and account for a scale if present, so
% that it converts m -> mm

% transform_fname = 'P:\Bob\Primate\E85281\S01R\rad15R01R.tfm';

fid = fopen(transform_fname,'r');
data = fscanf(fid,'%f',[4 4])';
unit = fscanf(fid,'Units: %s');
fclose(fid);


RT = data;
if (strcmp(unit,'mm'))
    %do nothing, its what we want
elseif (strcmp(unit,'m'))
    RT(1:3,4) = RT(1:3,4)*1000;
end


