function model_wrl = GB_read_vrml_fast_BA(filename)
% BA: I created this for single output as structure

MIMICS_10_KEY = '#coordinates written in 1mm / 10000';
MIMICS_13_KEY = '#coordinates written in 1mm / 0';

fid = fopen(filename, 'r');
if (fid == -1)
    fprintf(1, '\aError: unable to open input IV file %s.\n',filename);
    return,
end; 

input = fscanf(fid,'%c');
fclose(fid);

% so now lets find the middle points section
if ~isempty(regexpi(input,'(.*)texCoord(.*)')) % Bardiya Added this, in case the model has texture in it
    tokens = regexpi(input,'(.*)point\s+(.*)point\s+\[([^\]]*)\](.*)coordIndex\s+\[([\d\s\,\n\r\.-]*)\](.*)','tokens');
    %now lets get the points
    last_ind = strfind(tokens{1}{2},']');
    pts = sscanf(tokens{1}{2}(2:last_ind-1),'%f %f %f,',[3 inf])';
    connects = sscanf(tokens{1}{5},'%d, %d, %d, %d,',[4 inf])';
else
    tokens = regexpi(input,'(.*)point\s+\[([^\]]*)\](.*)coordIndex\s+\[([\d\s\,\n\r\.-]*)\](.*)','tokens');
    %now lets get the points
    pts = sscanf(tokens{1}{2},'%f %f %f,',[3 inf])';
    connects = sscanf(tokens{1}{4},'%d, %d, %d, %d,',[4 inf])';
end

% tokens = regexpi(input,'(.*)texCoord(.*)');



%check for the way mimics 10, outputs pts
if (isempty(pts)),
    tok = regexp(tokens{1}{2},'^[\s]*#[^\r\n]*(.*)','tokens');
    pts = sscanf(tok{1}{1},'%f %f %f,',[3 inf])';
end;

%check for alternate format, with no commas between the numbers, just at
%the end of the row
if (size(connects,1)==1 && size(connects,2)==1),
    connects = sscanf(tokens{1}{4},'%d %d %d %d',[4 inf])';    
end;

%% Look for scale
scale = [1 1 1];

%check in headers
for i=[1 3 5],
    scaleToken = regexpi(tokens{1}{i} ,'scale\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)','tokens');
    if (size(scaleToken,1)>0),
        scale = str2double(scaleToken{1,1});
    end;
end;

%if our scale is not exactly [1 1 1]
if (sum(scale~=1)>0),
    for i=1:3,
        pts(:,i) = pts(:,i)*scale(i);
    end;
end;

% check for the existance of the mimics 10 or 13 key, if so, the output was
% in m, and must be converted to mm
if (~isempty(strfind(input,MIMICS_10_KEY)) || ~isempty(strfind(input,MIMICS_13_KEY)))
    pts = pts .* 1000;
end;


%% BA Output
if abs(mean(mean(abs(pts)))) < 0.5
    disp('Converting Units: m to mm');
    pts = pts.*1e3;
end

model_wrl.pts = pts;
model_wrl.conn = connects;