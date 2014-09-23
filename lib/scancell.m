function [IND_found]=scancell(T,targets,found_count)

line_count=1;
IND_found=cell(numel(targets),1);

if isempty(found_count)
    while 1
        for i=1:1:numel(targets)
            if ~isempty(strfind(T{line_count},targets{i}));
                IND_found{i}=[IND_found{i} line_count];
            end
        end
        line_count=line_count+1;
        if line_count>numel(T)
            break
        end
    end
else
    FOUND_count=zeros(size(found_count));
    while 1
        for i=1:1:numel(targets)
            if ~isempty(strfind(T{line_count},targets{i}));
                IND_found{i}=[IND_found{i} line_count];
                FOUND_count(i)=FOUND_count(i)+1;
            end
        end
        line_count=line_count+1;
        if line_count>numel(T)
            break
        end
        if all(FOUND_count==found_count)
            break
        end
    end
end

end