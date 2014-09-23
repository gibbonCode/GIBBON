function cell2txtfile(filename,T,skipOpt)

fid=fopen(filename,'w');
for q=1:size(T,1);
    l=T{q,:};
    %     if ~ischar(l) %If this isn't a string then assume its a number
    %         l=num2str(l);
    %     end
    
    if skipOpt==0 || ~isempty(deblank(l))
        fprintf(fid,'%s\r\n',l);
    end
    
end
fclose(fid);
