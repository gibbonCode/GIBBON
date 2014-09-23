function save_figure_image(hf,filename,r,imtype,opt)

set(hf,'InvertHardcopy','off');
switch opt
    case 1        
        print (imtype, ['-r',num2str(r)], filename);
    case 2
        I=getframe(hf);
        [I,Map] = frame2im(I);
        imwrite(I, filename,'Compression','none','Resolution',r);
end

end


