function C=makeSymbolicStiffnessTensor(opt)

C=sym(zeros([3 3 3 3]));
I=eye(3,3); 
II1=dyadicProduct(I,I,1);
II3=dyadicProduct(I,I,3);

for i=1:3; 
    for j=1:3;
        for k=1:3; 
            for l=1:3;
                cvar=['c',num2str(i),num2str(j),num2str(k),num2str(l)]; 
         
                
                C(i,j,k,l)=sym(cvar); 
            
            end;
        end; 
    end; 
end;

switch opt
    case 'iso'
        L=(II1+II3)==0;
    case 'transiso'
        
    case'ortho'
        
    case'full'
        L=false(size(C));
    case 'empty'
        L=true(size(C));
end
C(L)=0;

end