function FgradVisualize(Fbig,boxDim,F,V,plt,ang)

    %Creates an interactive visualization of a deformation or a series of
    %deformations, on a single element geometry. Final output plot has an
    %animation slider, and can be exported as a gif for external use. 
    
    %----MANDATORY ARGUMENTS----%
    %Fbig: 1xn Cell array containing n deformation gradients (3x3 matrices)
    
    %----OPTIONAL ARGUMENTS (leave as [] if not needed)----:
    %- Fibre angle (ang): A 1x3 or 3x1 matrix of directional cosines corresponding
    %to the element fibre direction. Include to have that visualized in the
    %plot
    
    %- Plot Handle (plt): A plot handle that the figure can optionally be
    %plotted on. Use to generate the interactive figure on an existing
    %figure/subplot. If left as [], interactive figure will be generated on a fresh
    %plot.
    
    %Geometry: Two options available:
    %  - If you have custom geometry, provide patch faces (F) and nodes
    %  (V).
    %  - If not, a generic rectangular element can be used. Provide
    %  dimensions of the box (boxDim) as a 1x3 vector [L, W, T].     
    %  - If none of the above is provided, a single-element 1x1 cube will be
    %  used.
    
    %NOTE: If using F,V, leave boxDim = [], and if using boxDim leave F & V
    %as []. If both are provided, F,V will be used.
     
    %Made by Darshan Senthil, 2026
        
    boxFlag=0; FVflag=0;
    %Determining which geometry to use.
    if(~isempty(F))
        if(isempty(V) && ~isempty(boxDim))
            warning('Entries for faces detected but no nodes given - switching to provided box dimensions');
            boxFlag=1;
        elseif(isempty(boxDim) && isempty(V))
            warning('Entries for faces detected but no nodes given - switching to single-element cube');
%             error("Input arguments specified incorrectly. Either both F & V or boxDim must be specified")
        elseif(~isempty(V))
            FVflag=1;            
        end
        
    elseif(isempty(F))
        if(isempty(boxDim))
%             error("Input arguments specified incorrectly. Either both F & V or boxDim must be specified")
        else
            boxFlag=1;
        end
    end
        
    
    %Loading geometry
    if(FVflag==1)
        Ff=F; 
        V=V;
    elseif(boxFlag==1 && FVflag==0)
        [meshStruct]=hexMeshBox(boxDim,[1,1,1]);
        Ff = meshStruct.F;
        V = meshStruct.V;
    else
%         error("Input arguments specified incorrectly. Either both F & V or boxDim must be specified")
        [meshStruct]=hexMeshBox([1,1,1],[1,1,1]);
        Ff = meshStruct.F;
        V = meshStruct.V;

    end

    centroid = mean(V);
    [n1mat, n2mat, n3mat, Vtrans] = applyFgrads(Fbig,V,centroid);

    if ~isempty(plt)
        hf = plt;
        figure(hf); hold on
    else
        hf = figure; hold on;
    end
    gpatch(Ff,V,'y','k',0.3);
    hpe=gpatch(Ff,V,'g','k',0.4);
    % hpf=quiverVec(centroid,ang,2,'r','k');
    if(~isempty(ang))
        hpf=quiver3(centroid(1,1),centroid(1,2),centroid(1,3),ang(1),ang(2),ang(3),2.5,'LineWidth',3,'Color','r','DisplayName','Fibre direction \(f)');
%         lgd=legend({'Deformed Element','Fibre Direction','n11','n22','n33'});
    else
%         lgd=legend({'Deformed Element','n11','n22','n33'});
    end
    % hp1=quiverVec(centroid,n1mat{1}',0.75,'m','m');
    % hp2=quiverVec(centroid,n2mat{1}',0.75,'g','g');
    % hp3=quiverVec(centroid,n3mat{1}',0.75,'b','b');
    hp1=quiver3(centroid(1,1),centroid(1,2),centroid(1,3),n1mat{1}(1),n1mat{1}(2),n1mat{1}(3),1.5,'LineWidth',3,'Color','m','DisplayName','n11');
    hp2=quiver3(centroid(1,1),centroid(1,2),centroid(1,3),n2mat{1}(1),n2mat{1}(2),n2mat{1}(3),1.5,'LineWidth',3,'Color','g','DisplayName','n22');
    hp3=quiver3(centroid(1,1),centroid(1,2),centroid(1,3),n3mat{1}(1),n3mat{1}(2),n3mat{1}(3),1.5,'LineWidth',3,'Color','b','DisplayName','n33');
    % hp1=quiver3(centroid(1,:),centroid(2,:),centroid(3,:),n1mat{1}(1,:),n1mat{2}(1,:),n1mat{3}(1,:),2);
    % title('Principal Directions Through Deformation')
    xlabel('x'); ylabel('y'); zlabel('z');
    axisGeom;

%     lgd=legend({'Element','Deformed Element','Fibre Direction','n11','n22','n33'});
    if(~isempty(ang))
        hpf=quiver3(centroid(1,1),centroid(1,2),centroid(1,3),ang(1),ang(2),ang(3),2.5,'LineWidth',3,'Color','r','DisplayName','Fibre direction \(f)');
         lgd=legend({'Element','Deformed Element','Fibre Direction','n11','n22','n33'});
    else
         lgd=legend({'Element','Deformed Element','n11','n22','n33'});
    end

    hpJ=plot3(centroid(1),centroid(2),centroid(3),'ko','DisplayName','J=1');

    for q=1:1:length(Vtrans)
        n1=n1mat{q}; 
        n2=n2mat{q}; 
        n3=n3mat{q};

        Jstr=strcat('J=',num2str(det(Fbig{q})));

        %Set entries in animation structure
        animStruct.Handles{q}=[hp1,hp1,hp1,hp2,hp2,hp2,hp3,hp3,hp3,hpe,hpJ]; %Handles of objects to animate
        animStruct.Props{q}={'UData','VData','WData','UData','VData','WData','UData','VData','WData','Vertices','DisplayName'}; %Properties of objects to animate
        animStruct.Set{q}={n1(1),n1(2),n1(3),n2(1),n2(2),n2(3),n3(1),n3(2),n3(3),Vtrans{q},Jstr}; %Property values for to set in order to animate

    end
    anim8(hf,animStruct);
    r = rotate3d(hf);
end


function [n1mat, n2mat, n3mat, Vtrans] = applyFgrads(Fbig,V,centroid)
%Applies each deformation gradient to the geometry, and calculates the
%resulting principal directions. Might also calculate fibre stretch so the
%fibre vector can scale accordingly but that might not be worth the effort
    for i=1:length(Fbig)
        F=Fbig{i};
        B=F'*F; J=det(F);
        [nn,lam]=eig(B);

        l3temp=lam(3,3);
        lam(3,3)=lam(1,1); lam(1,1)=l3temp;
        nn(:,[1,3]) = fliplr(nn(:,[1,3]));

        lam1(i)=(J^(-1/3))*sqrt(lam(1,1));
        lam2(i)=(J^(-1/3))*sqrt(lam(2,2));
        lam3(i)=(J^(-1/3))*sqrt(lam(3,3));

        n1=nn(:,1);
        n2=nn(:,2);
        n3=nn(:,3);

        n1mat{i}=n1; n2mat{i}=n2; n3mat{i}=n3;
        for x=1:length(V)
            Vtrans{i}(x,:)=(F*(V(x,:)-centroid)')'+centroid;
        end

    end


end