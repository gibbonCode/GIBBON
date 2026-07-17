classdef TEST_patchVolume < matlab.unittest.TestCase
    properties (TestParameter)
        geom = testGeometries();
    end
    methods (Test)
        function testPatchVolume(testCase, geom)
            makeAbsolute = 0;
            volEst=patchVolume(geom.F, geom.V, makeAbsolute);
            testCase.verifyEqual(volEst, geom.volTotalTrue, 'RelTol', 3.5e-3);
        end
    end
end

function C = testGeometries()

    C = struct();
    for idx = 1:9
        switch idx   
        case 1 
            name = 'trianglulatedSphere';
            r=2;
            ns=5;
            [F,V]=geoSphere(ns,r);         
            volTotalTrue=4/3*pi*r^3; %True theoretical volume
        case 2 %Trianglulated sphere
            name = 'quadrangulatedSphere';
            r=2;
            ns=6;
            [F,V]=quadSphere(ns,r);         
            volTotalTrue=4/3*pi*r^3; %True theoretical volume

            %Shift randomly to create non-zero centred patch data
            V=V+10.*randn(1,3);
        case 3
            name = 'triangulatedTorus';
            r=1; %Sphere radius
            R=2; %Central radius
            nr=76;
            nc=150;
            [F,V]=patchTorus(r,nr,R,nc,'tri');
            volTotalTrue=2*pi^2*r*R; %True theoretical volume
        case 4
            name = 'quadrangulatedTorus';
            r=1; %Sphere radius
            R=2; %Central radius
            nr=76;
            nc=150;
            [F,V]=patchTorus(r,nr,R,nc,'quad');
            volTotalTrue=2*pi^2*r*R; %True theoretical volume
        case 5
            name = 'honeycombTorus';
            r=1; %Sphere radius
            R=2; %Central radius        
            nr=76;
            nc=150;
            [F,V]=patchTorus(r,nr,R,nc,'honey');
            volTotalTrue=2*pi^2*r*R; %True theoretical volume
        case 6
            name = 'triangulatedBox';
            boxDim=[5 2 1];
            pointSpacing=1;
            [F,V]=triBox(boxDim,pointSpacing); 
            volTotalTrue=prod(boxDim); %True theoretical volume
        case 7
            name = 'quadrangulatedBox';
            boxDim=[5 2 1];
            [F,V]=quadBox(boxDim,[5 2 1]);
            volTotalTrue=prod(boxDim); %True theoretical volume
        case 8
            name = 'mixedMeshOfPentahedraAndHexahedra';
            r=2;
            ns=5;
            [Ft,Vt]=geoSphere(ns,r);
            [V,F]=patch_dual(Vt,Ft);        
            volTotalTrue=4/3*pi*r^3; %True theoretical volume
        case 9
            name = 'trianglulatedSphereWithFlippedFaces';
            r=2;
            ns=5;
            [F,V]=geoSphere(ns,r);         
            F=fliplr(F);
            volTotalTrue=-4/3*pi*r^3; %True theoretical volume
        end
        C.(name).F = F;
        C.(name).V = V;
        C.(name).volTotalTrue = volTotalTrue;
    end
end