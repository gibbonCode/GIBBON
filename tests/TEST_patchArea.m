classdef TEST_patchArea < matlab.unittest.TestCase
    properties (TestParameter)
        geom = testGeometries();
    end

    methods (Test)
        function testPatchArea(testCase, geom)
            areaEstimate = sum(patchArea(geom.F, geom.V), 'all');
            testCase.verifyEqual(areaEstimate, geom.areaTotalTrue, 'RelTol', 3.5e-3);
        end
    end
end

function C = testGeometries()

    C = struct();
    for idx = 1:6
    switch idx
    case 1 % Circle
        name = 'circle';
        r = 1;
        n = 250;
        t = linspace(-pi, pi, n)';
        t = t(1:end-1);
        V = [r*cos(t) r*sin(t) zeros(size(t))];
        F = 1:size(V,1);
        areaTotalTrue = pi*r^2;
    case 2 % Single element square 1x1
        name = 'singleElementSquare';
        V = [0 0 0; 1 0 0; 0 1 0; 1 1 0];
        F = [1 2 4 3];
        areaTotalTrue = 1;
    case 3 % Multi-element square rotated 1x1
        name = 'multiElementSquareRotated';
        V1 = 0.5*[-1 -1; -1 1; 1 1; 1 -1];
        regionCell = {V1};
        plotOn = 0;
        pointSpacing = 0.1;
        resampleCurveOpt = 1;
        [F,V] = regionTriMesh2D(regionCell, pointSpacing, resampleCurveOpt, plotOn);
        V(:,3) = 0;
        R = euler2DCM([-0.25*pi 0.25*pi 0]);
        V = V*R;
        areaTotalTrue = 1;
    case 4 % Cube
        name = 'cube';
        [V,F] = platonic_solid(2);
        V = V./max(V(:))/2;
        areaTotalTrue = 6;
    case 5 % Sphere quads
        name = 'sphereQuads';
        r = 1;
        n = 4;
        [F,V] = quadSphere(n, r, 2);
        areaTotalTrue = 4*pi*r^2;
    case 6 % Sphere triangles
        name = 'sphereTriangles';
        r = 1;
        n = 4;
        [F,V] = geoSphere(n, r);
        areaTotalTrue = 4*pi*r^2;
    end
    
    C.(name).F = F;
    C.(name).V = V;
    C.(name).areaTotalTrue = areaTotalTrue;
    end
end
