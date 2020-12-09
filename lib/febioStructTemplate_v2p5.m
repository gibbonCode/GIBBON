function [outputStruct]=febioStructTemplate_v2p5(varargin)

switch nargin
    case 0 
        inputStruct=[];
    case 1
        inputStruct=varargin{1};
    otherwise
        error('Wrong number of input arguments');
end

%% Set febio_spec version
febio_spec.ATTR.version='2.5';

%% Module section
febio_spec.Module.ATTR.type='solid'; %Use default set

%% Control section

febio_spec.Control.analysis.ATTR.type='static';
% febio_spec.Control.title='temp';
febio_spec.Control.time_steps=10;
febio_spec.Control.step_size=0.1;
febio_spec.Control.dtol=0.001;
febio_spec.Control.etol=0.01;
febio_spec.Control.rtol=0;
febio_spec.Control.lstol=0.9;
febio_spec.Control.time_stepper.dtmin=febio_spec.Control.step_size/3;
febio_spec.Control.time_stepper.dtmax=febio_spec.Control.step_size*3; %Constant
% febio_spec.Control.time_stepper.dtmax.ATTR.lc=1; %Load curve based
% febio_spec.Control.time_stepper.dtmax.VAL=0.1; %Load curve based
febio_spec.Control.time_stepper.max_retries=5;
febio_spec.Control.time_stepper.opt_iter=10;
febio_spec.Control.max_refs=15;
febio_spec.Control.max_ups=10;
febio_spec.Control.optimize_bw=0;
% febio_spec.Control.restart.ATTR.file='restartDumpFile.dmp'; %On with file specified
% febio_spec.Control.restart.VAL=0; %Off
febio_spec.Control.plot_level='PLOT_MAJOR_ITRS';
febio_spec.Control.plot_range=[0,-1];
febio_spec.Control.plot_stride=1;
febio_spec.Control.plot_zero_state=0;
febio_spec.Control.cmax=1e5;
febio_spec.Control.print_level='PRINT_MINOR_ITRS';
febio_spec.Control.min_residual=1e-20;
% febio_spec.Control.integration='N/A';
febio_spec.Control.output_level='OUTPUT_MAJOR_ITRS'; 

%% Globals section

febio_spec.Globals.Constants.T=0; 
febio_spec.Globals.Constants.R=0;
febio_spec.Globals.Constants.Fc=0; 

%% Material section
% febio_spec.Material.material{1}.ATTR.type='Ogden';
% febio_spec.Material.material{1}.ATTR.id=1;
% febio_spec.Material.material{1}.c1=1e-3;
% febio_spec.Material.material{1}.m1=6;
% febio_spec.Material.material{1}.c2=1e-3;
% febio_spec.Material.material{1}.m2=-6;
% febio_spec.Material.material{1}.k=1e-1;

%% LoadData section

febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; 1 1];

%% Output section

% Plot file
febio_spec.Output.plotfile.ATTR.type='febio';
febio_spec.Output.plotfile.var{1}.ATTR.type='displacement';
febio_spec.Output.plotfile.var{2}.ATTR.type='stress';
febio_spec.Output.plotfile.var{3}.ATTR.type='strain energy density';
febio_spec.Output.plotfile.var{4}.ATTR.type='relative volume';
febio_spec.Output.plotfile.var{5}.ATTR.type='reaction forces';
febio_spec.Output.plotfile.var{6}.ATTR.type='contact pressure';

%% Fill in missing if input structure is provided

switch nargin
    case 0 %Make output the full default 
        outputStruct=febio_spec;
    case 1 %Make output by complementing input with default
        [outputStruct]=structComplete(inputStruct,febio_spec,1); %Complement provided with default if missing
end
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.