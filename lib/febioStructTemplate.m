function [outputStruct]=febioStructTemplate(varargin)

% function [outputStruct]=febioStructTemplate(inputStruct)
% ------------------------------------------------------------------------
% This function creates a template for the febio_spec structure for FEBio.
% This structure is intended to be converted to an XML file using
% |febioStruct2xml|. 
%
%
% Change log: 
% 2020/12/01: KMM Changed to use febio_spec 3.0 by default
% 2021/04/29: KMM Removed strain energy output as this is not available for
% all materials (e.g. viscoelastic materials). 
% 2023/03/06: KMM Updating for FEBio spec 4.0. 
%
% To do: 
% * Check if still plot_stride is still needed
% * Check if both 0 and 1 and "symemtric" and presumably
% "assymetric" are options for the solver
% * No documentation for check_zero_diagonal
% 

% ------------------------------------------------------------------------

%%

switch nargin
    case 0 
        inputStruct=[];
    case 1
        inputStruct=varargin{1};
    otherwise
        error('Wrong number of input arguments');
end

%% Check for old FEBio version (febio2) 

if contains(lower(getFEBioPath),'febio2')
    [outputStruct]=febioStructTemplate_v2p5(inputStruct);     
elseif contains(lower(getFEBioPath),'febio3')
    [outputStruct]=febioStructTemplate_v3p0(inputStruct);     
    return    
end

%% Set febio_spec version
febio_spec.ATTR.version='4.0';

%% Module section
febio_spec.Module.ATTR.type='solid'; %Use default set
% febio_spec.Module.units='SI'; %Set unit system

%% Control section

febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=10;
febio_spec.Control.step_size=0.1;

febio_spec.Control.plot_zero_state=1; %Output initial state
febio_spec.Control.plot_range=[0,-1];
febio_spec.Control.plot_level='PLOT_MAJOR_ITRS';
febio_spec.Control.plot_stride=1; 
febio_spec.Control.output_level='OUTPUT_MAJOR_ITRS'; 

febio_spec.Control.adaptor_re_solve=1;

febio_spec.Control.time_stepper.ATTR.type='default';
febio_spec.Control.time_stepper.max_retries=5;
febio_spec.Control.time_stepper.opt_iter=6;
febio_spec.Control.time_stepper.dtmin=0.01;
febio_spec.Control.time_stepper.dtmax=0.1; 
febio_spec.Control.time_stepper.aggressiveness=0; 
febio_spec.Control.time_stepper.cutback=0.5; 
febio_spec.Control.time_stepper.dtforce=0; 

febio_spec.Control.solver.ATTR.type='solid';
febio_spec.Control.solver.symmetric_stiffness=1;%'symmetric';
febio_spec.Control.solver.equation_scheme='staggered';
febio_spec.Control.solver.equation_order='default';
febio_spec.Control.solver.optimize_bw=0;
febio_spec.Control.solver.lstol=0.9;

febio_spec.Control.solver.lsmin=0.01;
febio_spec.Control.solver.lsiter=5;

febio_spec.Control.solver.max_refs=25;
febio_spec.Control.solver.check_zero_diagonal=0;
febio_spec.Control.solver.zero_diagonal_tol=0;
febio_spec.Control.solver.force_partition=0;
febio_spec.Control.solver.reform_each_time_step=1;
febio_spec.Control.solver.reform_augment=0;
febio_spec.Control.solver.diverge_reform=1;
febio_spec.Control.solver.min_residual=1e-20;
febio_spec.Control.solver.max_residual=0;
febio_spec.Control.solver.dtol=0.001;
febio_spec.Control.solver.etol=0.01;
febio_spec.Control.solver.rtol=0;
febio_spec.Control.solver.rhoi=0;
febio_spec.Control.solver.alpha=1;
febio_spec.Control.solver.beta=0.25;
febio_spec.Control.solver.gamma=0.5;
febio_spec.Control.solver.logSolve=0;
febio_spec.Control.solver.arc_length=0;
febio_spec.Control.solver.arc_length_scale=0;

febio_spec.Control.solver.qn_method.ATTR.type='BFGS';
febio_spec.Control.solver.qn_method.max_ups=0; %Default 10
febio_spec.Control.solver.qn_method.max_buffer_size=0;
febio_spec.Control.solver.qn_method.cycle_buffer=1;
febio_spec.Control.solver.qn_method.cmax=1e5;

%% Globals section

febio_spec.Globals.Constants.R=8.314e-6; 
febio_spec.Globals.Constants.T=298;
febio_spec.Globals.Constants.Fc=96485e-9; 

%% Material section


%% LoadData section
% febio_spec.LoadData.load_controller{1}.ATTR.id=1;
% febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
% febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
% febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
% febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1];

%% Output section

% Plot file
febio_spec.Output.plotfile.ATTR.type='febio';
febio_spec.Output.plotfile.var{1}.ATTR.type='displacement';
febio_spec.Output.plotfile.var{2}.ATTR.type='stress';
febio_spec.Output.plotfile.var{3}.ATTR.type='relative volume';
febio_spec.Output.plotfile.var{4}.ATTR.type='reaction forces';
febio_spec.Output.plotfile.var{5}.ATTR.type='contact pressure';
febio_spec.Output.plotfile.compression=0;

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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
