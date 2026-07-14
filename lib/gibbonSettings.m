classdef (Abstract) gibbonSettings
% Static class to hold Gibbon settings,
% meant to be a single-point interface to set/get global GIBBON options.
%
% It currently uses MATLAB's `settings` for persistent storage,
% but the "back-end" could eventually be moved e.g. to a JSON file, for
% compatibility with Octave.
%
% Example:
%   s = gibbonSettings.get()
%
%   gibbonSettings.set('view','touch')
%   gibbonSettings.get('view')
%
%   gibbonSettings.reset()
%   gibbonSettings.get()
%
%   gibbonSettings.set(s)
%
% See also: settings

    properties
        FEBioPath (1,:) char {mustBeFileOrEmpty, warnIfFEBioVersion2} = ''
        ViewProfile (1,:) char = 'CAD'
        % FutureSetting (size) class {validators} = defaultValue
    end

    properties (Constant)
        gibbonPath = fileparts(fileparts(mfilename('fullpath')));

        % Used by pre-set validatestring in gibbonSettings.set
        validViewProfiles = {'CAD', 'febio', 'touchpad'}
    end

    properties (Constant, Hidden)
    % Option names, defaults and validators are picked up from the 
    % properties (and property attributes) above.
        names (1,:) cell = getNames()
        defaults (1,1) struct = getDefaults()
        validators (1,1) struct = getValidators()
    end

    methods(Static)

        function value = get(opt)
        % Get option(s) from settings.GIBBON, fall back to default(s)
        % Syntax:
        %   gibbonSettings.get() - get structure of settings
        %   gibbonSettings.get(NAME) - get a particular setting

            if nargin == 0
                opts = gibbonSettings.names;
                vals = cellfun(@gibbonSettings.get, opts, 'unif', 0);
                value = cell2struct(vals', opts');
                return
            end

            opt = validatestring(opt, gibbonSettings.names);
            value = gibbonSettings.defaults.(opt);
            
            s = settings;
            if hasGroup(s,'GIBBON') && hasSetting(s.GIBBON,opt)
                v = s.GIBBON.(opt).PersonalValue;
                if ~isempty(v)
                    value = v; 
                end
            end
        end

        function set(varargin)
        % Write option(s) to settings.GIBBON
        % Syntax:
        %   gibbonSettings.set(Name, Value, ...)
        %
        % Use gibbonSettings.get() to see current settings and available
        % names. Case-insensitive [short-hand] names are allowed.
        % 
        % Examples:
        %   gibbonSettings.set('febio','febio4');  % set FEBioPath
        %   gibbonSettings.set(ViewProfile='cad'); % set ViewProfile='CAD'
        %
        % See also: gibbonSettings.get, gibbonSettings.reset

            if nargin == 1 && isstruct(varargin{1})
                names = fieldnames(varargin{1});
                values = struct2cell(varargin{1});
            else
                assert(nargin > 0 && mod(nargin,2) == 0, ...
                    'gibbon:Settings:SetArgs','Expecting one or more name-value pairs')
                names = varargin(1:2:end);
                values = varargin(2:2:end);
            end

            for j = 1:length(names)
                opt = validatestring(names{j}, gibbonSettings.names);
                value = values{j};

                % "Massage" values before assignment
                switch opt
                case 'ViewProfile'
                    % Allow case insensitive match, e.g. 'cad' for 'CAD'
                    value = validatestring(value, gibbonSettings.validViewProfiles);
                end
    
                % Validate based on property attributes
                gibbonSettings.validators.(opt)(value);
                
                % Update settings
                s = settings;
                if ~hasGroup(s,'GIBBON'), addGroup(s,'GIBBON'); end
                if ~hasSetting(s.GIBBON,opt), addSetting(s.GIBBON,opt); end
                s.GIBBON.(opt).PersonalValue = value;
            end
        end

        function reset()
        % Revert to defaults (remove settings.GIBBON group)
            s = settings;
            if hasGroup(s,'GIBBON'), s.removeGroup('GIBBON'); end
        end

        function p = findFEBioPath(defaultPaths)
        % Try to find FEBio executable on system path
        
            if nargin == 0
                defaultPaths = {'febio4','febio3'};
            else
                defaultPaths = cellstr(defaultPaths);
            end

            current = gibbonSettings.get('FEBioPath');
            if ~isempty(current)
                defaultPaths = [{current}, defaultPaths];
            end

            if ispc, checkCmd = 'where'; else, checkCmd = 'which'; end

            for j = 1:numel(defaultPaths)

                if isfile(defaultPaths{j})
                    d = dir(defaultPaths{j});
                    p = fullfile(d.folder, d.name);
                    return;
                end

                [status, stdout] = system([checkCmd ' ' defaultPaths{j}]);
                if status ~= 0
                    continue
                end
                lines = strsplit(strtrim(stdout),'\n');
                p = lines{end};
                if ~isfile(p)
                    continue
                end
                return;
            end

            warning('gibbon:Settings:FEBioPath','Failed to find system FEBio')
            p = '';
        end
    end
end

function names = getNames()
    mc = ?gibbonSettings;
    props = mc.PropertyList;
    bad = arrayfun(@(p) p.Constant || p.Dependent || p.Abstract, props);
    names = {props(~bad).Name};
end

function d = getDefaults()
    mc = ?gibbonSettings;
    opts = gibbonSettings.names;
    d = struct();
    for j = 1:length(opts)
        mp = findobj(mc.PropertyList,'Name',opts{j});
        d.(opts{j}) = mp.DefaultValue;
    end
end

function s = getValidators()
    mc = ?gibbonSettings;
    opts = gibbonSettings.names;
    s = struct();
    for j = 1:length(opts)
        mp = findobj(mc.PropertyList,'Name',opts{j});
        s.(opts{j}) = @(x) mp.Validation.validateValue(x);
    end
end

function mustBeFileOrEmpty(path)
    if isempty(path), return; end
    mustBeFile(path)
end

function warnIfFEBioVersion2(path)
    if contains(lower(path),'febio2')
        warning('FEBio2 detected. FEBio2 support is deprecated. Please upgrade to FEBio3');
    end
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
% Copyright (C) 2006-2026 Kevin Mattheus Moerman and the GIBBON contributors
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
