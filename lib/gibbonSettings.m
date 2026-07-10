classdef gibbonSettings < matlab.mixin.SetGet
% Static-like* class to hold Gibbon settings,
% meant to be a single-point interface to set/get global GIBBON options.
%
% It currently uses MATLAB's `settings` for persistent storage,
% but the "back-end" could eventually be moved e.g. to a JSON file, for
% compatibility with Octave.
%
% (*) The only reason the class has a contstructor and (instance) get/set 
% methods is to take advantage of matlab.mixin.SetGet's auto-validation
% based on property attributes. For all practical purposes, however, the
% class should be considered to be a static interface to a global set of
% options. These two blocks are equivalent:
%
%   g = gibbonSettings;
%   set(g,'FEBiopath','febio3')
%   g.FEBioPath  % also get(g, 'FEBIOpath')
%
%   gibbonSettings('FEBiopath','febio3')
%   get(gibbonSettings,'FEBIOpath')
%
% See also: settings

    properties (Constant)
        gibbonPath = fileparts(fileparts(mfilename('fullpath')));
    end

    properties (SetAccess=protected)
        FEBioPath (1,:) char {checkFEBioPath} = findFEBio({'febio4','febio3',''})
        ViewProfile (1,:) char {validViewProfile} = 'CAD'
        % FutureOption ... = defaultValue
    end

    properties (Constant, Hidden)
        % constants for validation
        validViewProfiles = {'default', 'CAD', 'febio', 'touchpad'}

        % auto-populated from the main properties above
        defaults (1,1) struct = getDefaults()
        names (1,:) cell = getNames()
    end

    methods

        function obj = gibbonSettings(args)
        % (Re)create settings.GIBBON structure
        %
        % gibbonSettings() - will take values from settings.GIBBON,
        %       and complete with defaults.
        % gibbonSettings(Name, Value, ...) - override one or more settings.
        %       Case-insensitive [short-hand] names are allowed.
        % 
        % Examples:
        %   gibbonSettings('febio','febio4');   % overrides FEBioPath
        %   gibbonSettings(ViewProfile='cad');  % sets ViewProfile='CAD'
        %
        % See also: gibbonSettings.reset, gibbonSettings.restore

            arguments
                args.FEBioPath
                args.ViewProfile
            end

            % Modify existing handle, if available
            persistent lastObj;
            if isempty(lastObj) || ~isvalid(lastObj)
                lastObj = obj;
            end
            obj = lastObj;

            opts = gibbonSettings.names;
            for j = 1:length(opts)
                val = get(obj, opts{j});
                if isfield(args, opts{j})
                    val = args.(opts{j});
                end
                set(obj, opts{j}, val);
            end
        end

        function value = get(~, opt)
        % Get option from settings.GIBBON, fall back to default
            opt = validOptionName(opt);
            value = gibbonSettings.defaults.(opt);
            
            s = settings;
            if hasGroup(s,'GIBBON') && hasSetting(s.GIBBON,opt)
                v = s.GIBBON.(opt).PersonalValue;
                if ~isempty(v)
                    value = v; 
                end
            end
        end

        function set(obj, opt, value)
        % Write option to settings.GIBBON
            opt = validOptionName(opt);

            % Massage value before assignment
            switch opt
            case 'ViewProfile'
                % Allows case insensitive match, e.g. 'cad' for 'CAD'
                value = validViewProfile(value);
            end

            % Validate based on property attributes
            set@matlab.mixin.SetGet(obj, opt, value);
            
            % Update settings
            s = settings;
            if ~hasGroup(s,'GIBBON'), addGroup(s,'GIBBON'); end
            if ~hasSetting(s.GIBBON,opt), addSetting(s.GIBBON,opt); end
            s.GIBBON.(opt).PersonalValue = value;
        end
    end

    methods(Static)

        function obj = reset()
        % Revert to defaults (remove settings.GIBBON group)
            s = settings;
            if hasGroup(s,'GIBBON'), s.removeGroup('GIBBON'); end
            obj = gibbonSettings();
        end

        function s = copy()
        % Convert settings to structure, e.g. for use with
        % gibbonSettings.restore later.

            h = gibbonSettings();
            vals = cellfun(@(n) get(h,n), h.names, 'unif', 0);
            s = cell2struct(vals, h.names);
        end

        function restore(s)
        % Apply settings from structure
            args = [fieldnames(s), struct2cell(s)]';
            gibbonSettings(args{:});
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
        d.(opts{j}) = mp.Validation.validateValue(mp.DefaultValue);
    end
end

function varargout = validViewProfile(value)
    try
        value = validatestring(value, gibbonSettings.validViewProfiles);
    catch
        error('gibbon:Settings:InvalidViewProfile','Invalid ViewProfile value "%s".', value);
    end
    if nargout > 0, varargout{1} = value; end
end

function opt = validOptionName(opt)
    try
        opt = validatestring(opt, gibbonSettings.names);
    catch
        error('gibbon:Settings:InvalidOption','Invalid option name "%s".', opt);
    end
end

function mustBeSystemPath(path)

    if ispc, checkCmd = 'where'; else, checkCmd = 'which'; end
    [status, ~] = system([checkCmd ' ' path]);
    if status ~= 0
        error('gibbon:Settings:SystemPath','path "%s" is not recognized by the system', path)
    end
end

function checkFEBioPath(path)
% Warns if path or app name is not known by the system
% An empty PATH will NOT throw warnings

    if isempty(path), return; end

    try
        mustBeSystemPath(path)
    catch err
        warning('gibbon:Settings:FEBioPath', 'FEBio %s', err.message)
    end
end

function p = findFEBio(paths)
% Called during class load, sets the first working path as default

    p = '';
    paths = cellstr(paths);
    for j = 1:numel(paths)
        try
            mustBeSystemPath(paths{j})
            p = paths{j};
            break;
        catch
            continue
        end
    end
end
