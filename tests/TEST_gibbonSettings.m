function tests = TEST_gibbonSettings
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)

    % Record settings.GIBBON before test
    s = settings;
    if s.hasGroup('GIBBON')
        testCase.TestData.originalSettings = backupSettings(s.GIBBON);
    else
        testCase.TestData.originalSettings = [];
    end

    % Record warnings state before test
    testCase.TestData.warnings = warning;
end

function setup(~)
    s = settings;
    if s.hasGroup('GIBBON')
        s.removeGroup('GIBBON')
    end
end

function teardownOnce(testCase)

    % Restore settings.GIBBON
    s = settings;
    if isempty(testCase.TestData.originalSettings)
        if s.hasGroup('GIBBON')
	        s.removeGroup('GIBBON')
        end
    else
        if ~s.hasGroup('GIBBON')
            s.addGroup('GIBBON')
        end
        restoreSettings(testCase.TestData.originalSettings, s.GIBBON)
    end

    % Restore warnings
    warning(testCase.TestData.warnings);
end

function s = backupSettings(obj)
    s = struct();
    p = properties(obj);
    for j = 1:numel(p)
        try
            v = obj.(p{j});
        catch
            continue
        end
        if isobject(v)
            s.(p{j}) = backupSettings(v);
        else
            s.(p{j}) = v;
        end
    end
end

function restoreSettings(backup, obj)
    f = fieldnames(backup);
    for j = 1:numel(f)
        v = obj.(f{j});
        b = backup.(f{j});
        if isobject(v)
            restoreSettings(b, v);
        else
            try %#ok<TRYNC>
                obj.(f{j}) = b;
            end
        end
    end
end

function testConstructor(testCase)

    warning('off','gibbon:Settings:FEBioPath')
    g = gibbonSettings();

    f = gibbonSettings.names;
    for j = 1:length(f)
        verifyEqual(testCase, g.(f{j}), gibbonSettings.defaults.(f{j}))
    end
end

function testGetSet(testCase)

    warning('off','gibbon:Settings:FEBioPath')
    g = gibbonSettings();

    set(g,'view','touch');
    set(g,'febio','febio3');

    verifyEqual(testCase, g.ViewProfile, 'touchpad')
    verifyEqual(testCase, g.FEBioPath, 'febio3')

    h = gibbonSettings();
    verifyEqual(testCase, h.ViewProfile, g.ViewProfile)
    verifyEqual(testCase, h.FEBioPath, g.FEBioPath)
end

function testReset(testCase)

    warning('off','gibbon:Settings:FEBioPath')
    g = gibbonSettings();

    set(g,'ViewProfile','touchpad');
    set(g,'FEBioPath','febio3');

    g.reset();

    f = gibbonSettings.names;
    for j = 1:length(f)
        verifyEqual(testCase, g.(f{j}), gibbonSettings.defaults.(f{j}))
    end
end

function testInvalidOption(testCase)

    g = gibbonSettings();
    verifyError(testCase, @() get(g, 'notAnOption'), 'gibbon:Settings:InvalidOption')
    verifyError(testCase, @() set(g, 'notAnOption', 'value'), 'gibbon:Settings:InvalidOption')
end

function testInvalidViewProfileValue(testCase)

    g = gibbonSettings();
    verifyError(testCase, @() set(g, 'ViewProfile', 'not-a-profile'), 'gibbon:Settings:InvalidViewProfile')
end

function testInvalidFEBioPathValue(testCase)

    warning('off','gibbon:Settings:FEBioPath')
    g = gibbonSettings();

    warning('error','gibbon:Settings:FEBioPath')
    verifyError(testCase, @() set(g, 'FEBioPath', 'definitely_not_a_real_febio_binary'), 'gibbon:Settings:FEBioPath')
end