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

    % Create a dummy FEBio executable 
    testCase.TestData.dummyFEBio = tempname;
    fclose(fopen(testCase.TestData.dummyFEBio,'w'));
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

    if isfile(testCase.TestData.dummyFEBio)
        delete(testCase.TestData.dummyFEBio);
    end
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

function testDefaults(testCase)

    g = gibbonSettings.get();

    f = gibbonSettings.names;
    for j = 1:length(f)
        verifyEqual(testCase, g.(f{j}), gibbonSettings.defaults.(f{j}))
    end
end

function testGetSet(testCase)

    gibbonSettings.set('view','touch');
    gibbonSettings.set('febio', testCase.TestData.dummyFEBio);

    g = gibbonSettings.get();
    verifyEqual(testCase, g.ViewProfile, 'touchpad');
    verifyEqual(testCase, g.FEBioPath, testCase.TestData.dummyFEBio);

    h = gibbonSettings.get();
    verifyEqual(testCase, g, h)
end

function testFindFEBioPicksFirstValid(testCase)

    options = {'definitely_not_a_file', testCase.TestData.dummyFEBio, 'foo'};
    p = gibbonSettings.findFEBioPath(options);
    verifyEqual(testCase, p, testCase.TestData.dummyFEBio);
end

function testFindFEBioKeepsExisting(testCase)

    gibbonSettings.set('FEBioPath', testCase.TestData.dummyFEBio);
    p = gibbonSettings.findFEBioPath();
    verifyEqual(testCase, p, testCase.TestData.dummyFEBio);
end

function testReset(testCase)

    gibbonSettings.set('view','touch','febio', testCase.TestData.dummyFEBio);
    gibbonSettings.set('view','touch');

    gibbonSettings.reset();
    g = gibbonSettings.get();

    f = gibbonSettings.names;
    for j = 1:length(f)
        verifyEqual(testCase, g.(f{j}), gibbonSettings.defaults.(f{j}))
    end
end

function testInvalidOption(testCase)

    verifyError(testCase, @() gibbonSettings.get('notAnOption'), 'MATLAB:unrecognizedStringChoice')
    verifyError(testCase, @() gibbonSettings.set('notAnOption', 'value'), 'MATLAB:unrecognizedStringChoice')
end

function testInvalidViewProfileValue(testCase)

    verifyError(testCase, @() gibbonSettings.set('ViewProfile', 'not-a-profile'), 'MATLAB:unrecognizedStringChoice')
end

function testInvalidFEBioPathValue(testCase)

    verifyError(testCase, @() gibbonSettings.set('FEBioPath', 'definitely_not_a_file'), 'MATLAB:validators:mustBeFile')
end