function tests = TEST_getSubPaths
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    rootDir = [tempname '.hasdot'];
    mkdir(rootDir);

    visible = {'foo', 'bar', fullfile('foo', 'nested')};
    hidden = {'.hidden', fullfile('foo', '.secrets'), '.git', fullfile('.git', 'stuff')};
    files = {'.gitignore', 'some.txt', fullfile('foo','bam.m')};

    visible = cellfun(@(d) fullfile(rootDir, d), visible, 'unif', 0);
    hidden = cellfun(@(d) fullfile(rootDir, d), hidden, 'unif', 0);
    cellfun(@mkdir, [visible, hidden]);

    for f = files
        fclose(fopen(fullfile(rootDir, f{1}), 'w'));
    end

    testCase.TestData.rootDir = rootDir;
    testCase.TestData.visible = visible;
    testCase.TestData.hidden = hidden;
end

function teardownOnce(testCase)
    if isfolder(testCase.TestData.rootDir)
	    rmdir(testCase.TestData.rootDir, 's');
    end
end

function testReturnsSubfoldersOnly(testCase)
    pathNames = getSubPaths(testCase.TestData.rootDir);

    verifyTrue(testCase, iscell(pathNames));
    expected = [testCase.TestData.visible, testCase.TestData.hidden];
    weird = setxor(pathNames, expected);
    verifyTrue(testCase, isempty(weird))
end

function testIgnoreHiddenFolders(testCase)
    pathNames = getSubPaths(testCase.TestData.rootDir, true);

    verifyTrue(testCase, iscell(pathNames));
    weird = setxor(pathNames, testCase.TestData.visible);
    verifyTrue(testCase, isempty(weird))
end
