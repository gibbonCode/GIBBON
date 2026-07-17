function tests = functionTestTemplate
% Example of function-based unit tests. See:
% https://www.mathworks.com/help/matlab/matlab_prog/function-based-unit-tests.html
% https://www.mathworks.com/help/matlab/matlab_prog/ways-to-write-unit-tests.html
%
% See also: ClassTestTemplate

    % In most cases this line shouldn't change
    % It tells MATLAB to look for tests in the local functions of this file
    tests = functiontests(localfunctions);
end

%% (Optional) fixtures
% Fixtures are functions that run before and after tests, to do preparation or cleanup.
% The names 'setupOnce', 'teardownOnce', 'setup', and 'teardown' are reserved,
% (that's how MATLAB recognizes them as fixtures).
% The testCase input argument will contain a TestCase object,
% most useful to pass data between fixtures and test functions.

function setupOnce(testCase)
% Runs once, at the beginning of the test session.
% A typical use case is e.g. creating a temporary directory

    tempDir = tempname;
    mkdir(tempDir);

    % Fields stored in testCase.TestData will be available
    % to test functions and to other fixtures
    testCase.TestData.myTempDir = tempDir;
end

function teardownOnce(testCase)
% Runs once, at the end of the test session.
% It should 'undo' anything done in setupOnce, e.g. delete our temporary directory

    if isfolder(testCase.TestData.myTempDir)
        rmdir(testCase.TestData.myTempDir, 's');
    end
end

function setup(testCase)
% Runs before each test function is executed.
% Useful to 'reset the stage' before each test, e.g. start from a known directory:

    testCase.TestData.oldPath = pwd;
    cd(testCase.TestData.myTempDir);
end

function teardown(testCase)
% Runs after each test function is executed. Should undo anything 'setup' (or the test itself)
% might have done and that we don't want interfering with other tests, e.g.:

    cd(testCase.TestData.oldPath);
end


%% Test functions
% Each local function starting with 'test' will be treated as a test function.
% The testCase input argument will contain a TestCase object, which has methods for verifying
% conditions and reporting results.

function testWritingAFileWorks(testCase)

    % Calls to the code you want to test
    cd(testCase.TestData.myTempDir);
    writelines('Unit-tests are great!', 'test.txt');

    % plain assertions can be used, but the TestCase.verify* methods are more flexible
    % and will report failures in a way that is compatible with MATLAB's test runner
    testCase.verifyTrue(isfile('test.txt'));
end

function testHelperFunction(testCase)

    % Another example, this time using a helper function defined below
    result = some_helper_function(6, 7);
    testCase.verifyEqual(result, 42, 'AbsTol', 1e-6);
end

%% (Optional) helper functions
% Local functions that don't start with 'test' or use the reserved names work normally.

function c = some_helper_function(a, b)
    c = a * b + 1e-10 * rand();
end