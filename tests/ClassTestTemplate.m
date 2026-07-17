classdef ClassTestTemplate < matlab.unittest.TestCase
% Example of class-based unit tests. See:
% https://www.mathworks.com/help/matlab/matlab_prog/class-based-unit-tests.html
% https://www.mathworks.com/help/matlab/matlab_prog/ways-to-write-unit-tests.html
%
% If you're new to tests, you may want to start with function-based unit tests
% See also: functionTestTemplate

    %% (Optional) Test Parameters
    % TestParameter properties are used to define parameterized tests,
    % i.e. "variations" of a single test method, ran with different inputs.
    % See the 'testParameterizedMultiply' method below, and:
    % https://www.mathworks.com/help/matlab/matlab_prog/use-parameters-in-class-based-tests.html

    properties (TestParameter)
        b = {3,4}
        a = struct('one', 1, 'two', 2)  % use a struct to name parameter values 
    end

    %% (Optional) Custom Properties
    % Define your own properties to pass data between fixtures and test methods.

    properties
        TempDir  % used to store a temporary directory for the tests
        oldPath  % used to store the current directory before each test
    end

    %% (Optional) Fixtures
    % Fixtures run before and after tests, to do preparation or cleanup.

    methods (TestClassSetup)
        function setupOnce(testCase)
            % Runs once for the whole class.
            % Use for expensive one-time setup (temp dirs, large test data, etc.).
            tempDir = tempname;
            mkdir(tempDir);
            testCase.TempDir = tempDir;
        end
    end

    methods (TestClassTeardown)
        function teardownOnce(testCase)
            % Runs once after all tests in this class.
            if isfolder(testCase.TempDir)
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods (TestMethodSetup)
        function setup(testCase)
            % Runs before each test method.
            testCase.oldPath = pwd;
            cd(testCase.TempDir);
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            % Runs after each test method.
            cd(testCase.oldPath);
        end
    end

    %% Test Methods

    methods (Test)
        function testWritingFilesWorks(testCase)
  
            writelines('Class tests are great!', 'test.txt');
            testCase.verifyTrue(isfile('test.txt'));
        end

        function testParameterizedMultiply(testCase, a, b)
        % a, b will take the values defined in the TestParameter properties above

            c = multiply(a,b);
            testCase.assertEqual(c, a*b);
        end
    end
end

function c = multiply(a, b)
    c = a*b;
end