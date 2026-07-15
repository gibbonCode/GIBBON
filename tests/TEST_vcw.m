function tests = TEST_vcw
    tests = functiontests(localfunctions);
end

function setup(testCase)
    testCase.TestData.hf = figure();
end

function teardown(testCase)
    if ishandle(testCase.TestData.hf)
        close(testCase.TestData.hf);
    end
end

function testViewProfiles(testCase)
% Check that VCW recognizes all gibbonSettings.validViewProfiles

    profiles = gibbonSettings.validViewProfiles;
    for j = 1:length(profiles)
        vcw(testCase.TestData.hf, profiles{j});
    end
end
