function tests = TEST_viewText
    tests = functiontests(localfunctions);
end

function setup(testCase)
    testCase.TestData.handle = [];
end

function teardown(testCase)
    if ishandle(testCase.TestData.handle)
        delete(testCase.TestData.handle)
    end
end

function testViewTextFromString(testCase)

    h = textView({'Hello'});
    testCase.TestData.handle = h;
    verifyInstanceOf(testCase,h,'matlab.ui.Figure')
end

function testViewTextFromFile(testCase)

    file = which('textView');
    h = textView(file);
    testCase.TestData.handle = h;
    verifyInstanceOf(testCase,h,'matlab.ui.Figure')
end

function testViewTakesStructure(testCase)

    s = struct('BackgroundColor',[1,0,0]);
    h = textView({'Hello'}, s);
    testCase.TestData.handle = h;
    verifyInstanceOf(testCase,h,'matlab.ui.Figure')
    verifyEqual(testCase, h.Children.Children.BackgroundColor, [1,0,0])
end

function testViewTakesNameValue(testCase)

    h = textView({'Hello'}, FontColor=[1,0,0]);
    testCase.TestData.handle = h;
    verifyInstanceOf(testCase,h,'matlab.ui.Figure')
    verifyEqual(testCase, h.Children.Children.FontColor, [1,0,0])
end
