classdef TestNewtonMethod < matlab.unittest.TestCase
    methods (Test)
        function convergesOnSimpleRoot(testCase)
            addpath(genpath('src'));
            [solution, history] = newton_m(@(x) x.^2 - 2, @(x) 2*x, 1, 10);
            testCase.verifyGreaterThanOrEqual(numel(history), 1);
            testCase.verifyEqual(solution, sqrt(2), 'AbsTol', 1e-8);
        end
    end
end
