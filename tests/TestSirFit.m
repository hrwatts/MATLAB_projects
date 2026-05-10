classdef TestSirFit < matlab.unittest.TestCase
    methods (Test)
        function returnsSimulationOutputs(testCase)
            addpath(genpath('src'));
            testCase.assumeFalse(isempty(ver('symbolic')), 'Symbolic Math Toolbox is required for sir_fit.');
            [tspan, curve, k, a] = sir_fit('sample-data/sir_covid/al_covid.csv', 227900);
            testCase.verifyGreaterThan(numel(tspan), 0);
            testCase.verifyEqual(size(curve, 2), 3);
            testCase.verifyGreaterThan(k, 0);
            testCase.verifyGreaterThan(a, 0);
        end
    end
end
