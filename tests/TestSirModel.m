classdef TestSirModel < matlab.unittest.TestCase
    methods (Test)
        function conservesPopulationDerivative(testCase)
            addpath(genpath('src'));
            dxdt = sir_model(0, [0.9; 0.1; 0], 0.5, 0.2);
            testCase.verifyEqual(sum(dxdt), 0, 'AbsTol', 1e-12);
            testCase.verifySize(dxdt, [3 1]);
        end
    end
end
