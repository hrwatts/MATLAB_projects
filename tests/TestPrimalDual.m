classdef TestPrimalDual < matlab.unittest.TestCase
    methods (Test)
        function returnsExpectedShapes(testCase)
            addpath(genpath('src'));
            A = [1 1 2; 2 1 3];
            b = [3; 5];
            c = [2; 1; 4];
            lam = [0; 0];
            [x_star, z_star, lam_star] = primal_dual(A, b, c, lam);
            testCase.verifySize(x_star, [3 1]);
            testCase.verifySize(lam_star, [2 1]);
            testCase.verifyTrue(isfinite(z_star));
        end
    end
end
