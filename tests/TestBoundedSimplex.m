classdef TestBoundedSimplex < matlab.unittest.TestCase
    methods (Test)
        function returnsColumnVectorAndScalar(testCase)
            addpath(genpath('src'));
            A = [1 1; 2 1];
            b = [6; 8];
            c = [2; 0];
            d = [1; 1];
            h = [5; 4];
            [x_star, z_star] = bounded_simplex(A, b, c, d, h);
            testCase.verifySize(x_star, [2 1]);
            testCase.verifyTrue(isfinite(z_star));
        end
    end
end
