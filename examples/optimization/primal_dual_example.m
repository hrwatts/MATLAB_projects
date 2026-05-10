% Example problems for primal_dual.
% Run from the repository root after calling addpath(genpath('src')).

A = [1 1 2; 2 1 3];
b = [3; 5];
c = [2; 1; 4];
lam = [0; 0];

[x_star_1, z_star_1, lam_star_1] = primal_dual(A, b, c, lam);
disp(x_star_1)
disp(z_star_1)
disp(lam_star_1)

A = [1 -2  1  3  2 -3  4  1  0  1  2  4  4 -1  0;
     0  2  4 -1  0  1  2 -1  1  2  0  2  1 -2  4;
     2  1  0  1  2 -1 -2  2  1  0 -1 -2  1  0  1;
    -2 -4  2  1 -2  1  0  1  2 -1  0  0  2  1 -1;
     1  0  2  4  1  0  1 -1  1  2  0  1 -1  2 -2];
b = [37; 30; 8; -4; 8];
c = [1; -10; 15; 20; 0; 0; 12; 3; 10; 8; 4; 12; 3; 10; -10];
lam = [1; 0; -1; 2; 4];

[x_star_2, z_star_2, lam_star_2] = primal_dual(A, b, c, lam);
disp(x_star_2)
disp(z_star_2)
disp(lam_star_2)
