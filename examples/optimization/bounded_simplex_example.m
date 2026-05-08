% Example problems for bounded_simplex.
% Run from the repository root after calling addpath(genpath('src')).

A = [1 0 1 -1 2;
     0 1 2 2 1];
b = [2; 3];
c = [2; 1; 3; -2; 10];
d = [-1; -1; -1; -1; -1];
h = [6; 9; 0; 4; 2];

[x_star_1, z_star_1] = bounded_simplex(A, b, c, d, h);
disp(x_star_1)
disp(z_star_1)

A = [1 -2  1  3  2 -3  4  1  0  1  2  4  4 -1  0;
     0  2  4 -1  0  1  2 -1  1  2  0  2  1 -2  4;
     2  1  0  1  2 -1 -2  2  1  0 -1 -2  1  0  1;
    -2 -4  2  1 -2  1  0  1  2 -1  0  0  2  1 -1;
     1  0  2  4  1  0  1 -1  1  2  0  1 -1  2 -2];
b = [37; 30; 8; -4; 8];
c = [1; -10; 15; 20; 0; 0; 12; 3; 10; 8; 4; 12; 3; 10; -10];
h = [4; 4; 1; 5; 6; 0; 4; 10; 4; 6; 4; 6; 4; 8; 8];
d = [0; -1; -2; 0; 1; -4; -1; 4; 0; 0; -4; 0; -2; -2; 1];

[x_star_2, z_star_2] = bounded_simplex(A, b, c, d, h);
disp(x_star_2)
disp(z_star_2)

A = [1 1; 2 1];
b = [6; 8];
c = [2; 0];
d = [1; 1];
h = [5; 4];

[x_star_3, z_star_3] = bounded_simplex(A, b, c, d, h);
disp(x_star_3)
disp(z_star_3)
