% % % % Ex. Use the Primal-Dual algorithm to get x* [and lam*].
% % % %  minimize 2x1 + x2 + 4x3
% % % %  subject to x1 + x2 + 2x3 = 3
% % % %  2x1 + x2 + 3x3 = 5
% % % %  x1 >= 0, x2 >= 0, x3 >= 0

A = [1 1 2; 2 1 3]

b = [3; 5]

c = [2;1;4]

lam = [0;0]


[x_star z_star lam_star] = primal_dual(A,b,c,lam);


'Press key to continue.'
pause

% Example 2
% 
% 
% 

A=[1 -2  1  3  2 -3  4  1  0  1  2  4  4 -1  0;
   0  2  4 -1  0  1  2 -1  1  2  0  2  1 -2  4;
   2  1  0  1  2 -1 -2  2  1  0 -1 -2  1  0  1;
  -2 -4  2  1 -2  1  0  1  2 -1  0  0  2  1 -1;
   1  0  2  4  1  0  1 -1  1  2  0  1 -1  2 -2]
b=[37  30  8  -4  8]'
c= [1  -10   15  20  0   0  12  3   10   8  4   12  3  10  -10]'              
lam=[1 0 -1 2 4]'

[x_star z_star lam_star] = primal_dual(A,b,c,lam);
