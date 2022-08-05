%Example 1. min z = 2x1 + x2 + 3x3 - 2x4 + 10x5
% s.t. x1 +x3 –x4 +2x5=2, 
%      x2 +2x3 +2x4 +x5=3,
% d=(-1, -1, -1, -1, -1), 
% h = (6, 9, 0, 4, 2).

A=[1 0 1 -1 2;
   0 1 2 2 1]
b=[2;3]
c= [2 1 3 -2 10]'
d=[-1,-1,-1,-1,-1]'          
h=[6,9,0,4,2]'   



[x_star z_star]= bounded_simplex(A,b,c,d,h);

'Press key to continue.'
pause

% Example 2
%
%
%
%
%

A=[1 -2  1  3  2 -3  4  1  0  1  2  4  4 -1  0;
   0  2  4 -1  0  1  2 -1  1  2  0  2  1 -2  4;
   2  1  0  1  2 -1 -2  2  1  0 -1 -2  1  0  1;
  -2 -4  2  1 -2  1  0  1  2 -1  0  0  2  1 -1;
   1  0  2  4  1  0  1 -1  1  2  0  1 -1  2 -2]             
b=[37;  30;  8;  -4;  8]
c= [1  -10   15  20  0   0  12  3   10   8  4   12  3  10  -10]'
h=[4  4  1  5  6  0  4 10  4  6  4  6  4  8  8]'          
d=[0 -1 -2  0  1 -4 -1  4  0  0 -4  0 -2 -2  1]'   

[x_star z_star]= bounded_simplex(A,b,c,d,h);

'Press key to continue.'
pause

% % % Example 3
% % % 
% % % max z = 4x1 + 3x2
% % % s.t. x1 + x2 = 6
% % % 2x1 + x2 = 8
% % % 1 <= x1 <= 5 
% % % 1 <= x2 <= 4

A = [1 1; 2 1]
b = [6; 8]
c = [2;0]
d = [1; 1]
h = [5; 4]


[x_star z_star]= bounded_simplex(A,b,c,d,h);
