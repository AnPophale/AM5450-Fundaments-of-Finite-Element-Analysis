clc; clear all; close all; format compact; format shortg;

syms x y z;
eqn1 = 2*x - 3*y - z == 5;
eqn2 = 4*x + 4*y - 3*z == 3;
eqn3 = -2*x + 3*y - z ==1;

sol = solve([eqn1, eqn2, eqn3], [x,y,z]);
x = sol.x
y = sol.y
z = sol.z