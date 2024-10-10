clc; clear all; close all; format compact; format shortg;

syms x y c1 c2 c3 c4 c5 c6;

%Constants
E = 80*10^9; %Young's modulus
nu = 0.3; %Poisson Ratio 
%E and nu for aluminium are taken from online data as they were not specified in the problem

q = 50; %Load
t = 10^-3; %Thickness
D = E*t^3/(12*(1-nu^2)); %Bending Stiffness

%Dimensions
x0 = 0; xL = 2; y0 = 0; yL = 2;

%Trial solution (Assumed in a way that EBCs are satisfied)
u = x*y*(x^2-4)*(y^2-4)*(c1 + c2*x + c3*y + c4*x^2 + c5*y^2 + c6*x*y);
eqns = sym(zeros(6, 1));

for i = 1:6
    %Weights for MGM
    weight = diff(u, eval(['c' num2str(i)]), 1);
    
    %Weighted residual equation for MGM
    term1 = int(int(diff(u,x,2) * diff(weight,x,2), x,x0,xL), y,y0,yL);
    term2 = 2*int(int(diff(diff(u,x,1),y,1) * diff(diff(weight,x,1),y,1), x,x0,xL), y,y0,yL);
    term3 = int(int(diff(u, y, 2) * diff(weight, y, 2), y, 0, 2), x, 0, 2);
    term4 = int(int(weight*(q/D), x,x0,xL), y,y0,yL);
    eqns(i) = term1 + term2 + term3 + term4;
end

%Solving the system of equations
sol = solve(eqns, [c1 c2 c3 c4 c5 c6]);
u = subs(u, [c1 c2 c3 c4 c5 c6], [sol.c1 sol.c2 sol.c3 sol.c4 sol.c5 sol.c6]);
disp(vpa(u,4)); %Displaying the solution

%Plotting the solution
figure;
fsurf(u, [0 2 0 2]);
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('u(x,y) vs x,y in 3D');
colormap(jet); colorbar;
grid on;

figure;
fcontour(u, [0 2 0 2], 'LineWidth', 2);
xlabel('x');
ylabel('y');
title('Contour plot for u(x,y) vs x,y');
colormap(jet); colorbar;
grid on;
