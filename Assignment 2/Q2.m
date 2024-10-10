clc; clear all; close all; format compact; format shortg;

%Parameters
syms x
k = 5;
ubar = 10;
x0 = 0; xL = 1;

%True Solution
x_vec = x0:0.005:xL;
u_true = 2*sin(5*(x0:0.05:xL))/cos(5);

%Galerkin Method
figure(1); hold on;
figure(2); hold on;

Rt1_mean = 1; %Mean residual 
N1 = 1;
while (Rt1_mean > 1e-6)
    N1 = N1+1;
    [u1,Rt1] = galerkin(N1, x0, xL, k, ubar);

    figure(1);
    plot(x_vec, subs(u1, x_vec), '-', 'DisplayName', num2str(N1), 'LineWidth', 1.3);

    figure(2);
    plot(x_vec, subs(Rt1, x_vec), '-', 'DisplayName', num2str(N1), 'LineWidth', 1.3);
    Rt1_mean = mean(abs(subs(Rt1,x_vec)));
    fprintf('Error for Galerkin Method = %.5e with an order %d polynomial\n', Rt1_mean,N1);
end
fprintf('\n');

figure(1);
plot(x0:0.05:xL, u_true, 'xk', 'DisplayName', 'Exact');

figure(1); legend("show"); xlabel('x'); ylabel('u(x)'); title('u(x) vs x using Galerkin Method')
figure(2); legend("show"); xlabel('x'); ylabel('Residual'); title('e(x) vs x using Galerkin Method')

%Least Squares Method
figure(3); hold on;
figure(4); hold on;

Rt2_mean = 1; %Mean Residual
N2 = 1;
while (Rt2_mean > 1e-6)
    N2 = N2+1;
    [u2,Rt2] = least_squares(N2, x0, xL, k, ubar);

    figure(3);
    plot(x_vec, subs(u2, x_vec), '-', 'DisplayName', num2str(N2), 'LineWidth', 1.3);

    figure(4);
    plot(x_vec, subs(Rt2, x_vec), '-', 'DisplayName', num2str(N2), 'LineWidth', 1.3);
    Rt2_mean = mean(abs(subs(Rt2,x_vec)));
    fprintf('Error for Least Squares Method = %.5e with an order %d polynomial\n', Rt2_mean,N2);
end

figure(3);
plot(x0:0.05:xL, u_true, 'xk', 'DisplayName', 'Exact');

figure(3); legend("show"), xlabel('x'); ylabel('u(x)'); title('u(x) vs x using Least Squares Method')
figure(4); legend("show"); xlabel('x'); ylabel('Residual'); title('e(x) vs x using Least Squares Method')

%Function Definitions
%Galerkin Method
function [u,Rt] = galerkin(N,x0,xL,k,ubar)   
    syms x 
    a = sym('a', [1,N+1]);
    b = zeros(1,N+1);

    %Trial Function 
    u = poly2sym(flip(a));
    dudx = diff(u,x,1);
    du2dx2 = diff(u,x,2);

    %Essential Boundary Condition
    BC1 = subs(u,x0) == 0;
    b(1) = solve(BC1,a(1));

    %Domain Residual
    Rd = du2dx2 + k^2*subs(u,a(1),b(1));
    %Boundary Residual
    Rb = subs(dudx - ubar,x,xL); %Natural BC
    
    %Weighted residual equations
    eq = sym(zeros(N,1));
    for i = 2:N+1
        w = diff(u,a(i),1);
        wb = subs(diff(u,a(i),1),x,xL);
        eq(i-1) = int(w*Rd,x,x0,xL) + wb*Rb == 0;
    end
    
    sol = solve(eq,a(2:end));
    for i = 2:N+1
        b(i) = sol.(['a',num2str(i)]);
    end
    
    u = subs(u,a,b);
    %Total Residual
    Rt = subs(Rd,a,b) + subs(Rb,a,b); 
end

%Least Squares Method
function [u,Rt] = least_squares(N,x0,xL,k,ubar)   
    syms x 
    a = sym('a', [1,N+1]);
    b = zeros(1,N+1);

    %Trial Function 
    u = poly2sym(flip(a));
    dudx = diff(u,x,1);
    du2dx2 = diff(u,x,2);

    %Essential BC
    BC1 = subs(u,x0) == 0;
    b(1) = solve(BC1,a(1));

    %Domain Residual
    Rd = du2dx2 + k^2*subs(u,a(1),b(1));    
    %Boundary Residual
    Rb = subs(dudx - ubar,x,xL); %Natural BC

    %Weighted residual equations
    eq = sym(zeros(N,1));
    for i = 2:N+1
        w = diff(Rd,a(i),1);
        wb = diff(Rb,a(i),1);
        eq(i-1) = int(w*Rd,x,x0,xL) + wb*Rb == 0;
    end

    sol = solve(eq,a(2:end));
    for i = 2:N+1
        b(i) = sol.(['a',num2str(i)]);
    end

    u = subs(u,a,b); 
    %Total Residual
    Rt = subs(Rd,a,b) + subs(Rb,a,b);
end