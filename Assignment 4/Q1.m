clc; clear all; close all; format compact; format shortg;

tic
syms x
L = 2; x0 = 0; xL = L; 
xvec = x0:0.01:xL;

figure(1); hold on;
figure(2); hold on;
figure(3); hold on;

N = 1;
Rt_mean = 1;
while Rt_mean>10^-3
    N = N+1;
    [u,uprime,Rt] = galerkin(N);
    
    %Plotting only for polynomials above order 15
    if(N>15)
        figure(1) 
        plot(xvec,subs(u,x,xvec),'-', 'DisplayName', num2str(N), 'LineWidth', 1.3);
        
        figure(2)
        plot(xvec,subs(uprime,x,xvec),'-', 'DisplayName', num2str(N), 'LineWidth', 1.3);
    
        figure(3)
        plot(xvec,subs(Rt,x,xvec),'-', 'DisplayName', num2str(N), 'LineWidth', 1.3);
    end

    Rt_mean = mean(abs(subs(Rt,x,xvec)));
    Rt_vec(N-1) = Rt_mean;
end

figure(1); legend("show"); xlabel('x'); ylabel('u(x)'); title('u(x) vs x using Galerkin Method')
figure(2); legend("show"); xlabel('x'); ylabel('u''(x)'); title('u''(x) vs x using Galerkin Method')
figure(3); legend("show"); xlabel('x'); ylabel('Residual'); title('e(x) vs x using Galerkin Method')

%Plotting the residual vs order of trial function
figure(4)
semilogy(2:1:N,Rt_vec,'-o')
xlabel("Order of polynomial")
ylabel("Total Residue")
title("Total Residue vs Order of trial function ")
toc

function [u,uprime,Rt] = galerkin(N)
    syms x
    a = sym('a', [1,N+1]);
    b = zeros(1,N+1);

    %Constants
    L = 2; x0 = 0; xL = L; 
    E = 200*10^9; 
    A0 = 3.1416*10^-4; AL = 1.9635*10^-5;
    P = 2*10^3;
    c = 10^3;

    %Linear variation of area
    A = A0 + (x/L)*(AL-A0);
    
    %Trial Function
    u = poly2sym(flip(a));
    dudx = diff(u,x,1);
    
    %Essential BC
    BC1 = subs(u,x,x0) == 0;
    b(1) = solve(BC1,a(1));
    u = subs(u,a(1),b(1));
    
    %Domain Residual
    Rd = diff(E*A*dudx,x,1) + c*x;
    %Boundary Residual
    Rb = E*AL*subs(dudx,x,xL) - P;
    
    %Weighted Residual Equations
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
    
    u = subs(u,a,b); %Displacement
    uprime = diff(u,x,1); %Strain

    %Total Residual (Normalized)
    Rt = subs(Rd,a,b)/(c*L) + subs(Rb,a,b)/P;
end