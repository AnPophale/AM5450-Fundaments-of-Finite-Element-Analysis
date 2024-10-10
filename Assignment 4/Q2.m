clc; clear all; close all; format compact; format shortg;

tic
% Constants
syms x;
x0 = 0; xL = 2; L = 2;
A0 = 3.1416e-4; AL = 1.9635e-5;
A = A0 + (x/L)*(AL-A0);
E = 200e9;
P = 2e3;
c = 1e3;

%Starting with 2 elements
Ne = 2; 
[u_old, uprime_old, xvec_old] = element(Ne);

Re1 = 1;
Re2 = 1;
Re1_vec = [];
Re2_vec = [];
Re3_vec = [];
Re4_vec = [];
while (Re1 > 1e-3 && Re2 > 1e-3)
    Ne = Ne + 1;
    [u_new, uprime_new, xvec_new] = element(Ne); 
    
    % Relative error in u at bar tip
    Re1 = abs((u_new(end) - u_old(end)) / u_new(end));
    Re1_vec(Ne-2) = Re1;
    
    % u at bar midpoint
    if mod(Ne, 2) == 0
        u_mid_old = u_old(Ne/2 + 1);
        u_mid_new = u_new(Ne/2 + 1);
    else %Linear interpolaton 
        u_mid_old = 0.5*(u_old((Ne+1)/2) + u_old((Ne+3)/2));
        u_mid_new = 0.5*(u_new((Ne+1)/2) + u_new((Ne+3)/2));
    end
    
    % Relative error in u at bar midpoint
    Re2 = abs((u_mid_new - u_mid_old) / u_mid_new);
    Re2_vec(Ne-2) = Re2;

    % Relative error in u' at bar tip
    Re3 = abs((uprime_new(end) - uprime_old(end)) / uprime_new(end));
    Re3_vec(Ne-2) = Re3;
    
    % u' at bar midpoint
    if mod(Ne, 2) == 0 %Linear interpolation
        uprime_mid_old = 0.5*(uprime_old(Ne/2) + uprime_old(Ne/2 + 1));
        uprime_mid_new = 0.5*(uprime_new(Ne/2) + uprime_new(Ne/2 + 1));
    else
        uprime_mid_old = uprime_old((Ne+1)/2);
        uprime_mid_new = uprime_new((Ne+1)/2);
    end

    %Relative error in u' at bar midpoint
    Re4 = abs((uprime_mid_new - uprime_mid_old) / uprime_mid_new);
    Re4_vec(Ne-2) = Re4;

    u_old = u_new;
    uprime_old = uprime_new;
    xvec_old = xvec_new;
end

fprintf('Relative Error in u at tip: %.5e with %d elements\n', Re1, Ne);
fprintf('Relative Error in u at midpoint: %.5e with %d elements\n', Re2, Ne);
fprintf('Relative Error in u'' at tip: %.5e with %d elements\n', Re3, Ne);
fprintf('Relative Error in u'' at midpoint: %.5e with %d elements\n', Re4, Ne);

figure(1);
hold on
plot(3:1:Ne,Re1_vec,'-o','DisplayName','u tip','LineWidth', 1.3)
plot(3:1:Ne,Re2_vec,'-o','DisplayName','u mid','LineWidth', 1.3)
plot(3:1:Ne,Re3_vec,'-o','DisplayName','u'' tip','LineWidth', 1.3)
plot(3:1:Ne,Re4_vec,'-o','DisplayName','u'' mid','LineWidth', 1.3)
set(gca,"YScale","log")
xlim([3,Ne])
legend('show')
xlabel('Number of elements')
ylabel('Relative Error')
title('Relative Error in displacement and strain at bar tip and midpoint')

figure(2); hold on;
figure(3); hold on;
for N = 2:3:Ne
     [u, uprime, xvec] = element(N);

     figure(2)
     plot(xvec,u,'-', 'DisplayName', num2str(N), 'LineWidth', 1.3)

     figure(3)
     stairs(xvec,[uprime;uprime(end)],'-', 'DisplayName', num2str(N), 'LineWidth', 1.3)
end
figure(2); legend("show"); xlabel('x'); ylabel('u(x)'); title('u(x) vs x using different number of elements')
figure(3); legend("show"); xlabel('x'); ylabel('u''(x)'); title('u''(x) vs x with different number of elements')
toc

function [u, uprime, xvec] = element(Ne)
    syms x;
    x0 = 0; xL = 2; L = 2;
    A0 = 3.1416e-4; AL = 1.9635e-5;
    A = A0 + (x/L)*(AL-A0);
    E = 200e9;
    P = 2e3;
    c = 1e3;

    xvec = x0:L/Ne:xL; % Node locations

    % Initialize matrices
    K = zeros(Ne+1, Ne+1);
    rq = zeros(Ne+1, 1);
    rb = zeros(Ne+1, 1);
    rb(end) = P;

    for i = 1:Ne
        A1 = subs(A, x, xvec(i));
        A2 = subs(A, x, xvec(i+1));
        Am = (A1 + A2)/2;
        l = xvec(i+1) - xvec(i);
    
        K(i, i) = K(i, i) + E*Am/l;
        K(i, i+1) = K(i, i+1) - E*Am/l;
        K(i+1, i) = K(i+1, i) - E*Am/l;
        K(i+1, i+1) = K(i+1, i+1) + E*Am/l;
    
        rq(i) = rq(i) + 0.25*l*c*(xvec(i) + xvec(i+1) - l/3);
        rq(i+1) = rq(i+1) + 0.25*l*c*(xvec(i) + xvec(i+1) - l/3);
    end
   
    % Solve for displacement
    u = zeros(Ne+1, 1);
    u(2:end) = K(2:end, 2:end)\(rq(2:end) + rb(2:end));

    % Calculate strain
    uprime = zeros(Ne, 1);
    for i = 1:Ne
        l = xvec(i+1) - xvec(i);
        uprime(i) = (u(i+1) - u(i)) / l;
    end
end