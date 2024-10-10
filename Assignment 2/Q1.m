clc; clear all; close all; format compact; format shortg;

x0 = 0; xL = 3;
x_vec = x0:0.1:xL;

%Solutions using LSM and GM are the same for same order of trial function
%(When both essential BCs are satisfied by the trial function)

%First order trial function
u1 = 9*x_vec;
u1_prime = 9*ones(numel(x_vec),1);
r1 = -6*x_vec; %Residual

%Second order trial function
u2 = 4.5*(x_vec.^2 - x_vec);
u2_prime = 4.5*(2*x_vec - 1);
r2 = 9 - 6*x_vec;

%Third order trial function 
r3 = zeros(numel(x_vec),1);
u3 = x_vec.^3;
u3_prime = 3*x_vec.^2;

%Analytical Solution
u = x_vec.^3;
u_prime = 3*x_vec.^2;

%Plotting
%u vs x
figure;
hold on
plot(x_vec,u1,'-b')
plot(x_vec,u2,'-r')
plot(x_vec,u3,'-g')
plot(x_vec,u,'xk')
legend("Order 1","Order 2","Order 3","Exact")
xlabel("x")
ylabel("u")
title("u vs x")

%u' vs x
figure;
hold on
plot(x_vec,u1_prime,'-b')
plot(x_vec,u2_prime,'-r')
plot(x_vec,u3_prime,'-g')
plot(x_vec,u_prime,'xk')
legend("Order 1","Order 2","Order 3","Exact")
xlabel("x")
ylabel("u'")
title("u' vs x")

%e(x) vs x
figure;
hold on
plot(x_vec,r1,'-r')
plot(x_vec,r2,'-b')
plot(x_vec,r3,'-g')
legend("Order 1","Order 2","Order 3")
xlabel("x")
ylabel("Residual")
title("e(x) vs x")