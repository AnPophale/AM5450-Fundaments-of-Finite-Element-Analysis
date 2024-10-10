clc; clear all; close all; format compact; format shortg;

%Constants
L = 12;   
E = 200*10^9;
A = 300*500*10^-6;
I =  (500*300^3*10^-12)/12; %I_min = bh^3/12, h = 300 mm
syms P;

%Analytical solution for critical load
P_actual = pi^2*E*I/L^2;
err = 1;
N = 1;
P_vec = [];

%Increasing number of elements till error in P critial is below 10^-3
while err > 10^-3
    N = N+1; %Number of elements
    n = N+1; %Number of nodes
    l = L/N; %Length of each element
    
    %Constructing the global matrices
    ke_global = zeros(n,n);
    kp_global = sym(zeros(n,n));
    for i = 1:N
        [ke,kp] = buckling(P,l,E,I);
        ke_global(i:i+1,i:i+1) = ke_global(i:i+1,i:i+1) + ke;
        kp_global(i:i+1,i:i+1) = kp_global(i:i+1,i:i+1) + kp;
    end
    k_global = sym(ke_global) - kp_global;
    
    %Removing rows and columns for EBC
    k_global([1,end],:) = [];
    k_global(:,[1,end]) = [];
    
    %Solving the eigen value problem to find P_cr
    P_cr = vpa(min(solve(det(k_global) == 0)),5);
    P_vec(N-1) = P_cr;

    %Calculating the absolute error
    err = abs(P_cr-P_actual)/P_actual; 
end

fprintf("The converged value for the critical buckling load is %.3f kN \n",P_cr/1000)
fprintf("The analytical value for the critical buckling load is %.3f kN \n",P_actual/1000)
fprintf("The relative error is %f\n",err)

figure;
plot(2:1:N,P_vec,'-or',LineWidth=1.3)
xlabel("Number of elements")
ylabel("Critical Load")
title("Critical load vs Number of elements")

function [ke, kp] = buckling(P,L,E,I)
    ke = 1/L*[1, -1; -1, 1];
    kp = (P*L)/(E*I)*[1/3, 1/6; 1/6, 1/3];
end