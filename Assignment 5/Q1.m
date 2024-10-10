clc; clear all; close all; format compact; format shortg;

%Constants
E1 = 200*10^9;
E2 = 70*10^9;
L = 2;
q = -10*10^3; %q is in downward direction (-ve Y direction)
F = -18*10^3; %F is in downward direction (-ve Y direction)
I = 4*10^-4;

%Stiffness matrices for each element
k1 = stiffness(E1,I,L);
k2 = stiffness(E1,I,L);
k3 = stiffness(E2,I,L);

%Global Stiffness Matrix
Kglobal = zeros(8,8);
Kglobal(1:4,1:4) = Kglobal(1:4,1:4) + k1;
Kglobal(3:6,3:6) = Kglobal(3:6,3:6) + k2;
Kglobal(5:8,5:8) = Kglobal(5:8,5:8) + k3;

%Global point and distributed load vectors
rF = zeros(8,1);
rF(3) = F;

rq = zeros(8,1);
rq(5:8) = [q*L/2; q*L^2/12;q*L/2; -q*L^2/12];

%Removing rows and columns with EBC (u and theta are known)
Kglobal_r = Kglobal;
Kglobal_r([1,2,7],:) = [];
Kglobal_r(:,[1,2,7]) = [];
rq_r = rq;
rq_r([1,2,7]) = [];
rF_r = rF;
rF_r([1,2,7]) = [];

%Solving for d vector
d_r = Kglobal_r\(rq_r + rF_r);
d = zeros(8,1);
d(3:6) = d_r(1:4);
d(end) = d_r(end);

fprintf("The deflection at the interface is %f m\n",d(5))
fprintf("The slope of the beam at the inteface is %f\n",d(6))

%Functon to calculate stiffness matrix of beam element
function k = stiffness(E,I,L)
    k = ((E*I)/(L^3))*[12, 6*L, -12, 6*L; 6*L, 4*L^2, -6*L, 2*L^2; -12, -6*L, 12, -6*L; 6*L, 2*L^2, -6*L, 4*L^2];
end
