clc; clear all; close all; format compact; format shortg;

%Constants
E = 200*10^9;
A = 6400*10^-6;
I = 1.2*10^5*10^-9;
P = 10*10^3;
L = 10; 
M = 5*10^3;
EA = E*A; EI = E*I;

nodes = [0,0;0,L;L,L;L,0];
conn =[1,2;2,3;3,4];
n = size(nodes,1); %Number of nodes
N = size(conn,1); %Number of elements

%Assembling the global stiffness matrix
k_global = zeros(3*n,3*n);
for i = 1:N
    node1 = conn(i,1); node2 = conn(i,2);
    x1 = nodes(node1,1); x2 = nodes(node2,1);
    y1 = nodes(node1,2); y2 = nodes(node2,2);

    k = PlaneFrameElement(L,EA,EI,x1,y1,x2,y2);
    k_global(3*node1-2:3*node1, 3*node1-2:3*node1) = k_global(3*node1-2:3*node1, 3*node1-2:3*node1) + k(1:3,1:3);
    k_global(3*node1-2:3*node1, 3*node2-2:3*node2) = k_global(3*node1-2:3*node1, 3*node2-2:3*node2) + k(1:3,4:6);
    k_global(3*node2-2:3*node2, 3*node1-2:3*node1) = k_global(3*node2-2:3*node2, 3*node1-2:3*node1) + k(4:6,1:3);
    k_global(3*node2-2:3*node2, 3*node2-2:3*node2) = k_global(3*node2-2:3*node2, 3*node2-2:3*node2) + k(4:6,4:6);
end

%External force vector
F =  zeros(3*n,1);
F(4) = P;
F(9) = M;

%Removing rows and columns for EBCs
k_global_r = k_global;
k_global_r([1:3,end-2:end],:) = [];
k_global_r(:,[1:3,end-2:end]) = [];

F_r = F;
F_r([1:3,end-2:end]) = [];

%Solving for the nodal displacement and slopes
d = zeros(3*n,1);
d_r = k_global_r\F_r;
d(4:end-3) = d_r;

fprintf("Displacement and rotation at nodes 2 and 3:\n")
fprintf("Displacement of node 2 in X direction is %f m\n",d(4));
fprintf("Displacement of node 2 in Y direction is %f m\n",d(5));
fprintf("Rotation at node 2 is is %f\n",d(6));
fprintf("Displacement of node 3 in X direction is %f m\n",d(7));
fprintf("Displacement of node 3 in Y direction is %f m\n",d(8));
fprintf("Rotation of node 3 is %f\n\n",d(9));

%Forces and moments in element 1
node1 = conn(1,1); node2 = conn(1,2);
x1 = nodes(node1,1); x2 = nodes(node2,1);
y1 = nodes(node1,2); y2 = nodes(node2,2);
d1 = d(1:6);
[F1,M1] = ForceMoment(L,EA,EI,x1,y1,x2,y2,d1);
M1 = vpa(M1,5);

fprintf("Force and moment in element 1:\n")
fprintf("Axial force in element 1 is %f N\n",F1);
fprintf("Bending moment in element 1 in terms of s is %s Nm\n",char(M1));

function k = PlaneFrameElement(L,EA,EI,x1,y1,x2,y2)
    %Direction cosines for the element
    ls = (x2-x1)/L; ms=(y2-y1)/L;

    %Local stiffness matrix
    kl = [EA/L, 0, 0, -EA/L, 0, 0;
          0, 12*EI/L^3, 6*EI/L^2, 0, -12*EI/L^3, 6*EI/L^2;
          0, 6*EI/L^2, 4*EI/L, 0, -6*EI/L^2, 2*EI/L;
          -EA/L, 0, 0, EA/L, 0, 0;
          0, -12*EI/L^3, -6*EI/L^2, 0, 12*EI/L^3, -6*EI/L^2;
          0, 6*EI/L^2, 2*EI/L, 0, -6*EI/L^2, 4*EI/L];

    %Transformation matrix
    T = eye(6);
    T(1:2,1:2) = [ls,ms;-ms,ls];
    T(4:5,4:5) = [ls,ms;-ms,ls];
 
    %Global stiffness matrix
    k = T'*kl*T;
end

function [F,M] = ForceMoment(L,EA,EI,x1,y1,x2,y2,d)
    %Direction cosines for the element
    ls = (x2-x1)/L; ms=(y2-y1)/L;
    
    %Transformation matrix
    T = eye(6);
    T(1:2,1:2) = [ls,ms;-ms,ls];
    T(4:5,4:5) = [ls,ms;-ms,ls];
    
    %Local d vector
    dl = T*d;
    %Force 
    F = EA/L*(dl(4)-dl(1));
    %Moment
    syms s;
    M = EI*[-6+12*s, L*(-4+16*s), 6-12*s, L*(-2 + 6*s)]*[dl(2);dl(3);dl(5);dl(6)];
end