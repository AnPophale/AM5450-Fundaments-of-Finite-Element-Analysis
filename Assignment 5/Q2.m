clc; clear all; close all; format compact; format shortg;

%Constants
L = 0.5; 
E = 200*10^7    ; 
d = 20*10^-3; A = (pi*d^2)/4;
F1 = 10*10^3; F2 = 12*10^3;

nodes = [0,0; 0,L; L,0; L,L; 2*L,0; 2*L,L]; %Node coordinates
conn = [1,2; 1,3; 2,3; 2,4; 1,4; 3,4; 3,5; 4,5; 3,6; 4,6; 5,6]; %Connectivity
N = size(conn,1); %Number of elements;
n = size(nodes,1); %Number of nodes;

Kglobal = zeros(2*n,2*n);
%Constructing the global stiffness matrix
for i = 1:N
    node1 = conn(i,1);
    node2 = conn(i,2);
    
    x1 = nodes(node1,1); y1 = nodes(node1,2);
    x2 = nodes(node2,1); y2 = nodes(node2,2);

    k = stiffness(x1,y1,x2,y2,E,A);
    
    Kglobal(2*node1-1:2*node1,2*node1-1:2*node1) = Kglobal(2*node1-1:2*node1,2*node1-1:2*node1) + k(1:2,1:2);
    Kglobal(2*node1-1:2*node1,2*node2-1:2*node2) = Kglobal(2*node1-1:2*node1,2*node2-1:2*node2) + k(1:2,3:4);
    Kglobal(2*node2-1:2*node2,2*node1-1:2*node1) = Kglobal(2*node2-1:2*node2,2*node1-1:2*node1) + k(3:4,1:2);
    Kglobal(2*node2-1:2*node2,2*node2-1:2*node2) = Kglobal(2*node2-1:2*node2,2*node2-1:2*node2) + k(3:4,3:4);
end

%Construncting the global force vector
Fglobal = zeros(2*n,1);
Fglobal(4) = -F1;
Fglobal(8) = -F2;
Fglobal(11) = F1;

%Removing part of the matrices which are not solved for (EBC)
Kglobal_r = Kglobal; Fglobal_r = Fglobal;
Kglobal_r([1,2,10],:) = [];
Kglobal_r(:,[1,2,10]) = [];
Fglobal_r([1,2,10]) = [];

%Solving for the displacements
d_r = Kglobal_r\Fglobal_r;
%Adding the known u values from EBC to the solution vector
d = zeros(n,1);
d(3:9) = d_r(1:7);
d(11:12) = d_r(8:9);

u = zeros(n,1); v = zeros(n,1);
for i = 1:n
    u(i) = d(2*i-1);
    v(i) = d(2*i);
end

fprintf("Nodal Displacements:\n")
T1 = table((1:1:n)',u,v,VariableNames=["Node","u (m)","v (m)"]);
disp(T1);

%Calculating the reactions at supports
R1x = Kglobal(1,:)*d; %Reaction at node 1 in X direction 
R1y = Kglobal(2,:)*d; %Reaction at node 1 in Y direction
R5y = Kglobal(10,:)*d; %Reaction at node 5 in Y direction

fprintf("\nReaction Forces:\n")
fprintf("Reaction force in X direction at node 1 = %f N\n",R1x)
fprintf("Reaction force in Y direction at node 1 = %f N \n",R1y)
fprintf("Reaction force in Y direction at node 5 = %f N \n\n",R5y)

%New node positions
nodes_new = zeros(n,2);
for i = 1:n
    nodes_new(i,1) = nodes(i,1) + d(2*i-1);
    nodes_new(i,2) = nodes(i,2) + d(2*i);
end

%Calculating stress and force in each element
stress_vec = zeros(N,1);
force_vec = zeros(N,1);
for i = 1:N
    node1 = conn(i,1);
    node2 = conn(i,2);

    x1 = nodes(node1,1); y1 = nodes(node1,2);
    x2 = nodes(node2,1); y2 = nodes(node2,2);
    
    u1 = d(2*node1-1); v1 = d(2*node1);
    u2 = d(2*node2-1); v2 = d(2*node2);

    stress_vec(i) = element_stress(x1,y1,x2,y2,u1,v1,u2,v2,E);
    force_vec(i) = A*stress_vec(i);
end

fprintf("Elemental Stress and Forces:\n")
T2 = table((1:1:N)',stress_vec,force_vec,VariableNames=["Element","Stress (Pa)","Force (N)"]);
disp(T2)

%Plotting the old and new truss structure
figure;
hold on;
for i = 1:N
    node1 = conn(i,1);
    node2 = conn(i,2);

    x1 = nodes(node1,1); y1 = nodes(node1,2);
    x2 = nodes(node2,1); y2 = nodes(node2,2);

    plot([x1,x2],[y1,y2],'-.r',LineWidth=1.5)

    x1_new = nodes_new(node1,1); y1_new = nodes_new(node1,2);
    x2_new = nodes_new(node2,1); y2_new = nodes_new(node2,2);
    
    plot([x1_new,x2_new],[y1_new,y2_new],'--.b',LineWidth=1.5)
end
title("Undeformed (Red) and Deformed (Blue) Truss Structure")
xlabel("x (m)")
ylabel("y (m)")

%function to calculate the stiffness matrix in GCS
function k = stiffness(x1,y1,x2,y2,E,A)
    l = sqrt((x1-x2)^2 + (y1-y2)^2);
    c = (x2-x1)/l;
    s = (y2-y1)/l;

    T = [c,s,0,0; 0,0,c,s]; %Transformation matrix
    kl = (E*A/l)*[1,-1; -1,1];
    k = T'*kl*T;
end

%function to calculate stress and strain in each element
function s = element_stress(x1,y1,x2,y2,u1,v1,u2,v2,E)
    l = sqrt((x1-x2)^2 + (y1-y2)^2);
    c = (x2-x1)/l;
    s = (y2-y1)/l;

    T = [c,s,0,0; 0,0,c,s]; %Transformation matrix
    d = [u1;v1;u2;v2];
    dl = T*d;
    s = E*(dl(2)-dl(1))/l;
end
