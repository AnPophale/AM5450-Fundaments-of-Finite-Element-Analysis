clc; clear all; close all; format compact; format shortg;

%Constants
P = 50*10^3; %Pressure in the vessel
E = 200*10^9; %Young's Modulus
nu = 0.3; %Poisson Ratio

%Importing Mesh data
connectivity = readmatrix("Q2_Mesh.xlsx","Sheet",1);
coords = readmatrix("Q2_Mesh.xlsx","Sheet",2);
X = coords(:,2); Y = coords(:,3);

%Boundary nodes
Xsym_nodes = [4, 36, 35, 34, 5];
Ysym_nodes = [1, 7, 8, 9, 2];
P_vertical_nodes = [1, 13, 12, 11, 10, 3];
P_horizontal_nodes = [4, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 3];

N = size(connectivity,1); %Number of elements;
n = size(coords,1); %Number of nodes

%Assembling global stiffness matrix
k_k_global = zeros(2*n,2*n); 
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4);
    nodes = [node1,node2,node3];
   
    x = X(nodes);
    y = Y(nodes);

    k_k = stiffness_k(x,y,E,nu);
    for j = 1:numel(nodes)
        k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes-1) = k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes-1) + k_k([2*j,2*j-1],[1,3,5]);
        k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes) = k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes) + k_k([2*j,2*j-1],[2,4,6]);
    end
end

%Horizontal Pressure Boundary
r_p_global = zeros(2*n,1);
for i = 1:numel(P_horizontal_nodes)-1
    node1 = P_horizontal_nodes(i); node2 = P_horizontal_nodes(i+1);
    x1 = X(node1); x2 = X(node2);
    y1 = Y(node1); y2 = Y(node2);
    L = ((x1-x2)^2 + (y1-y2)^2)^0.5;

    r_p_global(2*node1-1) = r_p_global(2*node1-1) + (L/2)*P; 
    r_p_global(2*node2-1) = r_p_global(2*node2-1) + (L/2)*P; 
end

%Vertical Pressure Boundary
for i = 1:numel(P_vertical_nodes)-1
    node1 = P_vertical_nodes(i); node2 = P_vertical_nodes(i+1);
    x1 = X(node1); x2 = X(node2);
    y1 = Y(node1); y2 = Y(node2);
    L = ((x1-x2)^2 + (y1-y2)^2)^0.5;

    r_p_global(2*node1) = r_p_global(2*node1) - (L/2); 
    r_p_global(2*node2) = r_p_global(2*node2) - (L/2); 
end

%Applying symmetry conditions by removing corresponding equations
k_r = k_k_global;
k_r([2*Ysym_nodes-1,2*Xsym_nodes],:) = [];
k_r(:,[2*Ysym_nodes-1,2*Xsym_nodes]) = [];

r_r = r_p_global;
r_r([2*Ysym_nodes-1,2*Xsym_nodes]) = [];

%Solving for d vector
d_r = k_r\r_r;

d = zeros(2*n,1);
idx = 1:1:2*n;
idx([2*Ysym_nodes-1,2*Xsym_nodes]) = [];
d(idx) = d_r;

%For plotting we are multiplying the displacememtns by 1000
%Otherwise the displacements are very small to be seen on the plot
X_new = X + 1000*d(1:2:end);
Y_new = Y + 1000*d(2:2:end);

%Plotting the displacement field
figure(1); 
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4);
    nodes = [node1,node2,node3];

    x = X_new(nodes); y = Y_new(nodes);
    subplot(1,2,1), fill(x,y,d(2*nodes-1)); hold on;
    subplot(1,2,2), fill(x,y,d(2*nodes)); hold on;
    
end
subplot(1,2,1), xlabel("X");  ylabel("Y"); title("U displacement"); colormap("jet"); colorbar;   
subplot(1,2,2), xlabel("X");  ylabel("Y"); title("V displacement"); colormap("jet"); colorbar; 

%Strain Calculation
e11 = zeros(N,1); e22 = zeros(N,1); gamma12 = zeros(N,1);
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4);
    nodes = [node1,node2,node3];
    x = X(nodes);
    y = Y(nodes);

    dl = zeros(6,1);
    dl([1,3,5]) = d(2*nodes-1);
    dl([2,4,6]) = d(2*nodes);

    e = strain(x,y,dl);
    e11(i) = e(1); e22(i) = e(2); gamma12(i) = e(3);
end

%Plotting strain field
figure;
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4);
    nodes = [node1,node2,node3];

    x = X_new(nodes); y = Y_new(nodes);
    subplot(1,3,1), fill(x,y,e11(i)); hold on;
    subplot(1,3,2), fill(x,y,e22(i)); hold on;
    subplot(1,3,3), fill(x,y,gamma12(i)); hold on;
end
subplot(1,3,1), xlabel("X");  ylabel("Y"); title("\epsilon_{11}"); colormap("jet"); colorbar
subplot(1,3,2), xlabel("X");  ylabel("Y"); title("\epsilon_{22}"); colormap("jet"); colorbar
subplot(1,3,3), xlabel("X");  ylabel("Y"); title("\gamma_{12}"); colormap("jet"); colorbar

%Calculating and plotting stress field
sigma11 = zeros(N,1); sigma22 = zeros(N,1); tau12 = zeros(N,1);
figure;
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4);
    nodes = [node1,node2,node3];
    x = X_new(nodes); y = Y_new(nodes);

    C = (E/(1-nu^2))*[1, nu, 0;
                     nu, 1, 0;
                     0, 0 , (1-nu)/2];
    stress = C*[e11(i);e22(i);gamma12(i)];
    sigma11(i) = stress(1); sigma22(i) = stress(2); tau12(i) = stress(3);

    subplot(1,3,1), fill(x,y,sigma11(i)); hold on;
    subplot(1,3,2), fill(x,y,sigma22(i)); hold on;
    subplot(1,3,3), fill(x,y,tau12(i)); hold on;
end
subplot(1,3,1), xlabel("X");  ylabel("Y"); title("\sigma_{11}"); colormap("jet"); colorbar
subplot(1,3,2), xlabel("X");  ylabel("Y"); title("\sigma_{22}"); colormap("jet"); colorbar
subplot(1,3,3), xlabel("X");  ylabel("Y"); title("\tau_{12}"); colormap("jet"); colorbar

hoop = max(abs(sigma11));
long = max(abs(sigma22));

d = 0.4; t = 0.075;
hoop_thin = P*d/(2*t);
long_thin = P*d/(4*t);

fprintf("The hoop stress calculated is %f N\n",hoop)
fprintf("The longitudinal stress calculated is %f N\n", long);
fprintf("The analytical hoop stress assuming thin walled vessel is %f N\n",hoop_thin)
fprintf("The analytical longitudinal stress assuming thin walled vessel is  %f N\n", long_thin);

%Function to calculate stiffness matrix k_k
function k_k = stiffness_k(x,y,E,nu)
    x1 = x(1); x2 = x(2); x3 = x(3);
    y1 = y(1); y2 = y(2); y3 = y(3);

    f1 = x2*y3-x3*y2; f2 = x3*y1-x1*y3; f3 = x1*y2-x2*y1;
    A = 0.5*(f1+f2+f3);

    b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
    c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;

    B = (1/(2*A))*[b1, 0, b2, 0, b3, 0;
                0, c1, 0, c2, 0, c3;
                c1, b1, c2, b2, c3, b3]';
    
    C = (E/(1-nu^2))*[1, nu, 0;
                      nu, 1, 0;
                      0, 0 , (1-nu)/2];

    k_k = A*(B*C*B');
end

function e = strain(x,y,d)
    x1 = x(1); x2 = x(2); x3 = x(3);
    y1 = y(1); y2 = y(2); y3 = y(3);

    f1 = x2*y3-x3*y2; f2 = x3*y1-x1*y3; f3 = x1*y2-x2*y1;
    A = 0.5*(f1+f2+f3);

    b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
    c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;

    B = (1/(2*A))*[b1, 0, b2, 0, b3, 0;
                0, c1, 0, c2, 0, c3;
                c1, b1, c2, b2, c3, b3]';
    
    e = B'*d;
end

