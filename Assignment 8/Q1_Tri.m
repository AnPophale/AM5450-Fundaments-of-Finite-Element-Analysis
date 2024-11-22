clc; clear all; close all; format compact; format shortg;

%Constants
Fs = 500*10^3; %Surface Shear Force
Ft = 600; %Tip load
E = 80*10^9; %Young's Modulus
nu = 0.3; %Poisson Ratio
L = 1; w = 0.1; %Dimensions of domain

%Importing Mesh data
connectivity = readmatrix("Q1_TriMesh.xlsx","Sheet",1);
coords = readmatrix("Q1_TriMesh.xlsx","Sheet",2);
X = coords(:,2); Y = coords(:,3);

%Boundary nodes
EBC_nodes = [1, 52, 103, 154, 205, 256];
Ft_node = 51;
Fs_nodes = [256; 257; 258; 259; 260; 261; 262; 263; 264; 265; 266; 267; 268; 269; 270; 271; 272; 273; 274; 275; ...
         276; 277; 278; 279; 280; 281; 282; 283; 284; 285; 286; 287; 288; 289; 290; 291; 292; 293; 294; 295; ...
         296; 297; 298; 299; 300; 301; 302; 303; 304; 305; 306];


N = size(connectivity,1); %Number of elements;
n = size(coords,1); %Number of nodes
A = L*w/N; %Area of each element

%Assembling global stiffness matrix
k_k_global = zeros(2*n,2*n); 
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4);
    nodes = [node1,node2,node3];
   
    x = X(nodes);
    y = Y(nodes);

    k_k = stiffness_k(x,y,A,E,nu);
    for j = 1:numel(nodes)
        k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes-1) = k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes-1) + k_k([2*j,2*j-1],[1,3,5]);
        k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes) = k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes) + k_k([2*j,2*j-1],[2,4,6]);
    end
end

%Shear Force boundary
r_q_global = zeros(2*n,1);
for i = 1:numel(Fs_nodes)-1
    node1 = Fs_nodes(i); node2 = Fs_nodes(i+1);
    x1 = X(node1); x2 = X(node2);
    y1 = Y(node1); y2 = Y(node2);
    L = ((x1-x2)^2 + (y1-y2)^2)^0.5;

    r_q_global(2*node1-1) = r_q_global(2*node1-1) + (L/2)*Fs; 
    r_q_global(2*node2-1) = r_q_global(2*node2-1) + (L/2)*Fs; 
end

%Applying EBC by removing corresponding equations
k_r = k_k_global;
k_r(2*EBC_nodes,:) = [];
k_r(2*EBC_nodes-1,:) = [];
k_r(:,2*EBC_nodes) = [];
k_r(:,2*EBC_nodes-1) = [];

r_r = r_q_global;
%Adding point load to r vector
r_r(2*Ft_node) = r_r(2*Ft_node) + Ft;

r_r(2*EBC_nodes) = [];
r_r(2*EBC_nodes-1) = [];

%Solving for d vector
d_r = k_r\r_r;

d = zeros(2*n,1);
idx = 1:1:2*n;
idx([2*EBC_nodes,2*EBC_nodes-1]) = [];
d(idx) = d_r;

%For plotting we are multiplying the displacememtns by 100 
%Otherwise the displacements are very small to be seen on the plot
X_new = X + 100*d(1:2:end);
Y_new = Y + 100*d(2:2:end);

%Plotting the displacement field
figure(1); 
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4);
    nodes = [node1,node2,node3];

    x = X_new(nodes); y = Y_new(nodes);
    subplot(2,1,1), fill(x,y,d(2*nodes-1)); hold on;
    subplot(2,1,2), fill(x,y,d(2*nodes)); hold on;
end
subplot(2,1,1), xlabel("X");  ylabel("Y"); title("U displacement"); colormap("jet"); colorbar
subplot(2,1,2), xlabel("X");  ylabel("Y"); title("V displacement"); colormap("jet"); colorbar

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

    e = strain(x,y,dl,A);
    e11(i) = e(1); e22(i) = e(2); gamma12(i) = e(3);
end

%Plotting strain field
figure;
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4);
    nodes = [node1,node2,node3];

    x = X_new(nodes); y = Y_new(nodes);
    subplot(3,1,1), fill(x,y,e11(i)*ones(3,1)); hold on;
    subplot(3,1,2), fill(x,y,e22(i)*ones(3,1)); hold on;
    subplot(3,1,3), fill(x,y,gamma12(i)*ones(3,1)); hold on;
end
subplot(3,1,1), xlabel("X");  ylabel("Y"); title("\epsilon_{11}"); colormap("jet"); colorbar
subplot(3,1,2), xlabel("X");  ylabel("Y"); title("\epsilon_{22}"); colormap("jet"); colorbar
subplot(3,1,3), xlabel("X");  ylabel("Y"); title("\gamma_{12}"); colormap("jet"); colorbar

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

    subplot(3,1,1), fill(x,y,sigma11(i)); hold on;
    subplot(3,1,2), fill(x,y,sigma22(i)); hold on;
    subplot(3,1,3), fill(x,y,tau12(i)); hold on;
end
subplot(3,1,1), xlabel("X");  ylabel("Y"); title("\sigma_{11}"); colormap("jet"); colorbar; 
subplot(3,1,2), xlabel("X");  ylabel("Y"); title("\sigma_{22}"); colormap("jet"); colorbar
subplot(3,1,3), xlabel("X");  ylabel("Y"); title("\tau_{12}"); colormap("jet"); colorbar

%Function to calculate stiffness matrix k_k
function k_k = stiffness_k(x,y,A,E,nu)
    x1 = x(1); x2 = x(2); x3 = x(3);
    y1 = y(1); y2 = y(2); y3 = y(3);

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

function e = strain(x,y,d,A)
    x1 = x(1); x2 = x(2); x3 = x(3);
    y1 = y(1); y2 = y(2); y3 = y(3);

    b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
    c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;

    B = (1/(2*A))*[b1, 0, b2, 0, b3, 0;
                0, c1, 0, c2, 0, c3;
                c1, b1, c2, b2, c3, b3]';
    
    e = B'*d;
end

