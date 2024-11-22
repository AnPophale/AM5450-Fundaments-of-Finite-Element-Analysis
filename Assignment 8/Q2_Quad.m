clc; clear all; close all; format compact; format shortg;

%Constants
P = 50*10^3; %Pressure in the vessel
E = 200*10^9; %Young's Modulus
nu = 0.3; %Poisson Ratio

%Importing Mesh data
connectivity = readmatrix("Q2_QuadMesh.xlsx","Sheet",1);
coords = readmatrix("Q2_QuadMesh.xlsx","Sheet",2);
X = coords(:,2); Y = coords(:,3);

%Boundary nodes
Xsym_nodes = [4, 56, 55, 54, 53, 52, 5];
Ysym_nodes = [1, 7, 8, 9, 10, 11, 2];
P_vertical_nodes = [1, 18, 17, 16, 15, 14, 13, 12, 3];
P_horizontal_nodes =  [4, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, ...
         41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 3];

N = size(connectivity,1); %Number of elements;
n = size(coords,1); %Number of nodes

%Assembling global stiffness matrix
k_k_global = zeros(2*n,2*n); 
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4); node4 = connectivity(i,5);
    nodes = [node1,node2,node3,node4];

    x = X(nodes);
    y = Y(nodes);

    k_k = stiffness_k(x,y,E,nu);
    for j = 1:numel(nodes)
        k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes-1) = k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes-1) + k_k([2*j,2*j-1],[1,3,5,7]);
        k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes) = k_k_global([2*nodes(j),2*nodes(j)-1],2*nodes) + k_k([2*j,2*j-1],[2,4,6,8]);
    end
end

%Horizontal Pressure Boundary
r_p_global = zeros(2*n,1);
for i = 1:numel(P_horizontal_nodes)-1
    node1 = P_horizontal_nodes(i); node2 = P_horizontal_nodes(i+1);
    x1 = X(node1); x2 = X(node2);
    y1 = Y(node1); y2 = Y(node2);
    L = ((x1-x2)^2 + (y1-y2)^2)^0.5;

    %LCS at boundary
    syms t;
    Nc = 0.5*[1-t; 1+t];
    xc = Nc'*[x1;x2]; yc = Nc'*[y1;y2];
    Jc = L/2;
    
    %Calculating rq using 2 point Gauss Quadrature
    r_q = P*Jc*2*pi*(subs(xc*Nc,t,0.577) + subs(xc*Nc,t,-0.577));

    r_p_global(2*node1-1) = r_p_global(2*node1-1) + r_q(1); 
    r_p_global(2*node2-1) = r_p_global(2*node2-1) + r_q(2); 
end

%Vertical Pressure Boundary
for i = 1:numel(P_vertical_nodes)-1
    node1 = P_vertical_nodes(i); node2 = P_vertical_nodes(i+1);
    x1 = X(node1); x2 = X(node2);
    y1 = Y(node1); y2 = Y(node2);
    L = ((x1-x2)^2 + (y1-y2)^2)^0.5;
    
    %LCS at boundary
    syms s;
    Nc = 0.5*[1-s; 1+s];
    xc = Nc'*[x1;x2]; yc = Nc'*[y1;y2];
    Jc = L/2;
    
    %Calculating rq using 2 point Gauss Quadrature
    r_q = P*Jc*2*pi*(subs(xc*Nc,s,0.577) + subs(xc*Nc,s ,-0.577));

    r_p_global(2*node1) = r_p_global(2*node1) + r_q(1); 
    r_p_global(2*node2) = r_p_global(2*node2) + r_q(2); 
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
X_new = X + 10^5*d(1:2:end);
Y_new = Y + 10^5*d(2:2:end);

%Plotting the displacement field
figure(1); 
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4); node4 = connectivity(i,5);
    nodes = [node1,node2,node3,node4];

    x = X_new(nodes); y = Y_new(nodes);
    subplot(1,2,1), fill(x,y,d(2*nodes-1)); hold on;
    subplot(1,2,2), fill(x,y,d(2*nodes)); hold on;
    
end
subplot(1,2,1), xlabel("X");  ylabel("Y"); title("U displacement"); colormap("jet"); colorbar;   
subplot(1,2,2), xlabel("X");  ylabel("Y"); title("V displacement"); colormap("jet"); colorbar; 

%Strain Calculation
e11 = zeros(N,1); e22 = zeros(N,1); e33 = zeros(N,1); gamma12 = zeros(N,1);
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4); node4 = connectivity(i,5);
    nodes = [node1,node2,node3,node4];
    x = X(nodes);
    y = Y(nodes);

    dl = zeros(8,1);
    dl([1,3,5,7]) = d(2*nodes-1);
    dl([2,4,6,8]) = d(2*nodes);
    
    %Calculating strain at the centroid (0,0)
    e = strain(x,y,dl,0,0);
    e11(i) = e(1); e22(i) = e(2); e33(i) = e(3); gamma12(i) = e(4);
end

%Plotting strain field
figure;
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4); node4 = connectivity(i,5);
    nodes = [node1,node2,node3,node4];

    x = X_new(nodes); y = Y_new(nodes);
    subplot(1,4,1), fill(x,y,e11(i)); hold on;
    subplot(1,4,2), fill(x,y,e22(i)); hold on;
    subplot(1,4,3), fill(x,y,e33(i)); hold on;
    subplot(1,4,4), fill(x,y,gamma12(i)); hold on;
end
subplot(1,4,1), xlabel("X");  ylabel("Y"); title("\epsilon_{11}"); colormap("jet"); colorbar
subplot(1,4,2), xlabel("X");  ylabel("Y"); title("\epsilon_{22}"); colormap("jet"); colorbar
subplot(1,4,3), xlabel("X");  ylabel("Y"); title("\epsilon_{33}"); colormap("jet"); colorbar
subplot(1,4,4), xlabel("X");  ylabel("Y"); title("\gamma_{12}"); colormap("jet"); colorbar

%Calculating and plotting stress field
sigma11 = zeros(N,1); sigma22 = zeros(N,1); sigma33 = zeros(N,1); tau12 = zeros(N,1);
figure;
for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4); node4 = connectivity(i,5);
    nodes = [node1,node2,node3,node4];
    x = X_new(nodes); y = Y_new(nodes);

    %C for axisymmetric case
    C = (E/((1+nu))*(1-2*nu))*[1-nu, nu, nu, 0;
                      nu, 1- nu, nu, 0;
                      nu, nu, 1-nu, 0;
                      0, 0, 0, (1-2*nu)/2];

    stress = C*[e11(i);e22(i);e33(i);gamma12(i)];
    sigma11(i) = stress(1); sigma22(i) = stress(2); sigma33(i) = stress(3); tau12(i) = stress(3);

    subplot(1,4,1), fill(x,y,sigma11(i)); hold on;
    subplot(1,4,2), fill(x,y,sigma22(i)); hold on;
    subplot(1,4,3), fill(x,y,sigma33(i)); hold on
    subplot(1,4,4), fill(x,y,tau12(i)); hold on;
end
subplot(1,4,1), xlabel("X");  ylabel("Y"); title("\sigma_{11}"); colormap("jet"); colorbar
subplot(1,4,2), xlabel("X");  ylabel("Y"); title("\sigma_{22}"); colormap("jet"); colorbar
subplot(1,4,3), xlabel("X");  ylabel("Y"); title("\sigma_{33}"); colormap("jet"); colorbar
subplot(1,4,4), xlabel("X");  ylabel("Y"); title("\tau_{12}"); colormap("jet"); colorbar

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
function k_k = stiffness_k(r_vec,z_vec,E,nu)
    syms s t;
    %Shape Functions
    N = 0.25*[(s-1)*(t-1), -(s+1)*(t-1), (s+1)*(t+1), -(s-1)*(t+1)];
    r = N*r_vec;

    %Derivative wrt to LCS
    dNdS = diff(N,s,1);
    dNdT = diff(N,t,1);
    
    %Calculating components of the Jacobian
    drdS = dNdS*r_vec;
    drdT = dNdT*r_vec;
    dzdS = dNdS*z_vec;
    dzdT = dNdT*z_vec;

    detJ = drdS*dzdT - drdT*dzdS;

    %B operator
    Br = (1/detJ)*(dzdT*dNdS - drdT*dNdT);
    Bz = (1/detJ)*(-dzdS*dNdS + drdS*dNdT);

    B = [Br(1), 0, Br(2), 0, Br(3), 0, Br(4), 0;
        0, Bz(1), 0, Bz(2), 0, Bz(3), 0, Bz(4);
        N(1)/r, 0, N(2)/r, 0, N(3)/r, 0, N(4)/r, 0;
        Bz(1), Br(1), Bz(2), Br(2), Bz(3), Br(3), Bz(4), Br(4)];

    B = B';
    
    %C matrix  for axisymmetric condition
    C = (E/((1+nu))*(1-2*nu))*[1-nu, nu, nu, 0;
                      nu, 1- nu, nu, 0;
                      nu, nu, 1-nu, 0;
                      0, 0, 0, (1-2*nu)/2];

   %Calculating stiffness matrix 
   %Using 2 point Gauss Quadrature
   s_vec = [0.577,-0.577];
   t_vec = [0.577,-0.577];
   k_k = zeros(8,8);
   for i = 1:2
       for j = 1:2
           k_k = k_k + subs((B*C*B'*detJ*2*pi*r),[s,t],[s_vec(i),t_vec(j)]);
       end
   end

end

function e = strain(r_vec,z_vec,d,s_val,t_val)
    syms s t;
    %Shape Functions
    N = 0.25*[(s-1)*(t-1), -(s+1)*(t-1), (s+1)*(t+1), -(s-1)*(t+1)];
    r = N*r_vec;

    %Derivative wrt to LCS
    dNdS = diff(N,s,1);
    dNdT = diff(N,t,1);
    
    %Calculating components of the Jacobian
    drdS = dNdS*r_vec;
    drdT = dNdT*r_vec;
    dzdS = dNdS*z_vec;
    dzdT = dNdT*z_vec;

    detJ = drdS*dzdT - drdT*dzdS;

    %B operator 
    Br = (1/detJ)*(dzdT*dNdS - drdT*dNdT);
    Bz = (1/detJ)*(-dzdS*dNdS + drdS*dNdT);

    B = [Br(1), 0, Br(2), 0, Br(3), 0, Br(4), 0;
        0, Bz(1), 0, Bz(2), 0, Bz(3), 0, Bz(4);
        N(1)/r, 0, N(2)/r, 0, N(3)/r, 0, N(4)/r, 0;
        Bz(1), Br(1), Bz(2), Br(2), Bz(3), Br(3), Bz(4), Br(4)];

    B = B';
    
    %Calculating strain at the given s,t in the element
    e = subs(B,[s,t],[s_val,t_val])'*d;
end

