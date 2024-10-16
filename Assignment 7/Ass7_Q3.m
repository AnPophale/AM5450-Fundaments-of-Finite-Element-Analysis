clc; clear all; close all; format compact; format shortg;

%Constants
% L = 0.015; %For 6 elements
% L = 0.005; %For 54 elements
L = 0.0025; %For 216 elements

a = L/2; b = L/2; %a and b for the square element
k = 45; %Thermal Conductivity
q0 = 8000; %Heat flux of left boundary
Q = 5*10^6; %Constant heat generation in the domain
h = 55; %Convective heat transfer coefficient
T0 = 110; %Essential BC at bottom boundary
Tamb = 20; %Ambient temperature
alpha1 = -h; %alpha for  convective NBC
beta1 = h*Tamb; %beta for convective NBC
beta2 = q0; %beta for constant heat flux NBC

%Importing Mesh data from excel
% connectivity = readmatrix("Q3_Mesh1.xlsx","Sheet",1); %For 6 elements
% connectivity = readmatrix("Q3_Mesh2.xlsx"); %For 54 elements
connectivity = readmatrix("Q3_Mesh3.xlsx"); %For 216 elements

% coords = readmatrix("Q3_Mesh1.xlsx","Sheet",2);
% coords = readmatrix("Q3_Mesh2.xlsx","Sheet",2);
coords = readmatrix("Q3_Mesh3.xlsx","Sheet",2);
X = coords(:,2); Y = coords(:,3);

elements = connectivity(:,1);
N = numel(elements); %Number of elements
n = size(coords,1); %Number of nodes

%Boundary nodes
% % Using 6 elements
% left_NBC_nodes = [7, 11, 6]; %Constant heat flux q0
% right_NBC_nodes1 = [3, 4];  %Insulated
% top_NBC_nodes1 = [6, 10, 5]; %Convective
% top_NBC_nodes2 = [1, 9, 4]; %Convective
% right_NBC_nodes2 = [1, 5]; %Convective
% EBC_nodes = [7, 12, 2, 8, 3]; %Fixed temperature

% % Using 54 elements
% left_NBC_nodes = [7, 33, 32, 31, 30, 29, 6]; % Constant heat flux q0
% right_NBC_nodes1 = [3, 15, 16, 4];           % Insulated
% top_NBC_nodes1 = [6, 28, 27, 26, 25, 24, 5]; % Convective
% top_NBC_nodes2 = [1, 21, 20, 19, 18, 17, 4]; % Convective
% right_NBC_nodes2 = [1, 22, 23, 5];           % Convective
% EBC_nodes = [7, 34, 35, 36, 37, 38, 2, 10, 11, 12, 13, 14, 3]; % Fixed temperature

%Using 216 elements
left_NBC_nodes = [7, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 6];  % Constant heat flux q0
right_NBC_nodes1 = [3, 24, 25, 26, 27, 28, 4];                        % Insulated
top_NBC_nodes1 = [6, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 5];  % Convective
top_NBC_nodes2 = [1, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 4];  % Convective
right_NBC_nodes2 = [1, 40, 41, 42, 43, 44, 5];                        % Convective
EBC_nodes = [7, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 2, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 3]; % Fixed temperature

%Assembling global k_k matrix and r_q vector
k_k_global = zeros(n,n);
r_q_global = zeros(n,1);

for i = 1:N
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4); node4 = connectivity(i,5);
    nodes = [node1,node2,node3,node4];
    k_k = stiffness_k(k,k,a,b);
    
    for j = 1:4
    k_k_global(nodes(j),nodes) = k_k_global(nodes(j),nodes) + k_k(j,:);
    r_q_global(nodes(j)) =  r_q_global(nodes(j)) + a*b*Q; 
    end
end

%Assembling k_alpha matrix and r_beta vector for NBCs
k_alpha_global = zeros(n,n);
r_beta_global = zeros(n,1);

%NBC on left boundary
%This method will work only when the boundary nodes are given sequentially,
%in the same order as they lie in the domain 
for i = 1:numel(left_NBC_nodes)-1
    node1 = left_NBC_nodes(i);
    node2 = left_NBC_nodes(i+1);

    r_beta = [b*beta2; b*beta2];
    r_beta_global([node1,node2]) = r_beta_global([node1,node2]) + r_beta;
end

%NBC on top boundary part 1
for i = 1:numel(top_NBC_nodes1)-1
    node1 = top_NBC_nodes1(i);
    node2 = top_NBC_nodes1(i+1);

    k_alpha = [-2*b*alpha1/3, -b*alpha1/3;
               -b*alpha1/3, -2*b*alpha1/3];
    k_alpha_global([node1,node2],[node1,node2]) = k_alpha_global([node1,node2],[node1,node2]) + k_alpha;

    r_beta = [a*beta1; a*beta1];
    r_beta_global([node1,node2]) = r_beta_global([node1,node2]) + r_beta;
end

%NBC for top boundary part 2
for i = 1:numel(top_NBC_nodes2)-1
    node1 = top_NBC_nodes2(i);
    node2 = top_NBC_nodes2(i+1);

    k_alpha = [-2*b*alpha1/3, -b*alpha1/3;
               -b*alpha1/3, -2*b*alpha1/3];
    
    k_alpha_global([node1,node2],[node1,node2]) = k_alpha_global([node1,node2],[node1,node2]) + k_alpha;

    r_beta = [a*beta1; a*beta1];
    r_beta_global([node1,node2]) = r_beta_global([node1,node2]) + r_beta;
end

%NBC for convective right side boundary
for i = 1:numel(right_NBC_nodes2)-1
    node1 = right_NBC_nodes2(i);
    node2 = right_NBC_nodes2(i+1);

    k_alpha = [-2*b*alpha1/3, -b*alpha1/3;
                -b*alpha1/3, -2*b*alpha1/3];    

    k_alpha_global([node1,node2],[node1,node2]) = k_alpha_global([node1,node2],[node1,node2]) + k_alpha;

    r_beta = [b*beta1; b*beta1];
    r_beta_global([node1,node2]) = r_beta_global([node1,node2]) + r_beta;
end

%Adding all the terms
k_global = k_k_global + k_alpha_global;
r_global = r_q_global + r_beta_global;

k_r = k_global;
r_r = r_global;

%Removing equations (rows) corresponding to EBCs
k_r(EBC_nodes,:) = [];
%Removing the rows corresponding to the EBCs
r_r(EBC_nodes) = []; 
%Adjusting the RHS using the EBCs
r_r = r_r - k_r(:,EBC_nodes)*T0*ones(numel(EBC_nodes),1);

%Removing the columns corresponding to the EBCs
k_r(:,EBC_nodes) = []; 

%Solving for the temperature field
T = zeros(n,1);
T_r = k_r\r_r;

%Combining the calculated T values and known T values from EBCs
idx = 1:1:n;
idx(EBC_nodes) = [];
T(idx) = T_r;
T(EBC_nodes) = T0;

%Plotting Temperature distributions
figure; hold on;
for i = 1:N
    %Using Fill to color each element with the temperature
    node1 = connectivity(i,2); node2 = connectivity(i,3); node3 = connectivity(i,4); node4 = connectivity(i,5);
    nodes = [node1,node2,node3,node4];
    fill(X(nodes),Y(nodes),T(nodes))
end
colormap("jet")
colorbar
scatter(X,Y,10,"filled",'ok')
xlabel("X (cm)"); ylabel("Y (cm)");
title("Temperature Distribution")
axis equal;
xlim([0,6]); ylim([0,3]);

%Plotting T at (3,3) with different number of elements
T_vec = [148.67,149.07,149.23];
N_vec = [6,54,216];

%Function to calculate k_k
function k_k = stiffness_k(kx,ky,a,b)
    k_k = [ky*a/(3*b) + kx*b/(3*a),   ky*a/(6*b) - kx*b/(3*a),    -ky*a/(6*b) - kx*b/(6*a),    kx*b/(6*a) - ky*a/(3*b);
          ky*a/(6*b) - kx*b/(3*a),    ky*a/(3*b) + kx*b/(3*a),    kx*b/(6*a) - ky*a/(3*b),    -ky*a/(6*b) - kx*b/(6*a);
         -ky*a/(6*b) - kx*b/(6*a),    kx*b/(6*a) - ky*a/(3*b),    ky*a/(3*b) + kx*b/(3*a),    ky*a/(6*b) - kx*b/(3*a);
          kx*b/(6*a) - ky*a/(3*b),    -ky*a/(6*b) - kx*b/(6*a),    ky*a/(6*b) - kx*b/(3*a),    ky*a/(3*b) + kx*b/(3*a)];
end