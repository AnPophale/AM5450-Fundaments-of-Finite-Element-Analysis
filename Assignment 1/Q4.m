clc; clear all; close all; format compact; format shortg;
 
%Constants
E = 200*10^9; d = 20*10^-3; A = pi*d^2/4; P = 2*10^3; L = 2; C = 10^3;

x = 0:0.01:L;

%Calculating displacement and strain fields
u = -(C/(6*E*A))*x.^3 + ((2*P+C*L^2)/(2*E*A))*x; 
u_prime = -(C/(2*E*A))*x.^2 + (2*P+C*L^2)/(2*E*A); 

%Plotting
figure;
plot(x,u,'-r')
title("u vs x")
xlabel("Length (m)")
ylabel('Displacement (m)')

figure;
plot(x,u_prime,'-b')
title("u' vs x")
xlabel("Length (m)")
ylabel('Strain')

%Changing the values of c
C_vec = 10^3*[1,3,5,7,9];

figure;
for i = 1:numel(C_vec)
    C = C_vec(i);
    u = -(C/(6*E*A))*x.^3 + ((2*P+C*L^2)/(2*E*A))*x;
    plot(x,u,DisplayName=num2str(C)); hold on       
end
title("u vs x")
xlabel("Length (m)")
ylabel('Displacement (m)')
legend("show")

figure;
for i = 1:numel(C_vec)
    C = C_vec(i);
    u_prime = -(C/(2*E*A))*x.^2 + (2*P+C*L^2)/(2*E*A); 
    plot(x,u_prime,DisplayName=num2str(C)); hold on
end
title("u' vs x")
xlabel("Length (m)")
ylabel('Strain')
legend("show")
