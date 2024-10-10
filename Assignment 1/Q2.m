clc; clear all; close all; format compact; format shortg;

syms x a0 a1 a2 a3 a4 a5
    
%Degree 5 polynomial
p = poly2sym([a5,a4,a3,a2,a1,a0])

%1st to 5th order derivatives
d1 = diff(p,x,1)
d2 = diff(p,x,2)
d3 = diff(p,x,3)
d4 = diff(p,x,4)
d5 = diff(p,x,5)

%Partial derivatives wrt to ai
da0 = diff(p,a1,1)
da1 = diff(p,a1,1)
da2 = diff(p,a2,1)
da3 = diff(p,a3,1)
da4 = diff(p,a4,1)
da5 = diff(p,a5,1)