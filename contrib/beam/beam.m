function [delta]= beam(p)

F= -10000;
L = 1;
x = 1;

h1 = p(1,:);
b1 = p(2,:);
b2 = p(3,:);
winkel = p(4,:);
E = p(5,:);

[delta] = biegung_vektor(x,b1,h1,b2,winkel,E,F,L);
