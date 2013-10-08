function [delta] = biegung_vektor(x,b1,h1,b2,winkel,E,F,L)
%mit winkel in GRAD
%F= Kraft, L= Laenge des Stabs, x= Punkt zur Berechnung, E= EModul, 
%I= Traegheitsmoment
%r= kruemmungsradius
%alle Angaben in METERN!!!!

[I]= traegheitsmoment_vektor(b1,h1,b2,winkel);
for i=1:size(x,2)
    delta(i,:)= (F*L*x(i)*x(i)*(3-x(i)/L))./ (6*E.*I);
    %ddw= - (2*F*x)./(3*E.*I) - (F*L*(x/L - 3))./(3*E.*I);  %2. Ableitung Biegelinie
    %r= (1+ delta.^2).^(3/2) ./ abs(ddw)
end

end