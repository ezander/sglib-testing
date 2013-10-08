function [I] = traegheitsmoment_vektor(b1,h1,b2,winkel)
% mit Winkelvariation in Grad
    %A= zeros(1,3);
    %I_= zeros(1,3);
    %z= zeros(1,3);
    %z_S= zeros(1,3);
    %I_S= zeros(1,3);
    
    
    A1 = b1.* h1;
    A2 = b2.* b1;
    A3 = A1;
    
    I_1 = (A1.*(h1.^2))/12;
    I_2 = (A2.*(b1.^2))/12;
    I_3 = I_1;
   
    h = cosd(winkel).*h1;             % ist die tatsaechliche hoehe des balken (verrechnung des winkels)
    z0 = (h1+ b1)/2;                 %ist der schwerpunkt von flaechen 1,2, OHNE winkeleinfluss
    winkel = 270- winkel;
      
    z1 = h1+b1+ z0.*sind(winkel);
    z2 = b1./2 +h;
    z3 = z1;
    
    zS = (2*z1.*A1+ z2.*A2)./ (2*A1+ A2);
            
    z_S1 = abs(zS- z1);
    z_S2 = abs(zS- z2);
    z_S3 = abs(zS- z3);
    I_S1 = I_1+ (z_S1).^2.*A1;
    I_S2 = I_2+ (z_S2).^2.*A2;
    I_S3 = I_3+ (z_S3).^2.*A3;

    I = I_S1+ I_S2+ I_S3;
    %z = abs(zS- (h+b1)); %neutrale Faser
end