function hermite_p103
%HERMITE_P103 Generates the plots of 2D Hermite polynomials from page 103.
%   HERMITE_P103 Generates the plots of 2D Hermite polynomials from A.
%   Keese's diss. (page 103). 
%
% Example (<a href="matlab:run_example hermite_p103">run</a>)
%   hermite_p103
%
% See also HERMITE_VAL_MULTI, HERMITE

%   Elmar Zander
%   Copyright 2009, Institute of Scientific Computing, TU Braunschweig.
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.


% properties 
I=[ 1, 1; 2, 2; 2, 3; 3, 3];
%I=[ 3, 4; 2, 5; 7, 3; 6, 2];
V=gpcbasis_create('H', 'I', I);

clf
set( gcf, 'Renderer', 'zbuffer' );
if ~strcmp( get( gcf, 'WindowStyle' ), 'docked' )
    set( gcf, 'Position', [0, 0, 900, 900])
end
colormap( 'jet' );


range=2;
numpoints=50;
R=linspace(-range,range,numpoints);
[X,Y]=meshgrid( R );
ZI=gpcbasis_evaluate(V, [X(:), Y(:)]');

for i=1:gpcbasis_size(V,1)
    subplot(2,2,i);
    colormap('jet');
    shading('faceted');
    Z=reshape( ZI(i,:), size(X) );
    surf( X, Y, Z );
    title(sprintf('H_{(%d,%d)}(\\theta_1,\\theta_2)',I(i,1), I(i,2)));
end
