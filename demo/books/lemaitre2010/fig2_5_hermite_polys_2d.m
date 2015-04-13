function fig2_5_hermite_polys_2d(disp_1_to_6)
% Fig 2.5 page 33 (lemaitre & knio spectral methods)

V = gpcbasis_create('H', 'm', 2, 'p', 3);

N = 30;
x = linspace(-3, 3, N);
y = linspace(-3, 3, N);
[X, Y] = meshgrid(x, y);
z = gpcbasis_evaluate(V, [X(:), Y(:)]');

H = gpcbasis_polynomials(V, 'symbols', {'\xi_1', '\xi_2'});


if nargin<1 || disp_1_to_6
    ind=1:6;
else
    ind=7:10;
end

multiplot_init(length(ind)/2, 2, 'ordering', 'row');
set(gcf, 'defaulttextinterpreter', 'latex');
for i=ind
    multiplot;
    Z = reshape(z(i,:), size(X));
    mesh(X, Y, Z, 'EdgeColor', 'k', 'FaceColor', 'none');
    view([35, 30])
    xlabel('$\xi_1$');
    ylabel('$\xi_2$');
    xlim([-3,3]); 
    ylim([-3,3]); 
    title(['$\Psi_' num2str(i) '=' H{i} '$']);
end
