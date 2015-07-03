function test_nd_exact

clf

% set maximum number of stages for the integration rules
maxp = 6;
% set maximum polynomial degree to test the rules with
pt = 33;

for p=1:maxp
    subplot(3,maxp,p)
    test_rule(pt, p, @gauss_legendre_rule, 'smolyak', 'Smolyak/Legendre')
    
    subplot(3,maxp,p+maxp)
    %test_rule(pt, p, @gauss_legendre_rule, 'full_tensor', 'Tensor/Legendre')
    %test_rule(pt, p, @(n)(clenshaw_curtis_nested(n,'mode',1)), 'smolyak', 'Fejer1')
    test_rule(pt, p, @(n)(clenshaw_curtis_nested(n,'mode',2)), 'smolyak', 'Fejer2')
    
    subplot(3,maxp,p+2*maxp)
    test_rule(pt, p, @clenshaw_curtis_nested, 'smolyak', 'Smolyak/CC')
end


function test_rule(pt, p, rule_func, grid, title_str)
% Create all bivariate monomials up to order pt
V = gpcbasis_create('M', 'm', 2, 'p', pt, 'full_tensor', true);
I = V{2};
% Get integration rule with p stages and transform to [0,1]^2
[x,w] = integrate_nd([], rule_func, 2, p, 'grid', grid);
x = 0.5 * (x+1);
w = 0.5^2 * w;
% and integrate numerically.
q = gpcbasis_evaluate(V, x) * w;
% Compute the integral analytically x^m y^n integrates to 1/(m+1)/(n+1)
q_ex = 1./prod(I+1,2);
% We consider the rule is ok if the numerical and the analytical results
% differ by no more than 1e-10
ok = (abs(q-q_ex)<1e-10);
plot_exactness(I(ok,:), I(~ok,:), pt);
N = length(w);
title(sprintf('%s (p=%d, N=%d)', title_str, p, N));

function plot_exactness(I1, I2, pt)
cla;
edge_color = 'k';
plot_ind(I1, 'b', edge_color, pt);
plot_ind(I2, 'r', edge_color, pt);


function plot_ind(I, face_color, edge_color, pt)
% get coordinates from multiindex
x = I(:,1);
y = I(:,2);

% create arrays of rectangle vertices
h=1/2;
X=[x-h,x+h,x+h,x-h,x-h]';
Y=[y-h,y-h,y+h,y+h,y-h]';
patch( X, Y, face_color, 'EdgeColor', edge_color );

% make it look nice
axis( 'image' );
ylim([-h, pt+h]);
xlim([-h, pt+h]);
box( 'on' );
axis square
