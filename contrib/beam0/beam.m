function delta = beam(p)

h1 = p(1,:);
b1 = p(2,:);
b2 = p(3,:);

delta = 0.01 * h1 .* b2 .^ 2 + 300 * b1 .^ 4;
