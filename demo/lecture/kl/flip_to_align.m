function f=flip_to_align(f, g)
% FLIP_TO_ALIGN Flips functions to be positively aligned

% Determine maximum number of functions in both F and G
n = min(size(f,2), size(g,2));

% Determines whether each of the 1 to n functions are positively aligned by
% computing their pairwise inner products <f_i, g_i>
factor=sign(sum(g(:,1:n).*f(:,1:n),1));

% Flips the not positively aligned functions
ind=factor==-1;
f(:,ind) = -f(:,ind);