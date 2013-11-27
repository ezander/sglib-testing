function C=structure_coefficients( p )

if ~exist('p')
    p=orthpoly_from_moments( {@normal_raw_moments, {0, 1}, {1, 2}}, ...
       5, 'monic', true );
end

% bring p into some sensible form (i.e. polynomials are column vectors in
% ascending order of coefficients).
q=fliplr(p');

C=zeros(n,n,n);
n=size(p,1);
for i=1:n
    for j=1:n
        p_ij=conv( p(i,:), p(j,:));
        %fprintf( '(%s) * (%s) = %s\n', format_poly(p(i,:)), format_poly(p(j,:)), format_poly(conv( p(i,:), p(j,:)) ) );
    end
end
