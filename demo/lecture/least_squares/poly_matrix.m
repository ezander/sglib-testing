function A=poly_matrix(x, n, order)

if nargin<3
    order = -1;
end

x=x(:);
A=zeros(length(x), n+1);
for i=0:n
    if order<0
        A(:,n-i+1) = x .^ i;
    else
        A(:,i+1) = x .^ i;
    end
end
