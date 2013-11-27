clear
n=5;

H1=hermite(n,true);
H1=H1(:,end:-1:1)

[i,j,s]=find(spones(ones(n+1)));
s(:)=0;
fun=@(m)(mod(m+1,2).*factorial(m)./(2.^(m/2).*factorial(floor(m/2))));
fun(1:10);
s=fun(i+j-2);
hh=full(sparse(i,j,s));
S=diag(sqrt(factorial(0:n)));
L=chol(hh)';
H2=(L'\S)';

[I,J]=meshgrid(0:n,0:n);
mi=0:2*n;
moments=mod(mi+1,2).*factorial(mi)./(2.^(mi/2).*factorial(floor(mi/2)));
hh=moments(I+J+1);
H3=gram_schmidt( eye(n+1), hh );
H3=diag(1./diag(H3))*H3'

round(H1*hh*H1')
round(H2*hh*H2')
round(H3*hh*H3')
norm(H1-H2)
norm(H1-H3)
