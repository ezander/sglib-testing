clear classes
clear variables
format compact
format short g

clc
top1 = TensorOperator({rand(5), rand(2)});
top2 = TensorOperator({rand(5), rand(2)});
x = rand(5,2);
size(top1)
size(top1,1)
size(top1,[],false)
size(top1,2,false)


top1*x
reshape(asmatrix(top1)*x(:), size(x)) - top1*x

top1*(top2*x) - (top1*top2)*x


top=(top1*top2);
top.asmatrix()



mop1 = MatrixOperator(rand(5))
mop2 = MatrixOperator(rand(5))
x = rand(5,2)
mop=(mop1*mop2);

y=mop*x
norm(mop2\(mop1\y)-x)
norm(mop\y-x)
norm(inv(mop)*y - x)
norm(inv(mop2)*inv(mop1)*y - x)

inv(mop.asmatrix()) - asmatrix(inv(mop))


clc
A = rand(5);
wsop = WoodburySolveOperator(MatrixOperator(A));


x = rand(5,2);
wsop*x - A*x

u = rand(5,2);
v = rand(5,2);
wsop=wsop.update(u,v);
norm(wsop * x - (A+u*v')*x)
norm(wsop\x - (A+u*v')\x)
