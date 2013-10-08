function test_deriv
% Experimentally test formulas for trace derivatives

N=20; M=30;
A=randn(N);
B=randn(N);
X=randn(N);
DX=1e-7*randn(N);

clc
format compact
format short g

test_frob(N, M)
test_trace_AXXB(A, B, X, DX)
test_trace_AX(A, X, DX)

function test_frob(N, M)
A = rand(N, M);
B = rand(N, M);
assert(abs(trace(A'*B) - frobenius_inner(A, B)) < 1e-10);
assert(abs(trace(A*B') - frobenius_inner(A, B)) < 1e-10);

function test_trace_AXXB(A, B, X, DX)
% Derivative of trace(A'*X'*X*B) is XBA'+XAB'
X2=X+DX;
DAXXB1 = trace( A'*X2'*X2*B )-trace( A'*X'*X*B )
D=X*B*A'+X*A*B';
DAXXB2 = frobenius_inner(D, DX)

function test_trace_AX(A, X, DX)
% derivative of trace(A*X) is A'
DAX1 = trace( A*(X+DX) ) - trace( A*X )
Y=A';
DAX2 = frobenius_inner(Y, DX)
