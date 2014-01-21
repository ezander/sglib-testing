%%
% Lets look at the following boundary value problem 
%
% $$-\frac{\partial}{\partial x}\left(a(x) \frac{\partial}{\partial x}\right)=1 \textrm{ on } 
% \mathcal{D}=[0,1]$$
% 
% where $a$ is a coefficient field given by 
%
% $$a(x) = \begin{cases} a_1, & x<0.5 \\ a_2, & x\geq0.5 \end{cases} $$
%%
% We define 
N=11;
N=8;
[pos,els,bnd]=create_mesh_1d(0, 1, N);
k_func=@(x,k1, k2)(k1*double(x<0.5) + k2 * double(x>=0.5));

clc
format compact
format short g
k=funcall(k_func, pos, 10, 200)';
full(stiffness_matrix(pos, els, k))

full(stiffness_matrix(pos, els, funcall(k_func, pos, 1, 0)'))*2/7
full(stiffness_matrix(pos, els, funcall(k_func, pos, 0, 1)'))*2/7
