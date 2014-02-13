function smolyak_num_nodes

% initialise arrays for clenshaw curtis
ncc = [1, 1+2.^(1:20)];
del = [ncc(1), diff(ncc)];

% check the numbers in A. Keese diss
for d=[2,3]
    for stages=[4,7,9]
        k = stages - 1;
        q = k+d;
        n = smolyak_num_nodes1(q, d, del); %#ok<NASGU>
        strvarexpand('n(d=$d$,k=$k$) = $n$')
    end
end

% check table 1 in K. Petras
for k=3:4
    d=10;
    q = k+d;
    n = smolyak_num_nodes1(q, d, del); %#ok<NASGU>
    strvarexpand('n(d=$d$,k=$k$) = $n$')
    [x,w]=smolyak_grid(d, k+1, @clenshaw_curtis_nested); %#ok<ASGLU>
    n=length(w); %#ok<NASGU>
    strvarexpand('n(d=$d$,k=$k$) = $n$ (from grid)')
end



function n = smolyak_num_nodes1(q, d, del)
n1 = cumsum(del);
if d==1
    n = n1(q);
else
    n = 0;
    for nu = 0:q-d
        n = n + del(nu+1) * smolyak_num_nodes1(q-1-nu, d-1, del);
    end
end
