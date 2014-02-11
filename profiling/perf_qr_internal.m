n=1000;
k=400;
k2=300;
[K1,K2]=meshgrid( 100:100:2000, 100:100:1000 )
T=zeros(size(K1));
for j=1:numel(K1)
    fprintf( '%d/%d\n', j, numel(K1));
    k1=K1(j);
    k2=K2(j);
    %     n=100;
    %     k1=k1/10;
    %     k2=k1/10;
    
    A0=rand(n,k1);
    Q0=orth( A0 );
    k=size(Q0,2);
    korth=k;
    B=rand(n,k2);
    v=rand(k,1);
    
    QQ=rand(k,100);
    A=[Q0*diag(v) B];
    tic
    for xx=1:10
        [Q,R]=qr_internal( A, [], korth );
    end
    t1=toc;
    tic
    for xx=1:10
        [Q,R]=qr_internal( A, [], 0 );
    end
    t2=toc;
    T(j)=t2/t1;
    
    %assert_matrix( Q, {'unitary'}, 'qroc_orth_Q' );
    %assert_equals( A, Q*R, 'qroc_eq_A_QR' );
    
end
T
save qrfoobar.mat -ascii -double
