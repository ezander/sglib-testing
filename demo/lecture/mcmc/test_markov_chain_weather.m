% Markov chain demo

% Prob of going from bad weather to good weather
pGB = 0.3;
pBB = 0.7;

% Prob of going from good weather to bad weather
pGG = 0.5;
pBG = 0.5;

% States: 0 == bad, 1 good
% Now going along one sample path and computing the mean at the end
state = 0;
states = [];
N =10000;
for i=1:N
    r = rand(1);
    if state==0
        if r<pGB
            state = 1;
        end
    elseif state==1
        if r<pBG
            state = 0;
        end
    end
    %fprintf('%g', state);
    states(end+1) = state;
end
strvarexpand('Prob for good weather: $mean(states)$');

% Computing the equilibrium probability vector using matrix iterations
pB = 0.5;
pG = 0.5;
p = [pB; pG];
% pBn = pBB * pB + pBG * pG
% pGn = pGB * pB + pGG * pG
M = [pBB pBG; pGB pGG];
for i=1:100
    p = M * p;
end
strvarexpand('Prob for good weather: $p(2)$');

% Computing the equilibrium probability vector using eigenvalue decomp
[U,D]=eig(M);
u = U(:,1);
u = u / sum(u);
strvarexpand('Prob for good weather: $u(2)$');

% Computing the equilibrium probability vector using kernel of M-I
u=null(M-eye(size(M)));
u = u / sum(u);
strvarexpand('Prob for good weather: $u(2)$');


