function S = rats2(X)
%RATS   Rational output.
%   RATS(X,LEN) uses RAT to display rational approximations to
%   the elements of X.  The string length for each element is LEN.
%   The default is LEN = 13, which allows 6 elements in 78 spaces.
%   Asterisks are used for elements which can't be printed in the
%   allotted space, but which are not negligible compared to the other
%   elements in X.
%
%   The same algorithm, with the default LEN, is used internally
%   by MATLAB for FORMAT RAT.
%
%   Class support for input X:
%      float: double, single
%
%   See also FORMAT, RAT.

% Modified version so that the output can be neatly used as input to
% unittests. Should be rewritten from scratch. Differences to the normal
% rats function: output is not centered, there are brackets around each
% line, and commas between numbers, semicolons at the end of each line, no
% unnecessary space, and not *** if the numbers are too large.


%   Copyright 1984-2009 The MathWorks, Inc. 
%   $Revision: 5.17.4.2 $  $Date: 2009/03/16 22:18:40 $

lens=10;
if nargin < 2, lens = 13; end
lhalf = (lens-1)/2;
tol = min(10^(-lhalf) * norm(X(isfinite(X)),1),.1);
[N,D] = rat(X,tol);
[m,n] = size(X);
nform = ['%d'];
dform = ['%d'];
S = [];
for i = 1:m
    s = [];
    for j = 1:n
        if D(i,j) ~= 1
            sj = [sprintf(nform,N(i,j)) '/' ...
                  sprintf(dform,D(i,j))];
        else
            sj = sprintf('%d',N(i,j));
        end
        s = [s ', ' sj];
    end
    s = ['[', s(3:end), '];'];
    disp(s);
    S(i,:) = rpad(s, 50);
end
S = char(S);


%-----------------------------
function t = lpad(s,len)
%LPAD Left pad with blanks.

t = [blanks(len-length(s)) s];


%-----------------------------
function t = rpad(s,len)
%RPAD Right pad with blanks.

t = [s blanks(len-length(s))];


%----------------------------
function t = cpad(s,len)
%CPAD Pad and center string with blanks.

padding = len-length(s);
t = [blanks(floor(padding/2)) s blanks(ceil(padding/2))];
