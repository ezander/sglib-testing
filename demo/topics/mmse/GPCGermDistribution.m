classdef GPCGermDistribution < Distribution
    % See also DISTRIBUTION 
    
    %   Elmar Zander
    %   Copyright 2016, Inst. of Scientific Computing, TU Braunschweig
    %
    %   This program is free software: you can redistribute it and/or
    %   modify it under the terms of the GNU General Public License as
    %   published by the Free Software Foundation, either version 3 of the
    %   License, or (at your option) any later version.
    %   See the GNU General Public License for more details. You should
    %   have received a copy of the GNU General Public License along with
    %   this program.  If not, see <http://www.gnu.org/licenses/>.
    
    properties (SetAccess=protected)
        V
        % The germ.
    end
    
    methods
        function dist=GPCGermDistribution(V)
            dist.V=V;
        end
        
        function y=pdf(dist,x)
            % PDF computes the probability distribution function of the
            % gpcgerm distribution.
            y = gpcgerm_pdf(dist.V, x);
        end
        
        function y=cdf(dist,x)
            % CDF computes the cumulative distribution function of the
            % gpcgerm distribution.
            y = gpcgerm_cdf(dist.V, x);
        end
        
        function x=invcdf(dist,y)
            % INVCDF computes the inverse CDF (or quantile) function of the
            % beta distribution.
            error('sglib:gpcgerm_noinvcdf', 'There is no invcdf for the distribution induced by the gpc germ');
        end
        
        function [mean,var,skew,kurt]=moments(dist)
            % MOMENTS computes the moments of the beta distribution.
            m = {nan, nan, nan, nan};
            [m{1:nargout}] = gpcgerm_moments( dist.V );
            [mean,var,skew,kurt]=deal(m{:});
        end
        
        function xi=sample(dist,n)
            %   Draw random samples from Normal distribution.
            %   XI=SAMPLE(DIST,N) draws N random samples from the random
            %   distribution DIST. If N is a scalar value XI is a column vector of
            %   random samples of size [N,1]. If N is a vector XI is a matrix (or
            %   tensor) of size [N(1), N(2), ...].
            xi = gpcgerm_sample(dist.V, n);
        end
        
        function str=tostring(dist)
            % Displays the distribution type
            str=strvarexpand('GPCGerm($dist.V$)');
        end
    end
    
end
