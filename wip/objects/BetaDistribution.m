classdef BetaDistribution < Distribution
    properties 
        a
        b
    end
    methods
        function obj=BetaDistribution(a,b)
            obj.a=a;
            obj.b=b;
        end
        
        function y=pdf(obj,x)
            % PDF Compute the probability distribution function of the Beta distribution.
            y=beta_pdf( x, obj.a, obj.b );
        end
        function y=cdf(obj,x)
            y=beta_cdf( x, obj.a, obj.b );
        end
        function x=invcdf(obj,y)
            x=beta_invcdf( y, obj.a, obj.b );
        end
        function [var,mean,skew,kurt]=moments(obj)
            [var,mean,skew,kurt]=beta_moments( obj.a, obj.b );
        end
        function y=stdnor(dist, x)
            y=dist.invcdf( normal_cdf( x ) );
        end
        
    end
end