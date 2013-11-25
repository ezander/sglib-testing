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
            y=beta_pdf( x, obj.a, obj.b );
        end
        function y=cdf(obj,x)
            y=beta_cdf( x, obj.a, obj.b );
        end
        function y=invcdf(obj,x)
            y=beta_invcdf( x, obj.a, obj.b );
        end
        function [var,mean,skew,kurt]=moments(obj,x)
            [var,mean,skew,kurt]=beta_moments( obj.a, obj.b );
        end
    end
end