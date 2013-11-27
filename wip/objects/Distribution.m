classdef Distribution
    % pdf, cdf, moments, stdnor, invcdf
    methods (Abstract)
        y=pdf(x);
        y=cdf(x);
        y=invcdf(x);
        y=moments(x);
    end
    methods
        function y=stdnor( obj, x )
            y=obj.invcdf( normal_cdf( x ) );
        end
    end
end