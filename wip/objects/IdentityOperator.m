classdef IdentityOperator
    % IDENTITYOPERATOR Efficient implementation of the identity operator.
    properties
        n
    end
    methods (Abstract)
        function idop=IdentityOperator(n)
            idop.n = n(:);
        end
        
        function x=apply(idop, x) %#ok<INUSL>
        end
        
        function x=solve(idop, x) %#ok<INUSL>
        end
        
        function s=size_impl(op)
            s = [op.n, op.n];
        end
        
        function A=asmatrix(op)
            A = speye(prod(n));
        end
        
        function op = compose(idop, other_op)
            % COMPOSE Create the composition of two operators.
            op = other_op;
        end
        
        function op = inv(idop)
            % INV Invert the operator.
            op = idop;
        end
        
    end
end
