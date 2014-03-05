classdef QuasiNewtonOperator < UpdatableOperator
    properties
        B % the Hessian approximation
        H % Approximation of the inverse Hessian
        mode = 'bfgs'
    end

    methods(Static,Access=protected)
        function y=apply_or_solve(A, B, x)
            if ~isempty(A)
                y = A*x;
            else
                y = B\x;
            end
        end
    end
    methods
        function qnop=QuasiNewtonOperator(B0, H0, mode)
            qnop.B = B0;
            qnop.H = H0;
            qnop.mode = mode;
        end
        
        
        function y=apply(qnop, x)
            y = qnop.apply_or_solve(qnop.B, qnop.H, x);
        end
        
        function x=solve(qnop, y)
            x = qnop.apply_or_solve(qnop.H, qnop.B, y);
        end
        
        function s=size_impl(qnop)
            if isempty(qnop.B)
                s = size(qnop.B);
            else
                s = size(qnop.H);
            end
        end
        
        function A=asmatrix(qnop)
            [A, Ainv] = deal(qnop.B, qnop.H);
            if isempty(A)
                A=inv(Ainv);
            end
        end
        
        function qnop = update(qnop, y, s)
            [qnop.B, qnop.H] = qn_matrix_update(qnop.mode, qnop.B, qnop.H, y, s);
        end
    end
    
end
