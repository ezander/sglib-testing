classdef SpringModel < BaseModel
    
    properties
        T = 10
    end
    
    methods
        function model = SpringModel(varargin)
            options = varargin2options(varargin);
            [T, options] = get_option(options, 'T', model.T);
            check_unsupported_options(options, mfilename);
            
            model.T = T;
        end
        
        function n=response_dim(model)
            n = 2;
        end
        
        function u=compute_response(model, q)
            u0 = q(1,:);
            v0 = q(2,:);
            m  = q(3,:);
            d  = q(4,:);
            k  = q(5,:);
            T = model.T;
            [x,v] = spring_solve(u0, v0, m, d, k, T);
            u = [x; v];
        end
        
        function y = compute_measurements(model, u)
            y = u(1, end);
        end
    end
end
