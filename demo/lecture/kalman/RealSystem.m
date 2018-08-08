classdef RealSystem < handle
    properties
        x = 0 
        y = 0 
        u = 1
        v = 0 
    end
    
    properties
        t = 0 
        h = 0.01
        sigma_q = 0.02 % process noise
        sigma_r = 0.02/10 % observation noise
    end
    
    methods
        function sys = RealSystem()
        end
        
        function y = get_observation(sys)
            y = [sys.x; sys.y] + multi_normal_sample(1, 0, sys.get_observation_noise_matrix());
        end
        
        function R = get_observation_noise_matrix(sys)
            R = eye( 2 ) * sys.sigma_r^2;
        end
        
        function Q = get_process_noise_matrix(sys)
            Q = eye( 4 ) * (sys.h * sys.sigma_q)^2;
        end
        
        function H = get_observation_matrix(sys)
            H = [1 0 0 0; 0 1 0 0];
        end

        function F = get_transition_matrix(sys)
            hh = sys.h;
            F = [1 0 hh 0; 0 1 0 hh; 0 0 1 0; 0 0 0 1];
        end
        
        function do_step(sys)
            sys.t = sys.t + sys.h;
            F = sys.get_transition_matrix();
            xn = F * sys.get_state();
            xn = xn + multi_normal_sample(1, 0, sys.get_process_noise_matrix());
            xn(3) = xn(3) + sys.h * sin(10*sys.t);
            xn(4) = xn(4) + sys.h * cos(10*sys.t);
            sys.set_state(xn);
        end
        
        function set_state(sys, x)
            sys.x = x(1);
            sys.y = x(2);
            sys.u = x(3);
            sys.v = x(4);
        end
        
        function x=get_state(sys)
            x = [sys.x;sys.y;sys.u;sys.v];
        end

    end
    
end