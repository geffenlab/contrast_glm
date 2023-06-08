classdef ToyModelNeuron < handle
   
    properties
        simulation
        base_rate = 50 % firing rate at the operating point
        beta = 0.1 % b parameter in the documentation
        % gain control: 0 for no gain control, 1 perfectly compensates
        % contrast change so that average firing rate is kept constant
        % (this holds in a linearized approximation valid when the contrast
        % is not too large), intermediate values interpolate
        gain_control = 1 % xi parameter in the documentation
        operating_point % operating point for gain control
        
        h % linear drive
        l % firing rate
        g % gain
        y % spikes
    end
    
    methods
        
        function obj = ToyModelNeuron(simulation, operating_point_offset, gain_control)
            obj.simulation = simulation;
            obj.operating_point = obj.simulation.x0 + operating_point_offset;
            obj.gain_control = gain_control;
        end
        
        function g = gain(obj, contrast)
            g = obj.gain_control * obj.simulation.mean_contrast./contrast +...
                (1-obj.gain_control) * 1;
        end
        
        function [lambda_scaled, this_gain] = gain_scaling(obj, lambda, contrast)
            
            this_gain = obj.gain(contrast);
            lambda_scaled = log(obj.base_rate) + this_gain .* (lambda-log(obj.base_rate));

        end
        
        function [l, linear_drive, g] = lambda(obj, x, contrast)
            linear_drive = log(obj.base_rate) + obj.beta.*(x-obj.operating_point);
            [drive_with_gain, g] = obj.gain_scaling(linear_drive, contrast);
            l = exp(drive_with_gain);
        end
        
        function y = simulate(obj, x, contrast)
            [obj.l, obj.h, obj.g] = obj.lambda(x, contrast);
            obj.y = poissrnd(obj.l);
            y = obj.y;
        end
    end
    
end