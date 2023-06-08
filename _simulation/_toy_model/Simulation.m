classdef Simulation < handle
   
    properties
       contrast_low = 2 % sigma_L in the documentation
       contrast_high = 5 % sigma_H in the documentation
       x0 = 30 % mu parameter in the documentation
       neurons = ToyModelNeuron.empty();
       n_timesteps = 20 % number of timesteps in both low and high contrast conditions
       n_trials = 500
       
       stim
       contrast
       mean_contrast % harmonic mean of low and high contrast
       dm % design matrix for the glm
       glms = {}
        
    end
    
    methods
        function obj = Simulation(gain_control_values, operating_point_offset)
            %SIMULATION simulation of a gain control model, and GLM fits on
            %the data generated from the simulation.
            %
            % |gain_control_values|: array of xi parameter values. The
            % simulation will have one neuron per entry of the array; that
            % neuron will have xi set to the corresponding value. This
            % allows to simulate a bunch of neurons that all have the same
            % input, but different degrees of gain control. Default:
            % [1,0.5,0] (three neurons: one with optimal gain control, one
            % with no gain control, and one halfway).
            %
            % |operating_point_offset|: offset between the operating point
            % of the gain for the neurons in the simulation, and the mean
            % of the stimulus distribution. With the notation used in the
            % documentaion, this is c-mu. Note that this value is shared
            % among all neurons in the simulation. Default: 0 (c=mu).
            %
            % Note that for the GLM fitts to work, glmnet-for-matlab must
            % be on the matlab path.
            %
            %
            % ---Example---
            % % Instantiate a simulation with the default xi settings, but
            % % an offset of 10 between the center of the stimulus
            % % distribution and the operating point of gain control in the
            % % forward models.
            %
            % s = Simulation([], 10);
            %
            % % Run simulations, fit GLMs and plot results
            %
            % s.plot();
            
            if nargin < 1 || isempty(gain_control_values)
                gain_control_values = [1, 0.5, 0];
            end
            if nargin < 2
                operating_point_offset = 0;
            end  
            for n_id=1:length(gain_control_values)
                gain_control = gain_control_values(n_id);
                obj.neurons(n_id) = ToyModelNeuron(obj, operating_point_offset, gain_control);
            end
            obj.mean_contrast = 1/((obj.contrast_high+obj.contrast_low)/(2*obj.contrast_high*obj.contrast_low));
        end
        
        function get_simulation(obj)
            if isempty(obj.neurons(1).y)
                obj.run_simulation();
            end
        end
        
        function get_design_matrix(obj)
            obj.get_simulation();
            if isempty(obj.dm)
                % define DM
                obj.dm = reshape(obj.stim-obj.x0, [], 1);
                new_predictor = (obj.stim-obj.x0).*obj.mean_contrast ./ obj.contrast;
                obj.dm = [obj.dm, reshape(new_predictor, [], 1)];
                obj.dm = [obj.dm, reshape(obj.mean_contrast./obj.contrast, [], 1)];
            end
        end
        
        function get_model_fit(obj)
            obj.get_design_matrix;

            opts = glmnetSet;
            opts.alpha = 0.95;
            if isempty(obj.glms)
                for n_id = 1:length(obj.neurons)
                    neuron = obj.neurons(n_id);
                    obj.glms{n_id} = cvglmnet(obj.dm,...
                        reshape(neuron.y, [], 1), 'poisson', opts,...
                        [], []);
                end
            end
        end
        
        function predictions = get_model_predictions(obj)
            obj.get_model_fit;
            
            predictions = NaN(numel(obj.stim), length(obj.neurons));
            for n_id=1:length(obj.neurons)
                predictions(:, n_id) = cvglmnetPredict(obj.glms{n_id},...
                    obj.dm, 'lambda_1se', 'response');
            end
            predictions = predictions';
            predictions = reshape(predictions, length(obj.neurons), obj.n_trials, 2*obj.n_timesteps);
        end
        
        function w = gain_modulation_index(obj, beta1, beta2, contrast)
            w = 1+beta2/(beta1+beta2).*(obj.mean_contrast./contrast-1);
        end
            
        function plot(obj)
            
            obj.get_model_fit();
            predictions = obj.get_model_predictions();
            
            for n_id=1:length(obj.neurons)
                neuron = obj.neurons(n_id);
                
                coeffs = cvglmnetPredict(obj.glms{n_id}, [], 'lambda_min', 'coefficients');
                beta0 = coeffs(1);
                beta1 = coeffs(2);
                beta2 = coeffs(3);
                beta3 = coeffs(4);
                this_prediction = squeeze(predictions(n_id,:,:));

                
                figure('Position', [250 300 1480 620]);
                subplot(2,4,1)
                imagesc(obj.stim)
                colorbar
                title('Stimulus x')
                xlabel('Time (step)')
                ylabel('Trial')
                
                subplot(2,4,2)
                imagesc(neuron.h)
                colorbar
                title('Linear drive a+b(x-c)')
                xlabel('Time (step)')
                ylabel('Trial')
                
                subplot(2,4,3)
                xr = 0:1:(obj.x0+3*obj.contrast_high);
                hold on
                
                plot(xr, neuron.lambda(xr, obj.contrast_high), 'r')
                scatter(reshape(obj.stim(:,obj.n_timesteps+1:end), 1, []),...
                    reshape(neuron.y(:,obj.n_timesteps+1:end), 1, []),...
                    'MarkerEdgeColor', 'r', 'Marker', '.')                
                temp_dm = [(xr-obj.x0)',...
                    (xr-obj.x0)'*obj.mean_contrast./obj.contrast_high,...
                    ones(length(xr),1)*obj.mean_contrast./obj.contrast_high];
                plot(xr, cvglmnetPredict(obj.glms{n_id}, temp_dm, 'lambda_min', 'response'), 'r:')
                
                plot(xr, neuron.lambda(xr, obj.contrast_low), 'b')
                scatter(reshape(obj.stim(:,1:obj.n_timesteps), 1, []),...
                    reshape(neuron.y(:,1:obj.n_timesteps), 1, []),...
                    'MarkerEdgeColor', 'b', 'Marker', '.')
                temp_dm = [(xr-obj.x0)',...
                    (xr-obj.x0)'*obj.mean_contrast./obj.contrast_low,...
                    ones(length(xr),1)*obj.mean_contrast./obj.contrast_low];
                plot(xr, cvglmnetPredict(obj.glms{n_id}, temp_dm, 'lambda_min', 'response'), 'b:')
                
                xlabel('Stimulus')
                ylabel('Firing rate')
                legend({['High contrast' newline '(forward model)'], 'Data (spikes)',...
                    'Fit', 'Low contrast', 'Data (spikes)', 'Fit'},...
                    'Location', 'northwest')
                title('Nonlinearity exp[a+gb(x-c)]')
                
                subplot(2,4,4)
                imagesc(neuron.l)
                colorbar
                title('Firing rate')
                xlabel('Time (step)')
                ylabel('Trial')
                
                subplot(2,4,5)
                imagesc(neuron.y)
                colorbar
                title('Spike count')
                xlabel('Time (step)')
                ylabel('Trial')
                
                subplot(2,4,6)
                hold on
                plot(mean(neuron.y))
                errorbar(mean(neuron.y), std(neuron.y));% sqrt(var(neuron.y)-var(obj.stim)))
                title('Empirical firing rate')
                xlabel('Time (step)')
                
                subplot(2,4,7)
                
                imagesc(this_prediction)
                title('Firing rate prediction')
                colorbar
                xlabel('Time (step)')
                ylabel('Trial')                
                
                subplot(2,4,8)
                hold on
                w = obj.gain_modulation_index(beta1, beta2, obj.contrast);
                plot(mean(w), 'k-')

                plot(mean(neuron.g), ':',...
                    'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5)
                title('Gain')
                legend({'Est. gain w', 'True gain g'})
                xlabel('Time (step)')
                ylim([0, 2])
                
                fprintf("\nGain control (xi): %.2g. Parameter estimates:\n"+...
                    "Intercept (beta0) %.2g\nBeta1: %.2g\nBeta2: %.2g\nBeta3: %.2g\n",...
                    neuron.gain_control, beta0, beta1, beta2, beta3);
            end

        end
    end
    
    methods (Access=private)
        function run_simulation(obj)
            obj.stim = [...
                normrnd(obj.x0, obj.contrast_low, obj.n_trials, obj.n_timesteps),...
                normrnd(obj.x0, obj.contrast_high, obj.n_trials, obj.n_timesteps)
                ];
            obj.contrast = [...
                repmat(obj.contrast_low, obj.n_trials, obj.n_timesteps),...
                repmat(obj.contrast_high, obj.n_trials, obj.n_timesteps)
                ];
            for n_id = 1:length(obj.neurons)
                neuron = obj.neurons(n_id);
                neuron.simulate(obj.stim, obj.contrast);
            end
        end 
    end
    
    
        
    
end