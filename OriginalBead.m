classdef OriginalBead < Bead
    methods
        function obj = OriginalBead(varargin)
            
            defaults = {'position', [0, 0];...
                'speed', 0;...
                'theta', 0;...
                'mass', 1e-5;...
                'bead_size', 0.04/32;...
                'repulsion_coefficient', 5e-5;
                'inertia',1e-7;
                };

            % Parse inputs
            p = inputParser;
            p.CaseSensitive = false;
            for i = 1:1:size(defaults,1)
                addParameter(p, defaults{i,1}, defaults{i,2})
            end
            parse(p, varargin{:})

            % Call superclass constructor
            obj = obj@Bead(p.Results);

        end
        
        function  new_velocity = updateVelocityImpl(obj,varargin)
            
            input_set = {'concentration','dt','direction'};

            % Parse inputs
            p = inputParser;
            p.CaseSensitive = false;
            addParameter(p, 'concentration', [0 0], @isnumeric);
            addParameter(p, 'dt', 0, @isnumeric);
            addParameter(p, 'direction', [0 1], @isnumeric);
            parse(p, varargin{:});

            concentration = p.Results.('concentration');
            dt = p.Results.('dt');
            direction = p.Results.('direction');


            new_velocity = direction.*concentration./obj.mass + obj.acceleration*dt;
%             obj.velocity = new_velocity;
        end

        function new_angular_velocity = updateAngVelocityImpl(varargin)

            new_angular_velocity = 0;

        end

        function new_orientation = updateOrientationImpl(obj,dt)
        
            new_orientation = [0,0];

        end
        
    end
end