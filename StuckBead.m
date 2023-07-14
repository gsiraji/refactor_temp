classdef StuckBead < Bead
    methods
        function obj = StuckBead(varargin)
            % Parse inputs
            p = inputParser;
            p.CaseSensitive = false;
            addParameter(p, 'position', [0 0], @isnumeric);
            addParameter(p, 'theta', 0, @isnumeric);
            addParameter(p, 'speed', 0, @isnumeric);
            addParameter(p, 'mass', 1e-5, @isnumeric);
            addParameter(p, 'bead_size', 1, @isnumeric);
            addParameter(p, 'inertia', 0, @isnumeric);
            parse(p, varargin{:});

            % Call superclass constructor
            obj = obj@Bead();

            % Assign parsed values to object properties
            fields = fieldnames(p.Results);
            for i = 1:numel(fields)
                obj.(fields{i}) = p.Results.(fields{i});
            end

        end

        
        
        function new_velocity = updateVelocityImpl(obj,varargin)
            new_velocity = 0;
        end

        function new_angular_velocity = updateAngularVelocityImpl(obj, concentration_gradient,dt)
            
            new_angular_velocity = 0;

        end

        function new_orientation = updateOrientationImpl(obj,dt)
        
            new_orientation = [0,0];

        end

    end
end