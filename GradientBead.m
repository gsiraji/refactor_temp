classdef GradientBead < Bead
    methods
        function obj = GradientBead(varargin)
            
            defaults = {'position', [0, 0];...
                'speed', 0;...
                'theta', 0;...
                'mass', 1e-5;...
                'bead_size', 0.04/32;...
                'repulsion_coefficient', 5e-5;
                'inertia',1e-7*0.04/32;
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
        
        function new_velocity = updateVelocityImpl(obj,varargin)
            
            % Parse inputs
            p = inputParser;
            p.CaseSensitive = false;
            addParameter(p, 'concentration', [0 0], @isnumeric);
            addParameter(p, 'dt', 0, @isnumeric);
            addParameter(p, 'direction', [0 1], @isnumeric);
            parse(p, varargin{:});

            concentration = p.Results.('concentration');
            dt = p.Results.('dt');

            % u_j = alpha*g^f 
            new_velocity = (concentration/obj.mass...
                + dt*obj.acceleration).*obj.orientation;

            
        end

        function new_angular_velocity = updateAngularVelocityImpl(obj, concentration_gradient,dt)


            % calclate (rhatxneblaC) 
            cross_product = -cross2D(obj.orientation,concentration_gradient);
            new_angular_velocity = (obj.bead_size/obj.inertia)*cross_product;

        end

        function new_orientation = updateOrientationImpl(obj,dt)
            
            obj.theta = obj.theta + ...
                + obj.angular_velocity*dt ;

            %normalize the orientation
            new_orientation = [cos(obj.theta),sin(obj.theta)];

        end


    end
end

