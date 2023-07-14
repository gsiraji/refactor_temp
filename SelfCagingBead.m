classdef SelfCagingBead < Bead
    methods
        function obj = SelfCagingBead(varargin)
           defaults = {'position', [0, 0];...
                'speed', 0;...
                'theta', 0;...
                'mass', 5e-5/(3e4);...
                'bead_size', 0.04/32;...
                'repulsion_coefficient', 5e-5;
                'inertia',(0.04/32)/(7e3)*5e-5;
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
            new_velocity = obj.speed.*obj.orientation ...
                - concentration./obj.mass...
                + obj.acceleration*dt;
        end

        function new_angular_velocity = updateAngularVelocityImpl(obj, concentration_gradient, dt)
            % calclate rhatx(rhatxC)
            cross_product = cross2D(obj.orientation,concentration_gradient);
            cross_product = cross_product.*fliplr(obj.orientation).*[1,-1];

            % calculate the angular velocity
            new_angular_velocity = (1/obj.inertia)*cross_product;
        end

        function new_orientation = updateOrientationImpl(obj,dt)
            
            obj.orientation = obj.orientation + ...
                + obj.angular_velocity*dt ;

            %normalize the orientation
            new_orientation = obj.orientation./vecnorm(obj.orientation');

        end

    end
end

