classdef (Abstract) Bead
    %BEAD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        position = [0 0]
        orientation = [0 1]
        velocity = [0 0]
        acceleration = [0 0]
        speed = 0
        theta = 0
        angular_velocity = 0
        mass = 1e-5
        inertia = 0;
        repulsion_coefficient = 0;%Interaction/repulsion coefficient
        interaction_range
        force_type = 2;
        bead_size = 0.04/32
    end
    
    methods

        function obj = Bead(varargin)
            %BEAD 
            %   the paticles in the simulation
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
            
            % Assign parsed values to object properties
            fields = fieldnames(p.Results);
            for i = 1:numel(fields)
                obj.(fields{i}) = p.Results.(fields{i});
            end

            obj.orientation = [cos(obj.theta), sin(obj.theta)];
            obj.velocity = obj.speed * obj.orientation;
            obj.acceleration = 0;
            obj.interaction_range = obj.bead_size;
        end


        function obj = updateAngularVelocity(obj,concentration_gradient,dt)
            %UPDATEANGULARVELOCITY Update the bead angular velocity.
            %   This method updates the bead's angular velocity based on the provided input
            %   arguments. 
            %   Name-Value Pair Arguments:
            %
            %
            %   - Concentration Gradient: A scalar or vector that specifies the concentration
            %     gradient at the bead position. 
            %
            %   - dt: A scalar that specifies the time step to use for updating the
            %     bead's velocity. Default is 1e-3.
            %
            %
            %   Returns:
            %   -------
            %   - obj: An instance of the Bead class with the updated angular velocity.
            
            
            % the switch catches divided by 0 cases
            % in a future version, this may be replaced
            % by a try block
            switch obj.inertia 
                case 0
                    obj.angular_velocity = 0;
                otherwise
                    % calclate rhatx(rhatxC)
                    obj.angular_velocity = obj.updateAngularVelocityImpl(concentration_gradient,dt);
            end
        end
        
        function obj = updatePosition(obj, dt)
            
            obj.position = obj.position + obj.velocity*dt ;
        end

        function obj = updateOrientation(obj, dt)
            
            obj.orientation = obj.updateOrientationImpl(dt);

            
           
        end
        

        function obj = updateVelocity(obj,varargin)
            %UPDATEVELOCITY Update the bead velocity.
            %   This method updates the bead's velocity based on the provided input
            %   arguments. You can specify optional name-value pairs to customize the
            %   behavior of the method.
            %
            %   Name-Value Pair Arguments:
            %
            %
            %   - Concentration: A scalar or vector that specifies the concentration
            %     at the bead position. If a vector is provided, the function will
            %     use the value corresponding to the current direction of the bead.
            %
            %   - dt: A scalar that specifies the time step to use for updating the
            %     bead's velocity. Default is 1e-3.
            %
            %   - direction: A vector that specifies the direction of motion of the
            %     bead. Default is obj.orientation.
            %
            %   Returns:
            %   -------
            %   - obj: An instance of the Bead class with the updated velocity.

            % Set default values
            default_dt = 1e-3;
            default_direction = obj.orientation;

            % Parse name-value pairs
            p = inputParser;
            addParameter(p,'concentration',[1 1]);
            addParameter(p,'dt',default_dt);
            addParameter(p,'direction',default_direction);
            parse(p,varargin{:});
            concentration = p.Results.concentration;
            dt = p.Results.dt;
            direction = p.Results.direction;

            % Call implementation function
            obj.velocity = obj.updateVelocityImpl('concentration',concentration,'dt', dt,'direction', direction);
        end



    end

%     methods(Abstract)
%         new_velocity = updateVelocityImpl(obj, concentration, dt, direction)
%     end

end

