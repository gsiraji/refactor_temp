classdef Experiment
    %EXPERIMENT
    %   Run one of the following experiments:
    % stuckBeads, inertia, gradientBeads, original
    
    properties
        name
        trial_number
        variable_of_interest
        movie = 0
        simulation_parameter_set
        parameter_file = ''
    end
    
    methods
        function obj = Experiment(varargin)
            
            defaults = {'name', 'actin_beads';
                'trial_number', 1;
                'variable_of_interest', '';
                'movie',0;
                'simulation_parameter_set',struct();
                'parameter_file',''
                };

            % Parse inputs
            p = inputParser;
            p.CaseSensitive = false;
            for i = 1:1:length(defaults)
                addParameter(p, defaults{i,1}, defaults{i,2})
            end
            parse(p, varargin{:})
            
            % Assign parsed values to object properties
            fields = fieldnames(p.Results);

            for i = 1:numel(fields)
                    obj.(fields{i}) = p.Results.(fields{i});
            end
            
            % check if there is a parameter_file to read
            if ~isempty(p.Results.('parameter_file'))
                obj.simulation_parameter_set = readParametersCSV(p.Results.('parameter_file'));
            end


        end
        

        function inertia(obj)

            % initiate the beads           
            obj.simulation_parameter_set.beadList = InertiaBead();

            % create the simulation object using the simulation parameter
            % set
            mySimulation = BeadSimulation(obj.simulation_parameter_set);%
            [position, theta] = generateBeadPosition(mySimulation);
            
            % % The simulation beads
            for i = 1:mySimulation.bead_count
                % initalize beads
                mySimulation.beadList(i) = InertiaBead('mass',mySimulation.mass, 'position'...
                    ,position(i,:),'theta',theta(i),'inertia',mySimulation.inertia);

            end

            % initialize the movieMaker
            myMovieMaker=MovieMaker(obj,mySimulation);
            myMovieMaker = myMovieMaker.makeMovie;

            %initialize the logger
            myLogger = Logger(obj,mySimulation);
            mySimulation.run(myLogger,myMovieMaker);

        end


        function gradientBeads(obj)

                       
            obj.simulation_parameter_set.beadList = GradientBead();
            mySimulation = BeadSimulation(obj.simulation_parameter_set);
            [position, theta] = generateBeadPosition(mySimulation);
            
            for i = 1:mySimulation.bead_count
                % initalize beads
                mySimulation.beadList(i) = GradientBead('mass',mySimulation.mass, 'position'...
                    ,position(i,:),'theta',theta(i),'inertia',mySimulation.inertia);

            end

            myMovieMaker=MovieMaker(obj,mySimulation);
            myMovieMaker = myMovieMaker.makeMovie;

            myLogger = Logger(obj,mySimulation);
            mySimulation.run(myLogger,myMovieMaker);

        end

        function stuckBeads(obj)
            
            obj.simulation_parameter_set.beadList = StuckBead();
            mySimulation = BeadSimulation(obj.simulation_parameter_set);%
            [position, theta] = generateBeadPosition(mySimulation);
            
            for i = 1:mySimulation.bead_count
                % initalize beads
                mySimulation.beadList(i) = StuckBead('mass',mySimulation.mass, 'position'...
                    ,position(i,:),'theta',theta(i),'inertia',mySimulation.inertia);

            end

            myMovieMaker=MovieMaker(obj,mySimulation);
            myMovieMaker = myMovieMaker.makeMovie;

            myLogger = Logger(obj,mySimulation);
            mySimulation.run(myLogger,myMovieMaker);
            
        end

        function original(obj) 

             obj.name = 'actin_original';

        


        end

    end
end