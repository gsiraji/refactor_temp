classdef Experiment
    %EXPERIMENTSETUP Summary of this class goes here
    %   Detailed explanation goes here
    
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

            if ~isempty(p.Results.('parameter_file'))
                obj.simulation_parameter_set = readParametersCSV(p.Results.('parameter_file'));
            end


        end
        
        function selfcaging(obj)
            N = 32; L = 0.04; dt = 0.01; time_max = 200;h = L/N;
            clockmax = ceil(time_max/dt); bead_count =3;mass = 1e-7;
            inertia = 1e-12;repulsion_coefficient = 25e-5;
            sample_count = 30; beadDistance = 0.25;
%             simulation_parameter_set.depolymerization_constant = 1000;
            
            simulation_parameter_set = struct;
            simulation_parameter_set.model = 'selfcaging_beads';
%             obj.simulation_parameter_set.delta_function_mode = 6;

            simulation_parameter_set.repulsion_coefficient = repulsion_coefficient;
            simulation_parameter_set.time_max = time_max; % Run until time

            logger_parameter_set=struct;
            logger_parameter_set.time_max = time_max; % Run until time
            logger_parameter_set.dt = dt; % Time step instability in speed when dt=>0.04
            logger_parameter_set.clockmax = clockmax;
            logger_parameter_set.sample_count = sample_count;
            logger_parameter_set.N = N;
            logger_parameter_set.bead_count = bead_count;
            movie_parameters.frame_count = sample_count;% number of video frames
            movie_parameters.step4video = ceil(clockmax/sample_count);% Video frame frequency^-1
            movie_parameters.movie = 0;
            movie_parameters.grid = repmat((0:1:(N-1)).*h,[N,1]); % ygrid
            movie_parameters.video_ID = strjoin({obj.name,num2str(obj.trial_number),obj.variable_of_interest},'_');

            simulation_parameter_set.bead_count = bead_count;
            simulation_parameter_set.mass = mass;
            simulation_parameter_set.inertia = inertia;
            simulation_parameter_set.N = N;
            simulation_parameter_set.L = L;
            
            
            
            myMovieMaker=MovieMaker(movie_parameters);
            myMovieMaker = myMovieMaker.makeMovie;
            


            position = zeros(simulation_parameter_set.bead_count,2);
            for i=1:simulation_parameter_set.bead_count
%                 myrand = rand(1,2)-0.5;
                position(i,:) = L/2*ones(1,2);
                
            end 
            position(1,:) = L/2+L*[beadDistance,0];
            switch simulation_parameter_set.bead_count
                case {2,3}
                    position(2,:) = L/2-L*[beadDistance,0];
            end
            theta = 2*pi.*rand([bead_count,1]); %[pi/2+pi/4,pi/4];
            
            for i = 1:simulation_parameter_set.bead_count
                % initalize beads
                simulation_parameter_set.beadList(i) = SelfCagingBead('mass',simulation_parameter_set.mass, 'position'...
                    ,position(i,:),'theta',theta(i),'inertia',inertia);

            end
            sim_parameters = simulation_parameter_set;
%             mySimulation = BeadSimulation(parameterSet,myMovieMaker);
            myLogger = Logger(obj,logger_parameter_set);
%             mySimulation = BeadSimulation('beadList', obj.simulation_parameter_set.beadList,'bead_count',obj.simulation_parameter_set.bead_count);%obj.simulation_parameter_set
            mySimulation = BeadSimulation(sim_parameters);
            mySimulation.run(myLogger,myMovieMaker);

        end
        function inertia(obj)

                       
            obj.simulation_parameter_set.beadList = InertiaBead();
            mySimulation = BeadSimulation(obj.simulation_parameter_set);%
            [position, theta] = generateBeadPosition(mySimulation);
            
            for i = 1:mySimulation.bead_count
                % initalize beads
                mySimulation.beadList(i) = InertiaBead('mass',mySimulation.mass, 'position'...
                    ,position(i,:),'theta',theta(i),'inertia',mySimulation.inertia);

            end

            myMovieMaker=MovieMaker(obj,mySimulation);
            myMovieMaker = myMovieMaker.makeMovie;

            myLogger = Logger(obj,mySimulation);
            mySimulation.run(myLogger,myMovieMaker);

        end

        function gradientBeads(obj)

                       
            obj.simulation_parameter_set.beadList = GradientBead();
            mySimulation = BeadSimulation(obj.simulation_parameter_set);%
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





%             obj.simulation_parameter_set.bead_distance = ...
%                 obj.simulation_parameter_set.bead_distance*0.025;
%             obj.simulation_parameter_set.bead_count = 2;
%             obj.simulation_parameter_set.beadList = StuckBead();
%             mySimulation = BeadSimulation(obj.simulation_parameter_set);%
%             [position, theta] = generateBeadPosition(mySimulation);
%             obj.movie = 0;
%             for i = 1:mySimulation.bead_count
%                 % initalize beads
%                 mySimulation.beadList(i) = StuckBead('mass',mySimulation.mass, 'position'...
%                     ,position(i,:),'theta',theta(i),'inertia',mySimulation.inertia);
% 
%             end
% 
%             myMovieMaker=MovieMaker(obj,mySimulation);
%             myMovieMaker = myMovieMaker.makeMovie;
% 
%             myLogger = Logger(obj,mySimulation);
%             mySimulation.run(myLogger,myMovieMaker);
            
        end

        function original(obj) 
%             obj.name = 'actin_original';
           simulation_parameter_set = readParametersCSV('inertia_sim_parameters.csv');
            simulation_parameter_set.model = 'original';
            
            % initalize beads
            bead_distance = 0.25;
            [position, ~] = generateBeadPosition(simulation_parameter_set,bead_distance);

            for i = 1:simulation_parameter_set.bead_count
                simulation_parameter_set.beadList(i) = OriginalBead('position',...
                    position(i,:));
            end

            myMovieMaker=MovieMaker(obj,simulation_parameter_set);
            myMovieMaker = myMovieMaker.makeMovie;

            myLogger = Logger(obj,simulation_parameter_set);

            mySimulation = BeadSimulation(simulation_parameter_set);%

            mySimulation.run(myLogger,myMovieMaker);



        end

    end
end

                        %         save(strcat('c_f',num2str(simnum),'t',num2str(iii),'.mat'),"vorticity",'-mat')
                        %         save(strcat('c_p',num2str(simnum),'t',num2str(iii),'.mat'),"vorticity2",'-mat')
