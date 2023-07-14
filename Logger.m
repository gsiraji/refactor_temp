classdef Logger
    %LOGs data and stuff
    %   Detailed explanation goes here
    
    properties
        simulation_ID = 0
        slice_file_name
        concentration_data_file_name
        simulation_parameters_file_name
        concentration_well_data_file_name
        concentration_field_file_name
        position_data_file_name

        concentration_field
        
        folder_address
        % total G- and F-actin concentrations
        concentration_table 
        % the difference in concentration from left to right of the beads
        well_depth_vector
        % 1d x-slice of the G-actin concentration
        slice_array
        log_step = 1
        sample_count =30 % how often save data/how many video frames
        sampling_step
        bead_position
        bead_velocity

    end
    
    methods
        function obj = Logger(experiment,beadSimulation)

%             obj.simulation_ID = beadSimulation.simulation_ID;
%             obj.log_step = 1;
            
            obj.folder_address = strcat(pwd, '/actin_data/', experiment.name);

            obj.slice_file_name = strcat(obj.folder_address,'/u_slice/u_slice_',num2str(experiment.trial_number),...
                '_',experiment.variable_of_interest,'.mat');
            
            obj.concentration_field_file_name = strcat(obj.folder_address,'/u_slice/concentration_',num2str(experiment.trial_number),...
                '_',experiment.variable_of_interest,'.mat');

            obj.concentration_data_file_name = strcat(obj.folder_address,'/total_concentration_',experiment.variable_of_interest,...
                '_',num2str(experiment.trial_number),'.csv');
            obj.concentration_well_data_file_name = strcat(obj.folder_address,'/delta_C_data_',experiment.variable_of_interest,...
                '_',num2str(experiment.trial_number),'.csv');
            obj.simulation_parameters_file_name = strcat(obj.folder_address,'/simulation_parameters_',...
                experiment.variable_of_interest,'_',num2str(experiment.trial_number),'.csv');
            obj.position_data_file_name = strcat(obj.folder_address,'/position_velocity_',...
                experiment.variable_of_interest,'_',num2str(experiment.trial_number),'.csv');
            
            obj.sample_count = beadSimulation.sample_count;
            obj.sampling_step = ceil(beadSimulation.clockmax/obj.sample_count);
            obj.concentration_table = zeros(obj.sample_count+1,2);
            obj.bead_position = zeros(obj.sample_count,2*beadSimulation.bead_count);
            obj.bead_velocity = zeros(obj.sample_count,2*beadSimulation.bead_count);
            obj.well_depth_vector = zeros(obj.sample_count,beadSimulation.bead_count*10);
            obj.slice_array = zeros(beadSimulation.N,obj.sample_count);
            obj.concentration_field = zeros(beadSimulation.N,beadSimulation.N,obj.sample_count);

            
        end

        function logInfo(obj,beadSimulation)
            
            mkdir(obj.folder_address)
            mkdir(strcat(obj.folder_address, '/u_slice/'))

            parameters = fieldnames(beadSimulation);
            T1 = string(parameters);
            T1 = T1(1:22);
            T2 = ones(1,22);
            for i = 1:22
                T2(1,i) = beadSimulation.(parameters{i});
            end
            info_table = array2table(T2, 'VariableNames',T1');

            writetable(info_table, obj.simulation_parameters_file_name,'WriteVariableNames',1);

        end


        function obj = logData(obj,beadSimulation)

            % calculate and store the total actin concentration 
            obj.concentration_table(obj.log_step,:) = ...
                reshape(sum(sum(beadSimulation.concentration.field)),[1,2])*beadSimulation.h^2;

            % save the 1d slice array
            obj.slice_array(:,obj.log_step) = beadSimulation.concentration.field(:,beadSimulation.N/2,1)';
            CC = calculateDeltaC(beadSimulation);
            obj.well_depth_vector(obj.log_step,:) = CC;

            % save X and U
            for indx = 1:beadSimulation.bead_count
            obj.bead_position(obj.log_step,2*indx-1:2*indx) = beadSimulation.beadList(indx).position;

            obj.bead_velocity(obj.log_step,2*indx-1:2*indx) = beadSimulation.beadList(indx).velocity;
            end

            % save u(:,:,1)
            obj.concentration_field(:,:,obj.log_step) = beadSimulation.concentration.field(:,:,1);
            

        end

        function obj = nextEntry(obj)

            obj.log_step = obj.log_step +1;

        end

        
        function saveData(obj)
            
            slice_array = obj.slice_array;
            concentration_field = obj.concentration_field;
            concentration_table = table(obj.concentration_table);
            delta_C_table = table(obj.well_depth_vector);

            position_velocity_table = table(obj.bead_position, obj.bead_velocity);

            writetable(concentration_table, obj.concentration_data_file_name,'WriteVariableNames',0)

            writetable(delta_C_table, obj.concentration_well_data_file_name,'WriteVariableNames',0)

            writetable(position_velocity_table, obj.position_data_file_name,'WriteVariableNames',1)

            save(obj.slice_file_name, 'slice_array')

            save(obj.concentration_field_file_name, 'concentration_field')

        end

        function progressMessage(obj)

            formatPerc = 'The simulation is %.1f %% done. \n';
            simulationProgress = obj.log_step/obj.sample_count*100;

            fprintf(formatPerc,simulationProgress)
                    


        end
    end
end

