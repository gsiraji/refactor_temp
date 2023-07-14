function s = readParametersCSV(csv_file)

% Open the file for reading
% Open file
    fileID = fopen(csv_file, 'r');
    
    % Read lines into a cell array
    lines = textscan(fileID, '%s', 'Delimiter', ',');
    lines = lines{1};
    
    % Initialize struct
    s = struct();
    
    % Loop over lines and parse field-value pairs
    for i = 3:2:numel(lines)
        % Split line into field and value
%         parts = strsplit(lines{i}, '=');
        field = lines{i};
        value = str2double(lines{i+1});
        
        % Convert value to number if possible
        if isnan(value)
            value = lines{i+1};
        end
     
        % Add field-value pair to struct
        s.(field) = value;
    end


    s.h = s.L/s.N; % Grid spacing

    s.dt = (0.01/(s.L/32)^2)*s.h^2; % Time step instability in speed when dt=>0.04
    s.clockmax = ceil(s.time_max/s.dt);

    % Close file
    fclose(fileID);
end