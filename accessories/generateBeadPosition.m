function [position, theta] = generateBeadPosition(simulation_parameter_set)




bead_distance = simulation_parameter_set.bead_distance;
position = zeros(simulation_parameter_set.bead_count,2);
for i=1:simulation_parameter_set.bead_count
    %                 myrand = rand(1,2)-0.5;
    position(i,:) = simulation_parameter_set.L/2*ones(1,2);

end
position(1,:) = simulation_parameter_set.L/2+simulation_parameter_set.L*[bead_distance,0];
switch simulation_parameter_set.bead_count
    case {2,3}
        position(2,:) = simulation_parameter_set.L/2-simulation_parameter_set.L*[bead_distance,0];
end
theta = pi/4.*ones([simulation_parameter_set.bead_count,1]); %[pi/2+pi/4,pi/4];