function delta_C_vec = calculateDeltaC(bead_simulation)
% uses bead position, concentration field, and a
% set distance to calculate difference in concentration
% between left and right of the beads for two stationary beads
%
% bead_simulation is a BeadSimulation class object 
% this script is 

position = zeros(bead_simulation.bead_count,2);

for i = 1:bead_simulation.bead_count; position(i,:) = bead_simulation.beadList(i).position;end

number_of_values_for_distance = 10;

delta_C_vec = zeros(1,number_of_values_for_distance*bead_simulation.bead_count);

measurement_distance = 0.025*bead_simulation.L;

measurement_distance_vector = linspace(measurement_distance/number_of_values_for_distance,measurement_distance,number_of_values_for_distance);

for i = 1:number_of_values_for_distance

X_shift = repmat([measurement_distance_vector(i),0],[bead_simulation.bead_count, 1]);
c_left = interp6forDeltaC(bead_simulation,position-X_shift);
c_right = interp6forDeltaC(bead_simulation,position+X_shift);
indices = bead_simulation.bead_count*(i-1)+1:bead_simulation.bead_count*i;
delta_C_vec(indices) = (c_left(:,1) - c_right(:,1))';

end

end

function U=interp6forDeltaC(bead_simulation,X)
u = bead_simulation.concentration.field(:,:,1);
U=zeros(size(X,1),1);
s=X/bead_simulation.h; % Get body position relative to grid
i=floor(s);
r=s-i;
w=phi6X(r(:,1)).*phi6Y(r(:,2));

w = permute(w, [1,3,2]); %Reogranize, this is quite fast

for k=1:size(X,1)
    i1=mod((i(k,1)-3):(i(k,1)+2),bead_simulation.N)+1; % Find adjacent fluid cells
    i2=mod((i(k,2)-3):(i(k,2)+2),bead_simulation.N)+1;

    ww = w(:,:,k);
    UU = sum(sum(ww.*u(i1,i2,1))); %Interpolate
    U(k,:)= UU;
end
end