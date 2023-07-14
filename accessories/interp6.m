function U=interp6(bead_simulation,X,u)
%             u = obj.concentration.field;
U=zeros(size(X,1),2);
s=X/bead_simulation.h; % Get body position relative to grid
i=floor(s);
r=s-i;
w=phi6X(r(:,1)).*phi6Y(r(:,2));

w = permute(w, [1,3,2]); %Reogranize, this is quite fast

for k=1:size(X,1)
    i1=mod((i(k,1)-3):(i(k,1)+2),bead_simulation.N)+1; % Find adjacent fluid cells
    i2=mod((i(k,2)-3):(i(k,2)+2),bead_simulation.N)+1;

    ww = w(:,:,k);
    UU = [sum(sum(ww.*u(i1,i2,1))), sum(sum(ww.*u(i1,i2,2)))]; %Interpolate
    U(k,:)= UU;
end
end