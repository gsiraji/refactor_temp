function c = cross2D(a,b)
%CROSS  Vector cross product.
%   C = CROSS2D(A,B) returns the two-dimensional (2D) cross product...
%   of the 2D vectors A and B.
%   That is, C = a_1*b_2 - b1_*a_2.  A and B must be 2 element
%   
%
%
%   Class support for inputs A,B:
%      float: double, single

% Calculate cross product
c = a(:,1).*b(:,2)-a(:,2).*b(:,1);


