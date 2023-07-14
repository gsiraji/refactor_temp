function w = get_w2(w,k_c)
%creates a k_c*6xk_c*6 delta function spread 
%   k_c is the effective size of the bead relative to the grid
%   w is the original 6pt delta function
%   direction sets the x/y parameter


wsize = size(w,1);

w = repmat(w, [k_c k_c 1]);

w = reshape(w, [k_c*wsize 1 k_c*wsize]);


end

