function w = phi6Y(r)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%phi_6pt:  6pt delta function with K_delta=constant second moment

w=zeros(6,size(r,1),6);
K_delta=(59/60)*(1-sqrt(1-(3220/3481)));

Alpha=28;

beta=(9/4)-(3/2)*(K_delta+r.^2)+((22/3)-7*K_delta)*r-(7/3)*r.^3;
gamma=(1/4)*( ((161/36)-(59/6)*K_delta+5*K_delta^2)*(1/2)*r.^2 + (-(109/24)+5*K_delta)*(1/3)*r.^4 + (5/18)*r.^6 );
discr=beta.^2-4*Alpha*gamma;
flag=0;
if(any(discr<0))
   flag=1
end


pm3=(-beta+sign((3/2)-K_delta)*sqrt(discr))/(2*Alpha);
w(:,:,6)=repmat(pm3',[6,1,1]) ;
w(:,:,5)= repmat((-3*pm3 - (1/16) + (1/8)*(K_delta+r.^2) + (1/12)*(3*K_delta-1)*r + (1/12)*r.^3)',[6,1,1]) ;
w(:,:,4)=  repmat((2*pm3 + (1/4)                   +  (1/6)*(4-3*K_delta)*r -  (1/6)*r.^3)',[6,1,1]) ;
w(:,:,3)=  repmat((2*pm3 + (5/8)  - (1/4)*(K_delta+r.^2))',[6,1,1]) ;
w(:,:,2)= repmat((-3*pm3 + (1/4)                   -  (1/6)*(4-3*K_delta)*r +  (1/6)*r.^3)',[6,1,1]) ;
w(:,:,1)=    repmat((pm3 - (1/16) + (1/8)*(K_delta+r.^2) - (1/12)*(3*K_delta-1)*r - (1/12)*r.^3)',[6,1,1]) ;
end


