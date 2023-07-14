function w = phi6X(r)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%phi_6pt:  6pt delta function with K=constant second moment

w=zeros(6,size(r,1),6);
K=(59/60)*(1-sqrt(1-(3220/3481)));

alpha=28;

beta=(9/4)-(3/2)*(K+r.^2)+((22/3)-7*K)*r-(7/3)*r.^3;
gamma=(1/4)*( ((161/36)-(59/6)*K+5*K^2)*(1/2)*r.^2 + (-(109/24)+5*K)*(1/3)*r.^4 + (5/18)*r.^6 );
discr=beta.^2-4*alpha*gamma;
flag=0;
if(any(discr<0))
   flag=1
end

pm3=(-beta+sign((3/2)-K)*sqrt(discr))/(2*alpha);
w(6,:,:)=repmat(pm3,[1,1,6]) ;
w(5,:,:)= repmat(-3*pm3 - (1/16) + (1/8)*(K+r.^2) + (1/12)*(3*K-1)*r + (1/12)*r.^3,[1,1,6]) ;
w(4,:,:)=  repmat(2*pm3 + (1/4)                   +  (1/6)*(4-3*K)*r -  (1/6)*r.^3,[1,1,6]) ;
w(3,:,:)=  repmat(2*pm3 + (5/8)  - (1/4)*(K+r.^2),[1,1,6]) ;
w(2,:,:)= repmat(-3*pm3 + (1/4)                   -  (1/6)*(4-3*K)*r +  (1/6)*r.^3,[1,1,6]) ;
w(1,:,:)=    repmat(pm3 - (1/16) + (1/8)*(K+r.^2) - (1/12)*(3*K-1)*r - (1/12)*r.^3,[1,1,6]) ;
end

