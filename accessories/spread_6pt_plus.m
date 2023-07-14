function f=spread_6pt_plus(F,X)
global N h k_c;

%% test
% L = 0.08;
% N = 64;
% h = L/N;
% X = [32,64;32,32.5].*h;
% F = -ones(size(X,1),1);
%%
c=1/(h*h);
f=zeros(N,N);
s=X/h; % Get body position relative to grid
i=floor(s);
r=s-i;
w1 = get_w1(phi1_6pt(r(:,1)),k_c);
w2 = get_w2(phi2_6pt(r(:,2)),k_c);
w=w1.*w2;

w = permute(w, [1,3,2]); %Reogranize, this is quite fast

for k=1:size(X,1)
  i1=mod((i(k,1)-3*k_c):(i(k,1)+3*k_c-1),N)+1; % Find adjacent fluid cells
  i2=mod((i(k,2)-3*k_c):(i(k,2)+3*k_c-1),N)+1;

  ww = w(:,:,k);
  f(i1,i2)=f(i1,i2)+(c*F(k))*ww; 
end
