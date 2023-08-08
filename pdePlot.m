colors = ['k', 'r', 'c'];
i = 1;
for R = L.*[0.01,0.025,0.05]

L = 0.04;
gamma = 10^-5;
D2 = 1.3e-7;
g0 = 3.3e-8*2*L;
% R = L*0.025;
gL = g0*0.1;


a = g0-gL*gamma*L/D2;
b = g0 - gL*gamma*(L-R/2)/D2;
d = a;
N = 1000;
linewidth = 2.5;


x1 = linspace(-L,-R+R/3,N);
x2 = linspace(-(1.001)*R,(1.001)*R,N);
x3 = linspace(R-R/3,L,N);

g1 = a + (a - g0)*x1/L;
g2 = b.*ones(1,N);
g3 = d + (g0-d)*x3/L;

plot(x1, g1, LineWidth=linewidth ,Color=colors(i));
hold on;
plot(x2,g2, LineWidth=linewidth ,Color=colors(i))
plot(x3,g3, LineWidth=linewidth ,Color=colors(i))
set(gca,fontname='Times New Roman',fontsize = 18)
xlabel("x (cm)")
ylabel("g^f(x)")

i = i+1;
end