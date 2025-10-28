function deltaP=pressure_drop_seal(v,L)
dc=0.0045;%颗粒直径
AR=2;
epsilon=0.38;%空隙率
miu=1.86/100000;
rou=1.165;

Vp=3.14/4*dc^2*dc*AR;%颗粒体积
Sp=2*3.14/4*dc^2+3.14*dc^2*AR;%颗粒表面积
phi=4*3.14*(3*Vp/4/3.14)^(2/3)/Sp;%球形度

k11=31.563*(1+0.419*AR)/dc^2.291/AR;
k12=(1-epsilon)^2/epsilon^3;%拆解粘性项系数
k21=3.926*phi^2.804/dc^1.174;
k22=(1-epsilon)*(0.587/phi^1.087-epsilon)/epsilon^3;%拆解惯性项系数
deltaP1=k11*k12*miu*v;%
deltaP2=k21*k22*rou*v^2;
deltaP=(deltaP1+deltaP2)*L;
