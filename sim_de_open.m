function y=sim_de_open(c,h,v,q_initial,T_de,time)

global  N us length concentration  beishu beishu_Ds0 beishu_Dk beishu_Kf delta_t t_ambient concentration_condens
N=400;%每个床200节点
us=v;%表观速度,m/s
length=h;%吸附剂填充长度,m
concentration=c;%入口浓度ppm
t_ambient=20;%环境温度



beishu=1.8;
beishu_Dk=1;
beishu_Kf=9;%原来是9
beishu_Ds0=5;%原来是15
delta_t=60;%时间步长 s
% time=43200*2*2;%模拟总时长，单位s
y0=zeros(5*(N+1),1);
    
%床一初始化

y0(1)=concentration;
y0(2:N+1)=0.0000001;
% y0(N+2)=0;
% y0(N+3:2*(N+1))=2;
y0(N+2:2*(N+1))=q_initial;
y0(2*N+3)=T_de+273.15;
y0(2*N+4:3*(N+1))=t_ambient+273.15;
y0(3*(N+1)+1:4*(N+1))=t_ambient+273.15;
y0(4*(N+1)+1:5*(N+1))=t_ambient+273.15;

[t,y]=ode15s(@de_open,(0:delta_t:time),y0); %循环脱附