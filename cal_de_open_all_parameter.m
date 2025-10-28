clear
clc
tic
global  N us length concentration  beishu beishu_Ds0 beishu_Dk beishu_Kf delta_t t_ambient concentration_condens T_de T_condens

c_in=0.000001;%脱附入口浓度
N=400;%每个床200节点
t_ambient=20;%环境温度
time=43200*2*2;%模拟总时长，单位s

h_series=0.3;% 0.1:0.1:1;
v_series=0.1;% 0.05:0.05:0.5;
q_initial_series=2.305944889;%[1.583422901	1.903413968	2.15535571	2.305944889	2.355732823];
c_ad_series=0.001;%[0.0002	0.0003	0.0005	0.001	0.002];
T_de_series=160;%80:20:220;%80:10:230
T_condens_series=0;%-25:5:25;%-15:15:15;

number1=size(h_series,2);
number2=size(v_series,2);
number3=size(q_initial_series,2);
number4=size(T_de_series,2);
number5=size(T_condens_series,2);

number=number1*number2*number3*number4*number5;
result=zeros(number,16);
position=1;
for i=1:number1

    for j=1:number2

        for k=1:number3

            for m=1:number4

                for n=1:number5
%                     position=number2*number3*number4*number5*(i-1)+number3*number4*number4*number5*(j-1)...
%                         +number4*number5*(k-1)+number5*(m-1)+n;




h=h_series(i);
v=v_series(j);
q_initial=q_initial_series(k);
c_ad=c_ad_series(k);
T_de=T_de_series(m);
T_condens=T_condens_series(n);
% 
% h=0.3;
% v=0.1;
% T_de=150;%脱附温度
% q_initial=2.3;
% c_ad=0.001;%吸附阶段浓度
% T_condens=0;%基准工况的冷凝温度设置为0℃
try


%开始计算开式脱附
y=sim_de_open(c_in,h,v,q_initial,T_de,time);
energy_data_de_open=energy_consumption_de_open(N,y,v,h,T_de,T_condens,c_ad,q_initial);

Q_de_open_unit=energy_data_de_open(6);
T_de_open=energy_data_de_open(7);%脱附结束的时间
q_de_open=energy_data_de_open(8);%脱附量,也是吸附量
q_de_recovery=energy_data_de_open(9);
q_de_open_left=energy_data_de_open(10);

%数据记录
%计算参数
result(position,1)=h;
result(position,2)=v;
result(position,3)=q_initial;
result(position,4)=c_ad;
result(position,5)=T_de;
result(position,6)=T_condens;
%能耗数据
% energy_data(1)=P_fan_de_unit;
% energy_data(2)=Q_heat_de_unit;
% energy_data(3)=Q_cool_sensible_unit;
% energy_data(4)=Q_cool_latent_unit;
% energy_data(5)=Q_cool_de_unit;
% energy_data(6)=Q_de_all_unit;
result(position,7)=energy_data_de_open(1);
result(position,8)=energy_data_de_open(2);
result(position,9)=energy_data_de_open(3);
result(position,10)=energy_data_de_open(4);
result(position,11)=energy_data_de_open(5);
result(position,12)=energy_data_de_open(6);%总单位回收能耗

result(position,13)=energy_data_de_open(7);%脱附时间
result(position,14)=energy_data_de_open(8);
result(position,15)=energy_data_de_open(9);%回收量
result(position,16)=energy_data_de_open(10);
end
position

position=position+1;
                end
            end
        end
    end
end
% xlswrite('计算结果.xlsx',result,'开始脱附首次计算')
toc






















