%修正了回收率的误差
function energy_data=energy_consumption_de_open(N,y,us,h,T_de,T_condens,c_ad,q_initial)

energy_data=zeros(1,10);%数据储存
concentration_condens=10^(11.53-2263/(273.15+T_condens))*10^-6;
[y_max,T1]=max(y(:,N+1));%找到峰值时间
[~,T2]=min(abs(y(T1:end,N+1)'-0.5*c_ad));%结束时间。开式脱附结束选择为吸附浓度的50%
T2=T2+T1-1;
concentrated_index=max(y(:,N+1))/c_ad;%浓缩倍数
q_de=q_initial-mean(y(T2,N+1+1:2*(N+1)));%脱附量
ratio_de=q_de/q_initial;%脱附率
depletion_curve=y(1:T2,N+1);
if y_max <= concentration_condens
    integration1=0;
else
        [~,position1]=min(abs(depletion_curve(1:T1)-concentration_condens));
        [~,position2]=min(abs(depletion_curve(T1:end)-concentration_condens));
        position2=position2+T1-1;
        integration1=trapz(position1:position2,depletion_curve(position1:position2)-concentration_condens);%回收量面积积分
end
        integration2=trapz(1:T2,depletion_curve);%脱附量面积积分
        q_recovery=integration1/integration2*q_de;
%         ratio1=integration1/integration2;%回收量占脱附的比率。就叫冷凝回收率吧
%         ratio2=ratio1*op_Statistics(j,5);%回收量占吸附量的比率，就叫脱附回收率吧

%向量数据记录
q_distribution1=y(T2,N+1+1:2*(N+1));%吸附量分布行向量
t_distribution1=y(T2,2*(N+1)+1:3*(N+1));%温度分布行向量
depleth_curve1=y(:,N+1);%耗竭曲线 列向量
t_out_curve1=y(:,3*(N+1));%出口温度曲线 列向量
q_de_curve1=q_initial-mean(y(:,N+1+1:2*(N+1)),2);%脱附量随时间变化 列向量

% %写数据
% position=(i-1)*number2+j;
% result_data1(position,1)=T1;
% result_data1(position,2)=T2;
% result_data1(position,3)=concentrated_index;
% result_data1(position,4)=q_de;
% result_data1(position,5)=ratio_de;

% q_distribution_data1(:,position)=q_distribution1';%末态吸附量分布记录
% t_distribution_data1(:,position)=t_distribution1';%末态温度分布记录
% depleth_curve_data1(:,position)=depleth_curve1;%脱附曲线记录
% q_de_curve_data1(:,position)=q_de_curve1;%脱附量记录
% non_condens_data(:,position)=min(depleth_curve1,concentration_condens);%不凝气浓度曲线记录
% net_c_de_data(:,position)=depleth_curve1-non_condens_data(:,position);%净脱附速率
% t_out_curve_data1(:,position)=t_out_curve1;%出口温度记录
delta_t=60;
%开始算能耗
%还是计算单位流通面积的能耗，吸附量啥的
density_N2=101325*28/(8.314*(T_de+273.15))/1000;
% Cpg_N2=(28.901-1.5713E-3*(T_de+273.15)+8.0805E-6*(T_de+273.15)^2-28.7256E-9*(T_de+273.15)^3)/28;%kJ/(kg*K)
Cpg_N2_1=4.63e-07*(T_de+273.15)^3+0.0001744*(T_de+273.15)^2+0.01004*(T_de+273.15)+1039;
Cpg_N2_2=4.63e-07*(T_condens+273.15)^3+0.0001744*(T_condens+273.15)^2+0.01004*(T_condens+273.15)+1039;
Cpg_N2=(Cpg_N2_1+Cpg_N2_2)/2/1000;%kJ/(kg*K) 取个过程的平均比热容
flowrate_m=us*1*density_N2;%氮气质流量 kg/s
m_AC=h*1*534;%活性炭质量
DH=48-6*(T_condens+50)/100;%纯组分汽化焓 kJ/mol
% X=(T_condens+273.15-10)/(308.15-T_condens-273.15);
% COP=4.1*10^-6*X^3-0.0025*X^2+0.7057*X-0.4691;
COP=((T_condens+273.15)-130.1)/(320.9-(T_condens+273.15));
deltaP=pressure_drop_seal(us,h);
%下面是闭式过程的能耗
P_fan_de=deltaP*flowrate_m/density_N2/0.5*T2*delta_t/1000;%风机能耗kJ.风机效率50%
Q_heat_de=Cpg_N2*(flowrate_m*T2*delta_t)*(T_de-T_condens)/1;%热负荷kJ,加热效率是1
Q_cool_sensible=Cpg_N2*(flowrate_m*T2*delta_t)*(mean(t_out_curve1(1:T2))-T_condens-273.15)/COP;%冷却显热kJ
Q_cool_latent=q_recovery*m_AC*DH/COP;%冷却潜热kJ
Q_cool_de=Q_cool_sensible+Q_cool_latent;%冷负荷
Q_de_all=P_fan_de+Q_cool_de+Q_heat_de;

P_fan_de_unit=P_fan_de/(q_recovery*m_AC);%转成单位的能耗，就是除以mol
Q_heat_de_unit=Q_heat_de/(q_recovery*m_AC);
Q_cool_sensible_unit=Q_cool_sensible/(q_recovery*m_AC);
Q_cool_latent_unit=Q_cool_latent/(q_recovery*m_AC);
Q_cool_de_unit=Q_cool_de/(q_recovery*m_AC);
Q_de_all_unit=Q_de_all/(q_recovery*m_AC);


energy_data(1)=P_fan_de_unit;
energy_data(2)=Q_heat_de_unit;
energy_data(3)=Q_cool_sensible_unit;
energy_data(4)=Q_cool_latent_unit;
energy_data(5)=Q_cool_de_unit;
energy_data(6)=Q_de_all_unit;
energy_data(7)=T2;
energy_data(8)=q_de;
energy_data(9)=q_recovery;
energy_data(10)=q_initial-q_de;























