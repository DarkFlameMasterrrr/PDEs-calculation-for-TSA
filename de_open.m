function dy=de_open(t,y)

global  N us length beishu beishu_Ds0 beishu_Dk beishu_Kf delta_t t_ambient T_de

P=101325;%气体压力，Pa
R=8.314;%气体常数，J/mol/K

db=0.215;%吸附床直径，m
thic=0.005;%吸附床厚度，m
e=0.34;%床孔隙率

u_set=us/e;%间隙速度,m/s T_de下的设定流速
density_wall=7800;%吸附床墙体密度,kg/m^3
Cpw=500;%墙体热容，J/(kg*K)原来是340
din=0.05;%绝热层厚度。
hs=88;%固体与气体之间的传热系数，W/(m^2*K) 100
hw=60;%气体与壁面之间的传热系数，W/(m^2*K) 60
U=1;%壁面与大气之间的传热系数，W/(m^2*K) 1.6
% U=0;%壁面与大气之间的传热系数，W/(m^2*K)
ep=0.38;%吸附剂孔隙率%压汞实验
dp=0.0045;%颗粒物直径,m
% rpore=2E-8;%颗粒物孔半径，m
rpore=6.8E-9;%颗粒物孔半径，m %这个值是不是需要改一下呀，正常来说Dk算的是哪里的扩散呀，应该不是微孔中的吧
density_particle=534;%颗粒物密度，kg/m^3堆积密度
deltH=-50000;%吸附热，J/mol
Cps=820;%吸附剂热容，J/(kg*K) 厂家给的
as=5/dp;%颗粒物外表面积与体积的比值，1/m，AR=2
aw=db/(thic*(thic+db));%内表面积与墙体体积的比值，1/m
% aa=din/(thic*(thic+db));%绝热平均对数面积与墙体体积的比值，1/m
aa=2*(din+thic)/log((2*(din+thic)+db)/db)/(thic*(thic+db));%绝热平均对数面积与墙体体积的比值，1/m Analysis of purge gas temperature in cyclic TSA process

%等温线参数
%脱附等温线公式采用sips q=M0*exp(M1/T)*B0*exp(B1/T)*y^n/(1+B0*exp(B1/T)*y^n)
M0=2.314;%根据等温线sips拟合得到的系数
M1=44.94;
B0=7.36E-05;
B1=2754;
n=0.5339;
% M0=2.889;%这一组数据是yun论文中的数据
% M1=79.11;
% B0=5.429E-07;
% B1=5292;
% n=1;
y_min_value=0.0001;


Ma=106.165;%邻二甲苯的分子质量，g/mol；
Mg=29;%空气的分子质量，g/mol
%Mg=28;%氮气的分子质量，g/mol
% Va=92.8;%乙酸乙酯的分子体积cm^3/mol
% Vg=20.1;%空气的分子体积cm^3/mol
Va=10;%乙酸乙酯的分子体积cm^3/mol
Vg=29.9;%空气的分子体积cm^3/mol
%Vg=17.9;%空气的分子体积cm^3/mol
z=length/N;%间隔长度,cm
% Ta=25.6+273.15;%温度，K
Ta=t_ambient+273.15;%环境温度，K

%%%%%%%%%%%%边界条件%%%%%%%%%%%%%%
% y(2*N+3)=157+273.15;
% y(2*N+3)=181.9*0.05175*(t/60)^1.227/(1+0.05175*(t/60)^1.227)+273.15;%自定义拟合，langmuir公式
% y(2*N+3)=163.4*exp(0.0009084*t/60) -159.6*exp(-0.07852*t/60)+273.15;%指数拟合
% y(2*N+3)=119.2*exp(-4.894e-05*t/60) -106.3*exp(-0.01329*t/60)+273.15;
% y(2*N+3)=124.4*exp(-0.0001468*t/60) -111*exp(-0.03986*t/60)+273.15;
%  y(2*N+3)=(0.0007046*(t/60)^3 - 0.3569*(t/60)^2 + 234.8*(t/60) + 207) / ((t/60) + 20.35)+273.15;
% y(2*N+3)=162*exp(4.103e-05*t/60) -156*exp(-0.09351*t/60)+273.15;%指数拟合

% y(2*N+3)=73.15+100*t/60/300+273.15;%设置入口温度,用来检测变边界条件是否可行，验证通过
%y0=0.00082;%初始浓度（摩尔分数，3500mg/m3）
% y(1)=concentration;
% if t<56.1*60
%     y(2*N+3)=-0.005715*(t/60)^2+0.623535*t/60+291.96;
% elseif  (56.1*60<t)&&(t<650.1*60)
%         y(2*N+3)=-8.01e-11*(t/60)^4+1.23435e-07*(t/60)^3-6.7569e-05*(t/60)^2+0.0139*t/60+308.485398624855;
%     else
%         y(2*N+3)=35.4+273.15;
% end

dy=zeros(5*(N+1),1);

%定义数据存储数组，前5*(N+1)项为床一的参数。1:N+1为VOC浓度，N+2:2*(N+1)为VOC吸附量，2N+3:3*(N+1)为气体温度，3N+4:4*(N+1)为固体体温度,4N+5:5*(N+1)墙体
%后5*(N+1)项为床二的参数。5N+6:6N+6为VOC浓度，6N+7:7N+7为VOC吸附量，7N+8:8N+8为气体温度，8N+9:9N+9为固体体温度,9N+1010N+10墙体
for ii=1:N+1
    if ii==1
        %床一
        u=u_set*y(ii+2*(N+1))/(T_de+273.15);
        niandu=(0.0335*(y(ii+2*(N+1))-273.15)+17.704)/1000000;%动力粘度是关于温度的函数,Pa*s
        %niandu=(1E-7*(y(ii+2*(N+1))+0.0002*(y(ii+2*(N+1))*(y(ii+2*(N+1))-0.0991*(y(ii+2*(N+1))+20.447)/1000000;%氮气动力粘度
        density_gas=P*Mg/(R*y(ii+2*(N+1)))/1000;%气体密度与温度有关,kg/m^3
        %         Dm=1.01E-3*y(ii+2*(N+1))^1.75*sqrt((Ma+Mg)/(Ma*Mg))/(Va^(1/3)+Vg^(1/3))^2;%分子扩散系数，m^2/s
        % Dm=1.01E-7*y(ii+2*(N+1))^1.75*sqrt((Ma+Mg)/(Ma*Mg))/(Va^(1/3)+Vg^(1/3))^2;%分子扩散系数，m^2/s
        Dm=435.7E-4*y(ii+2*(N+1))^1.5/P*sqrt((Ma+Mg)/(Ma*Mg))/(Va^(1/3)+Vg^(1/3))^2;%分子扩散系数，m^2/s 百度扩散系数，吉利兰Golloland半经验公式
        
        Re=density_gas*u*dp/niandu;%雷诺数计算
        Sc=niandu/density_gas/Dm;%斯密特数计算
        
        %         Kf=2.1/e*Re^0.45*Sc^0.33*Dm/dp;%膜传质系数计算，m/s
        Kf=beishu_Kf*(2.0+1.1*Re^0.6*Sc^(1/3))*Dm/dp;%膜传质系数计算，m/s wakao 雷诺数范围3-3000
        %         Kf=0.357/e*Re^0.64*Sc^0.33*Dm/dp*0.3048/60;%膜传质系数计算，m/s 博士论文
        %         Kf=0.357/e*Re^0.64*Sc^0.33*Dm/dp;%膜传质系数计算，m/s 改数据
        %Sh=2.0+1.8*Re^(1/2)*Sc^(1/3);%舍伍德数计算
        Dk=beishu_Dk*97*rpore*sqrt(y(ii+2*(N+1))/Ma);%克努森扩散系数，m^2/s
        Ds=3.871E-7*beishu_Ds0*exp(deltH*0.45/(R*y(ii+2*(N+1))));%表面扩散系数 文献中是1.05664E-6
        
        %         qm=2.889*exp(79.11/y(ii+3*(N+1)));%langmuir常数，温度T下的最大吸附量，mol/kg
        %         b=5.429E-7*exp(5292/y(ii+3*(N+1)));%langmuir常数，1/Pa
        qm=M0*exp(M1/y(ii+3*(N+1)));
        b=B0*exp(B1/y(ii+3*(N+1)));
        qe=qm*(b*(P*y(ii))^n)/(1+(b*(P*y(ii))^n));%温度T，浓度c下的平衡吸附量，mol/kg
        %         Cdiffq=qm/((qm-y(ii+N+1))^2*b*P);%平衡浓度对吸附量的导数，用于求解有效扩散系数
        Cdiffq=1/n*(y(ii+N+1)/((qm-y(ii+N+1))*b*P^n))^(1/n-1)*(y(ii+N+1)/((qm-y(ii+N+1))^2*b*P^n)+1/((qm-y(ii+N+1))*b*P^n));%yun的论文的基础上根据sips模型推导
        De=Ds+Dk*ep*e*density_gas/density_particle*Cdiffq;%有效扩散系数，m^2/s
        %         k_daoshu=1*density_particle/Kf/density_gas/as*qe/y(ii)+dp/10/De/as;
        y_xiuzheng=max(y(ii),y_min_value);
        k_daoshu=1*density_particle/Kf/density_gas/as*qe/y_xiuzheng+dp/10/De/as;
%         k_daoshu=1*density_particle/Kf/density_gas/as*qe/y(ii)+dp/10/De/as;
        
        k_daoshu=dp/10/De/as;
        % k_daoshu=dp/5/De/as;
        k=1/k_daoshu*beishu;
        %         t
        %         k_record(ii,round(t/delta_t))=k;
        %k=6*2*density_gas*Kf*y0*60*De/(6*density_gas*Kf*y0*dp^2+60*De*density_particle*dp*qe);%总质量传递系数，1/s
%                 k=0.00031;
        %qe_div_c=qm*b*P/(1+b*P*y(ii));%平衡吸附量除以浓度。
        %k_daoshu=density_particle/(Kf*density_gas*as)*qe_div_c+dp/(10*De*as);
        %k=1/k_daoshu;
        
        
        dy(ii)=0;%入口（第一个节点）浓度恒定，故导数为0
        dy(ii+N+1)=k*(qe-y(ii+N+1));%吸附速率方程
        dy(ii+2*(N+1))=0;%入口（第一个节点）温度恒定，故导数为0
        dy(ii+3*(N+1))=(1-e)*hs*as/(density_particle*Cps)*(y(ii+2*(N+1))-y(ii+3*(N+1)))...
            -1/Cps*deltH*k*(qe-y(ii+N+1));%吸附剂能量方程
        dy(ii+4*(N+1))=hw*aw/(density_wall*Cpw)*(y(ii+2*(N+1))-y(ii+4*(N+1)))...
            -U*aa/(density_wall*Cpw)*(y(ii+4*(N+1))-Ta);%吸附床体能量方程
        
        
    elseif ii==N+1
        %床一
         u=u_set*y(ii+2*(N+1))/(T_de+273.15);
        niandu=(0.0335*(y(ii+2*(N+1))-273.15)+17.704)/1000000;%动力粘度是关于温度的函数,Pa*s
        density_gas=P*Mg/(R*y(ii+2*(N+1)))/1000;%气体密度与温度有关,kg/m^3
        %         Dm=1.01E-3*y(ii+2*(N+1))^1.75*sqrt((Ma+Mg)/(Ma*Mg))/(Va^(1/3)+Vg^(1/3))^2;%分子扩散系数，m^2/s
        % Dm=1.01E-7*y(ii+2*(N+1))^1.75*sqrt((Ma+Mg)/(Ma*Mg))/(Va^(1/3)+Vg^(1/3))^2;%分子扩散系数，m^2/s
        Dm=435.7E-4*y(ii+2*(N+1))^1.5/P*sqrt((Ma+Mg)/(Ma*Mg))/(Va^(1/3)+Vg^(1/3))^2;%分子扩散系数，m^2/s 百度扩散系数，吉利兰Golloland半经验公式
        Re=density_gas*u*dp/niandu;%雷诺数计算
        Sc=niandu/density_gas/Dm;%斯密特数计算
        
        %         Kf=2.1/e*Re^0.45*Sc^0.33*Dm/dp;%膜传质系数计算，m/s
        Kf=beishu_Kf*(2.0+1.1*Re^0.6*Sc^(1/3))*Dm/dp;%膜传质系数计算，m/s wakao
        %        Kf=0.357/e*Re^0.64*Sc^0.33*Dm/dp*0.3048/60;%膜传质系数计算，m/s 博士论文
        %         Kf=0.357/e*Re^0.64*Sc^0.33*Dm/dp;%膜传质系数计算，m/s 改数据
        %Sh=2.0+1.8*Re^(1/2)*Sc^(1/3);%舍伍德数计算
        Dk=beishu_Dk*97*rpore*sqrt(y(ii+2*(N+1))/Ma);%克努森扩散系数，m^2/s
        Ds=3.871E-7*beishu_Ds0*exp(0.45*deltH/(R*y(ii+2*(N+1))));%表面扩散系数
        
        qm=M0*exp(M1/y(ii+3*(N+1)));
        b=B0*exp(B1/y(ii+3*(N+1)));
        qe=qm*(b*(P*y(ii))^n)/(1+(b*(P*y(ii))^n));%温度T，浓度c下的平衡吸附量，mol/kg
        %Cdiffq=qm/((qm-y(ii+N+1))^2*b*P);%平衡浓度对吸附量的导数，用于求解有效扩散系数
        Cdiffq=1/n*(y(ii+N+1)/((qm-y(ii+N+1))*b*P^n))^(1/n-1)*(y(ii+N+1)/((qm-y(ii+N+1))^2*b*P^n)+1/((qm-y(ii+N+1))*b*P^n));%yun的论文的基础上根据sips模型推导
        
        De=Ds+Dk*ep*e*density_gas/density_particle*Cdiffq;%有效扩散系数，m^2/s
        y_xiuzheng=max(y(ii),y_min_value);
        k_daoshu=1*density_particle/Kf/density_gas/as*qe/y_xiuzheng+dp/10/De/as;
        k_daoshu=1*density_particle/Kf/density_gas/as*qe/y(ii)+dp/10/De/as;
        k_daoshu=dp/10/De/as;
        % k_daoshu=dp/5/De/as;
        k=1/k_daoshu*beishu;
        %    k_record(ii,round(t/delta_t))=k;
        %k=6*2*density_gas*Kf*y0*60*De/(6*density_gas*Kf*y0*dp^2+60*De*density_particle*dp*qe);%总质量传递系数，1/s
%                 k=0.00031;
        %qe_div_c=qm*b*P/(1+b*P*y(ii));%平衡吸附量除以浓度。
        % k_daoshu=density_particle/(Kf*density_gas*as)*qe_div_c+dp/(10*De*as);
        %k=1/k_daoshu;
        
        dy(ii)=dy(ii-1);%为使方程封闭，假定出口处浓度变化率等于上一级节点
        dy(ii+N+1)=k*(qe-y(ii+N+1));%吸附速率方程
        dy(ii+2*(N+1))=dy(ii+2*(N+1)-1);%为使方程封闭，假定出口处气体温度变化率等于上一级节点
        dy(ii+3*(N+1))=(1-e)*hs*as/(density_particle*Cps)*(y(ii+2*(N+1))-y(ii+3*(N+1)))-deltH/Cps*k*(qe-y(ii+N+1));%吸附剂能量方程
        dy(ii+4*(N+1))=hw*aw/(density_wall*Cpw)*(y(ii+2*(N+1))-y(ii+4*(N+1)))...
            -U*aa/(density_wall*Cpw)*(y(ii+4*(N+1))-Ta);%吸附床体能量方程
        
        
    else
        %床一
         u=u_set*y(ii+2*(N+1))/(T_de+273.15);
        niandu=(0.0335*(y(ii+2*(N+1))-273.15)+17.704)/1000000;%空气动力粘度是关于温度的函数,Pa*s
        % niandu=(1E-7*(y(ii+2*(N+1))+0.0002*(y(ii+2*(N+1))*(y(ii+2*(N+1))-0.0991*(y(ii+2*(N+1))+20.447)/1000000;%氮气动力粘度
        density_gas=P*Mg/(R*y(ii+2*(N+1)))/1000;%气体密度与温度有关,kg/m^3
        %         Dm=1.01E-3*y(ii+2*(N+1))^1.75*sqrt((Ma+Mg)/(Ma*Mg))/(Va^(1/3)+Vg^(1/3))^2;%分子扩散系数，m^2/s
        % Dm=1.01E-7*y(ii+2*(N+1))^1.75*sqrt((Ma+Mg)/(Ma*Mg))/(Va^(1/3)+Vg^(1/3))^2;%分子扩散系数，m^2/s
        Dm=435.7E-4*y(ii+2*(N+1))^1.5/P*sqrt((Ma+Mg)/(Ma*Mg))/(Va^(1/3)+Vg^(1/3))^2;%分子扩散系数，m^2/s 百度扩散系数，吉利兰Golloland半经验公式
        Re=density_gas*u*dp/niandu;%雷诺数计算
        Sc=niandu/density_gas/Dm;%斯密特数计算
        
        % Dz=(20+0.5*Sc*Re)*Dm/e;%有效轴向扩散系数计算，m^2/s
        Dz=0.73*Dm+0.5*u*dp/(1+9.7*Dm/u/dp);%有效轴向扩散系数计算，m^2/s 吴红亮硕士论文
        % Dz=(20+0.5*Re*Sc)*Dm/e;%有效轴向扩散系数计算，m^2/s 闭环再生那篇论文
        %         Kf=2.1/e*Re^0.45*Sc^0.33*Dm/dp;%膜传质系数计算，m/s
        Kf=beishu_Kf*(2.0+1.1*Re^0.6*Sc^(1/3))*Dm/dp;%膜传质系数计算，m/s wakao
        % Kf=0.357/e*Re^0.64*Sc^0.33*Dm/dp*0.3048/60;%膜传质系数计算，m/s 博士论文
        %         Kf=0.357/e*Re^0.64*Sc^0.33*Dm/dp;%膜传质系数计算，m/s 改数据
        %Sh=2.0+1.8*Re^(1/2)*Sc^(1/3);%舍伍德数计算
        Dk=beishu_Dk*97*rpore*sqrt(y(ii+2*(N+1))/Ma);%克努森扩散系数，m^2/s
        Ds=3.871E-7*beishu_Ds0*exp(0.45*deltH/(R*y(ii+2*(N+1))));%表面扩散系数
        %Ds0=1.5E-7 Adsorption and thermal regeneration of methylene chloride vapor on an activated carbon bed
        %Ds0=7.547E-7 to 1.05664E-6 博士论文
        kg=5.4E-5*(y(ii+2*(N+1))-273.15)+0.024313;%空气的导热系数，W/(m*K)
        %kg=6E-5*y(ii+2*(N+1))+0.0062;%氮气导热系数
        %Cpg=31.15-1.357E-2*y(ii+2*(N+1))+2.68E-5*y(ii+(N+1))^2-1.168E-8*y(ii+2*(N+1))^3;%J/mol/K
        % Cpg=-1E-6*(y(ii+2*(N+1))-273.15)^2+0.1*(y(ii+2*(N+1))-273.15)+1009;%空气的比热容，J/(kg*K)
%         Cpg_O2=1000*(25.477+15.2022E-3*y(ii+2*(N+1))-5.0618E-6*y(ii+2*(N+1))*y(ii+2*(N+1))+1.3117E-9*y(ii+2*(N+1))*y(ii+2*(N+1))*y(ii+2*(N+1)))/32;
%         Cpg_N2=1000*(28.901-1.5713E-3*y(ii+2*(N+1))+8.0805E-6*y(ii+2*(N+1))*y(ii+2*(N+1))-28.7256E-9*y(ii+2*(N+1))*y(ii+2*(N+1))*y(ii+2*(N+1)))/28;
Cpg_N2=4.63e-07*y(ii+2*(N+1))^3+0.0001744*y(ii+2*(N+1))^2+0.01004*y(ii+2*(N+1))+1039;
        Cpg_O2=-1.295e-06*y(ii+2*(N+1))^3+0.0009329*y(ii+2*(N+1))^2+0.1084*y(ii+2*(N+1))+913.8;        
Cpg=0.7*Cpg_N2+0.3*Cpg_O2;%空气的比热容，J/(kg*K)
        %Cpg=(28.901-1.5713*y(ii+2*(N+1))+8.0805*y(ii+2*(N+1))*y(ii+2*(N+1))-28.7256*y(ii+2*(N+1))*y(ii+2*(N+1))*y(ii+2*(N+1)));%氮气比热容
        Pr=niandu*Cpg/kg;%普朗特数
        Kz=(7+0.5*Pr*Re)*kg;%气相有效轴向导热系数，W/(m*K) Heat transfer of gas flow through a packed bed
        
        qm=M0*exp(M1/y(ii+3*(N+1)));
        b=B0*exp(B1/y(ii+3*(N+1)));
        qe=qm*(b*(P*y(ii))^n)/(1+(b*(P*y(ii))^n));%温度T，浓度c下的平衡吸附量，mol/kg
        % Cdiffq=qm/((qm-y(ii+N+1))^2*b*P);%平衡浓度对吸附量的导数，用于求解有效扩散系数
        Cdiffq=1/n*(y(ii+N+1)/((qm-y(ii+N+1))*b*P^n))^(1/n-1)*(y(ii+N+1)/((qm-y(ii+N+1))^2*b*P^n)+1/((qm-y(ii+N+1))*b*P^n));%yun的论文的基础上根据sips模型推导
        De=Ds+Dk*ep*e*density_gas/density_particle*Cdiffq;%有效扩散系数，m^2/s
        % k=6*2*density_gas*Kf*y0*60*De/(6*density_gas*Kf*y0*dp^2+60*De*density_particle*dp*qe);%总质量传递系数，1/s
        y_xiuzheng=max(y(ii),y_min_value);
        k_daoshu=1*density_particle/Kf/density_gas/as*qe/y_xiuzheng+dp/10/De/as;
        k_daoshu=1*density_particle/Kf/density_gas/as*qe/y(ii)+dp/10/De/as;
        k_daoshu=dp/10/De/as;
        
        % k_daoshu=dp/5/De/as;
        k=1/k_daoshu*beishu;
        %    k_record(ii,round(t/delta_t))=k;
%                 k=0.00031;
        %qe_div_c=qm*b*P/(1+b*P*y(ii));%平衡吸附量除以浓度。
        %k_daoshu=density_particle/(Kf*density_gas*as)*qe_div_c+dp/(10*De*as);
        %k=1/k_daoshu;
        
        %         dy(ii)=Dz*(y(ii+1)-2*y(ii)+y(ii-1))/z^2-u*(y(ii+1)-y(ii))/z-Mg/1000/density_gas*((1-e)/e)*density_particle*k*(qe-y(ii+N+1));%传质方程
        dy(ii)=Dz*(y(ii+1)-2*y(ii)+y(ii-1))/z^2-u*(y(ii+1)-y(ii))/z-Mg/1000/density_gas/e*density_particle*k*(qe-y(ii+N+1));%传质方程
        dy(ii+N+1)=k*(qe-y(ii+N+1));%吸附速率方程
        dy(ii+2*(N+1))=-u*(y(ii+1+2*(N+1))-y(ii+2*(N+1)))/z+...
            Kz/density_gas/Cpg*(y(ii+1+2*(N+1))-2*y(ii+2*(N+1))+y(ii-1+2*(N+1)))/z^2-...
            hs*as*(1-e)/(e*density_gas*Cpg)*(y(ii+2*(N+1))-y(ii+3*(N+1)))...
            -hw*4/(db*e*density_gas*Cpg)*(y(ii+2*(N+1))-y(ii+4*(N+1)));%气体能量方程
        dy(ii+3*(N+1))=(1-e)*hs*as/(density_particle*Cps)*(y(ii+2*(N+1))-y(ii+3*(N+1)))-deltH/Cps*k*(qe-y(ii+N+1));%吸附剂能量方程
        dy(ii+4*(N+1))=hw*aw/(density_wall*Cpw)*(y(ii+2*(N+1))-y(ii+4*(N+1)))-U*aa/(density_wall*Cpw)*(y(ii+4*(N+1))-Ta);%吸附床体能量方程
        
    end
end