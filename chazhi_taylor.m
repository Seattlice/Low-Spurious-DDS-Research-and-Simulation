%% 插值法
%% 清理工作区
clc;            %清除命令行
clear all;      %清楚工作区变量，释放空间

%% ROM表
%% 数据参数 
%对二次项右移了十一位，使得从26位变成了15位，在与sin相加
F1=1e8;           %信号频率
Fs=200e6;        %采样频率
P1=0;           %信号初始相位
N_caiyang = 12;
N=2^N_caiyang;         %采样点数
t = linspace(0,1/(F1), N);
A = 2^0;
ADC = A-1; %f0 = fc * K / (2 ^N)

%% sin
ss=sin(2*pi*F1*t+pi*P1/180);
cs=cos(2*pi*F1*t+pi*P1/180);
% 
% ss =  round(A * ss + ADC);
% cs = round(A * cs + ADC);


ss =  (A * ss );
cs = (A * cs );

N_jiange = 4;
j = 0;
for i = 1 : 2^N_caiyang
    if mod(i,2^N_jiange) == 0
        j = j + 1;
        sina(j) = ss(i - 2^N_jiange + 2);
        cosa(j) = cs(i - 2^N_jiange  + 2);
    end
end

%% 参数设置

fre_add_weishu = 32;
romaddr_reg = 0;
dac_data = 0;
jieduan = 20; %截断位数
lfm_fre = 0;
fre_add = 0;
T = 2^20;
N = fre_add_weishu - jieduan;

Fc =1e8;
f0 = 0.01e6;

F_WORD = round(f0*2^fre_add_weishu/Fc);
% F_WORD = 429219;
P_WORD = 0;
%% 相位累加器
for i = 1:T     
    if fre_add + F_WORD > 2^fre_add_weishu - 1 %%累加判断是否溢出
%         fre_add = 0;
        fre_add = fre_add + F_WORD - 2^fre_add_weishu + 1 ;
    else
        fre_add = fre_add + F_WORD;
    end
    
    %fre_add = fre_add + F_WORD  + randi(2^(14));  % 进行数值运算
%     fre_add = fre_add + F_WORD;
    s1(i) = fre_add;
    
    % 相位截断
    romaddr_reg = bitshift(fre_add, -jieduan)+ P_WORD;
    if romaddr_reg >= 2^N
        romaddr_reg = romaddr_reg  - 2^N;
    end  

    s2(i) = romaddr_reg;
    romaddr_reg1  = bitshift(romaddr_reg, -(fre_add_weishu - jieduan - (N_caiyang - N_jiange)));
    x_dert(i) = (romaddr_reg - romaddr_reg1 * 2^(fre_add_weishu - jieduan - (N_caiyang - N_jiange)) )*2*pi/2^fre_add_weishu;
    cosb(i) =x_dert(i)* cosa(romaddr_reg1 + 1);
    %     disp(x_dert);
    %相幅转换器
    dac_data = sina(romaddr_reg1 + 1) + cosb(i) - x_dert(i)*x_dert(i)*sina(romaddr_reg1 + 1)/2;
    s3(i) = dac_data;

end

%% 结果进行验证
figure;
subplot(3,1,1);
T2 = (1 : T)/(Fc);
plot(T2,s3);
title('时域波形');
xlabel('时间 (秒)');
ylabel('幅度');



Y = fft(s3);  % 计算离散傅里叶变换
P2 = abs(Y/T);  % 计算幅度谱
P1 = P2(1:T/2+1);  % 取一半的频谱（单侧频谱）

P1(2:end-1) = 2*P1(2:end-1);  % 倍频谱幅度（除直流分量和Nyquist频率）
Z = 20*log10(P1);
f = F1*(0:(T/2))/T;  % 构建频率向量

subplot(3,1,2);
plot(f/1e6, Z);
title('频域波形');
xlabel('频率 (MHz)');
ylabel('幅度/db');

subplot(3,1,3);
plot(f/1e6, P1);
title('频域波形');
xlabel('频率 (MHz)');
ylabel('幅度');
