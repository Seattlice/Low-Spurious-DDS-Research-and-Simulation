%% CORDIC仿真程序
%% 清理工作区
clc;            %清除命令行
clear all;      %清楚工作区变量，释放空间

%% 数据参数
F1=1e8;           %信号频率
Fs=200e6;        %采样频率
P1=0;           %信号初始相位
N_caiyang = 12;
N=2^N_caiyang;         %采样点数
t = linspace(0,1/(F1), N);
A = 2^15;
ADC = A- 1 ;
T = 2 ^ 18;%仿真的实践
%% cordic ROM表
iterations = 16;  % 替换为所需的迭代次数

%% 参数设置

fre_weishu = 32; %累加器位数
fre_add = 0;
romaddr_reg = 0;
dac_data = 0;
jieduan = 12; %截断位数
P_jieduan = 7;%小数部分截断
lfm_fre = 0;
fuhao = 1;

Fc =1e8;
f0 = 5e6;

F_WORD = round(f0*2^fre_weishu/Fc);
P_WORD = 0;
%% 相位累加器
for i = 1:T     
    if fre_add + F_WORD > 2^fre_weishu -1 %%累加判断是否溢出
        fre_add = fre_add + F_WORD - 2^fre_weishu + 1;
    else
        fre_add = fre_add + F_WORD ;
    end
        s1(i) = fre_add;
      
    
    %fre_add = fre_add + F_WORD  + randi(2^(14));  % 进行数值运算
   
    s1(i) = fre_add;
    

    s2(i) = romaddr_reg;
    
     TP1 = bitshift(fre_add, -P_jieduan);
     TP2 =  bitshift(fre_add, -jieduan);   
     TP = TP1 - TP2 * 2^( jieduan - P_jieduan);
     TP = 6*TP/2^fre_weishu;
     p(i) = TP;   
  
    % 相位截断
    romaddr_reg = bitshift(fre_add, -jieduan)+ P_WORD;
    if romaddr_reg >= 2^(fre_weishu - jieduan)
        romaddr_reg = romaddr_reg  - 2^(fre_weishu - jieduan);
    end  
     
    %相幅转换器
    dac_data_sin =  sin_dds_cordic(pi*romaddr_reg/(2^(fre_weishu - jieduan)/2), iterations);
    dac_data_cos = sin_dds_cordic(pi*(romaddr_reg+262144)/(2^(fre_weishu - jieduan)/2),iterations);
    
%    dac_sin(i) = dac_data_sin + TP * dac_data_cos;
%    dac_cos(i) = dac_data_cos - TP * dac_data_sin;   
    s3(i) = dac_data_cos;
end



%% 结果进行验证
figure;
subplot(3,1,1);
T2 = (1 : T)/(Fc);
plot(T2,s3);
title('时域波形');
xlabel('时间 (秒)');
ylabel('幅度');
t3 = (1:1048576);
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
