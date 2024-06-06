%% 相位抖动法的仿真
%% 清理工作区
clc;            %清除命令行
clear all;      %清楚工作区变量，释放空间

%% ROM表
%% 数据参数
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

car =   ss ;
% for i = 1:N
%     if(car(i) < 0)
%         car(i) = 0;
%     end    
% end
%% 绘图验证
%figure;
%plot(4096*t,car)

%% 参数设置

fre_weishu = 32; %累加器位数
fre_add = 0;
romaddr_reg = 0;
dac_data = 0;
jieduan = 20; %截断位数
lfm_fre = 0;
fuhao = 1;%符号位，根据theta的象限来决定符号的正负
Fc =1e8;
f0 = 0.01e6;

F_WORD = round(f0*2^fre_weishu/Fc);
P_WORD = 0;
T = 2^20;
N = fre_weishu - jieduan;
%% 相位累加器
for i = 1:T     
    
    if fre_add + F_WORD > 2^fre_weishu -1 %%累加判断是否溢出
        fre_add = fre_add + F_WORD - 2^fre_weishu  + randi(2^(12)) + 1;
%         fre_add = 0;
    else
        fre_add = fre_add + F_WORD + randi(2^(12));
    end
        s1(i) = fre_add;
    
    %fre_add = fre_add + F_WORD  + randi(2^(12));  % 进行数值运算
    

    s1(i) = fre_add;
    
    % 相位截断
    romaddr_reg = bitshift(fre_add, -jieduan)+ P_WORD;
    if romaddr_reg >= 2^N
        romaddr_reg = romaddr_reg  - 2^N;
    end  

    s2(i) = romaddr_reg;
    
    %相幅转换器
    dac_data = car(romaddr_reg + 1);
    s3(i) = dac_data;
end



%% 结果进行验证
figure;
subplot(3,1,1);
T2 = (1 : T)/(Fc);
plot(T2,s3);grid on;
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
plot(f/1e6, Z);grid on;
title('频域波形');
xlabel('频率 (MHz)');
ylabel('幅度/db');

subplot(3,1,3);
plot(f/1e6, P1);
title('频域波形');
xlabel('频率 (MHz)');
ylabel('幅度');


