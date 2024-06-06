%% Cardarilli
%% 清理工作区
clc;            %清除命令行
clear all;      %清楚工作区变量，释放空间

%% ROM表
%% 数据参数

F1= 1e8;           %信号频率
Fs=2^12;        %采样频率
T4 = 1/(F1);   %周期的时间
N=2^12;         %累加器位数/采样点数
P1=0;           %信号初始相位
t = linspace(0, T4, N);%在2pi中的前二分之pi中生成四分之一个rom单位用来存储波形数据
T = 2 ^20;%仿真的实践

%% 参数设置

fre_weishu = 16; %累加器位数
fre_add = 0;
romaddr_reg = 0;
dac_data = 0;
jieduan = 0; %截断位数
lfm_fre = 0;
fuhao = 1;

Fc =1e8;
f0 = 5e6;

F_WORD = round(f0*2^fre_weishu/Fc);
% F_WORD = 429219;
P_WORD = 0;
%% ROM表格的生成
A = 2^0;
ADC = A - 1; 

N = 16;
for j = 1 : 2^(N/2) 
    i = j - 1;
%     sina(j) = round(A * sin(2*i*pi/(2^(N/2 ))));
%     cosa(j) = round(A * cos(2*i*pi/(2^(N/2 ))) );
%     sinb(j) = round(A *  sin(2*i*pi/(2^N)));
%     cosb(j) = round(A *  cos(2*i*pi/(2^N)) );
    
    sina(j) =  sin(2*i*pi/(2^(N/2 )));
    cosa(j) =  cos(2*i*pi/(2^(N/2 )));
    sinb(j) =  sin(2*i*pi/(2^N));
    cosb(j) =  cos(2*i*pi/(2^N));    
% % %     
end
%plot(sinb);


%% 绘图验证
%figure;
%plot(4096*t,car)


%% 相位累加器
for i = 1 : T
    
    if fre_add + F_WORD > 2^fre_weishu -1 %%累加判断是否溢出
%         fre_add = 0;
        fre_add = fre_add + F_WORD - 2^fre_weishu + 1;
    else
        fre_add = fre_add + F_WORD;
    end
        s1(i) = fre_add;
    
    % 相位截断
    %romaddr_reg = bitshift(fre_add, -(fre_weishu - jieduan))+ P_WORD; 
    romaddr_reg = fre_add;
    I = bitshift(romaddr_reg, -(N/2 ));
    s2(i) = romaddr_reg;
    F = romaddr_reg - (2^(fre_weishu - N/2))*I;

    if I == 2^(N / 2 )
        I = 0;

    end
 
    if F == 2^(N / 2 )
        F = 0;
    end

    dac_sin(i) =   (sina(I + 1) * cosb(F + 1) + cosa(I + 1)*sinb(F + 1))/A ;

    dac_cos(i) =   (cosa(I + 1) * cosb(F + 1) + sinb(F + 1)*sina(I + 1))/A ;
end
%% 结果进行验证
s3 = dac_sin;
figure;
subplot(3,1,1);
T2 = (1 : T)/(Fc);
plot(T2,s3);grid on;
title('时域波形');
xlabel('时间 (秒)');
ylabel('幅度');
Y = fft(s3);  % 计算离散傅里叶变换
P2 = abs(Y/T);  % 计算幅度谱
P1 = P2(1:T/2+1);  % 取一半的频谱（单侧频谱）

P1(2:end-1) = 2*P1(2:end-1);  % 倍频谱幅度（除直流分量和Nyquist频率）
Z = 20*log10(P1);
f = Fc*(0:(T/2))/T;  % 构建频率向量

subplot(3,1,2);
plot(f/1e6, Z);grid on;
title('频域波形');
xlabel('频率 (MHz)');
ylabel('幅度/db');

subplot(3,1,3);
plot(f/1e6, P1);grid on;
title('频域波形');
xlabel('频率 (MHz)');
ylabel('幅度');
ylabel('幅度');






