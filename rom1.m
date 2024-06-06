%% 正弦幅度相位差法波形的对称性
%% 清理工作区
clc;            %清除命令行
clear all;      %清楚工作区变量，释放空间

%% ROM表
%% 数据参数

F1= 1e8;           %信号频率
Fs=2^12;        %采样频率
T4 = 1/(4*F1);   %周期的时间
N=2^12;         %累加器位数/采样点数
P1=0;           %信号初始相位
t = linspace(0, T4, N/4);%在2pi中的前二分之pi中生成四分之一个rom单位用来存储波形数据

A = 2^0;
ADC = A-1; 
T = 2 ^ 20;%仿真的实践
%% sin
ss=sin(2*pi*F1*t+pi*P1/180);
%plot(ss);
car =   ss  ;
% for i = 1:N/4
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
jieduan = 18; %截断位数
lfm_fre = 0;
fuhao = 1;%符号位，根据theta的象限来决定符号的正负

Fc =1e8;
f0 = 0.01e6;

F_WORD = round(f0*2^fre_weishu/Fc);
% F_WORD = 8919876;
P_WORD = 0;

%% 相位累加器
for i = 1:T      
  
    if fre_add + F_WORD > 2^fre_weishu -1 %%累加判断是否溢出
        fre_add = fre_add + F_WORD - 2^fre_weishu + 1 ;
    else
        fre_add = fre_add + F_WORD;
    end
    %fre_add = fre_add + F_WORD + randi(2^(jieduan)) - 1;  % 进行数值运算,相位抖动注入
    s1(i) = fre_add;
    
    % 相位截断
    romaddr_reg = bitshift(fre_add, -jieduan)+ P_WORD; 

% pi/2 - N/4                   ,     pi - 2048  - > N/2 - theta 
%3/2 * pi - 3072 - > theta - N/2,   2*pi - 4096 - > N - theta

    if romaddr_reg >= N
        for j = 1 : 1000000
            romaddr_reg = romaddr_reg - N + 1;
               if romaddr_reg <= N
                  break;
               end
        end
    end
if(romaddr_reg) == 2^N/4
    romaddr_reg = 1;
end
    %第4象限
if (romaddr_reg > 3/4*N) && (romaddr_reg <= N)
        romaddr_reg = N - romaddr_reg + 1;
        dac_data = 0*ADC - car(romaddr_reg);
        fuhao = 4;
    %第3象限    
elseif romaddr_reg > N/2 && romaddr_reg <= 3/4*N 
        romaddr_reg = 3/4*N - romaddr_reg;
        dac_data = 0*ADC -car(N/4 - romaddr_reg);
        fuhao = 3;        
    %第2象限    
elseif romaddr_reg <= N/2 && romaddr_reg > N/4
        romaddr_reg =  N/2 - romaddr_reg ;
        dac_data = car( romaddr_reg +1) +0*ADC;
        fuhao = 2;
    %第1象限    
elseif romaddr_reg <= N/4  && romaddr_reg > 0
        romaddr_reg = romaddr_reg  ;
         dac_data =car(romaddr_reg)+0*ADC; 
        fuhao = 1;
end    
    
    s2(i) = romaddr_reg;
    %相幅转换器
%    disp(fuhao);
    s3(i) = dac_data ;
end

%绘图观察累加器中的每个变量的变化
%plot(s1);hold on;
%plot(s2);
%figure;
%plot(s2);


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


