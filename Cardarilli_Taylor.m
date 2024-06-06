%% Cardarilli
%% 清理工作区
clc;            %清除命令行
clear all;      %清楚工作区变量，释放空间

%% ROM表
%% 数据参数


T = 2 ^12;%仿真的实践

%% 参数设置

fre_weishu = 32; %累加器位数
fre_add = 0;
romaddr_reg = 0;
dac_data = 0;
jieduan = 14; %截断位数
P_jieduan =9;%小数部分截断
lfm_fre = 0;
fuhao = 1;

Fc =50e6;
f0 = 50;

F_WORD = round(f0*2^fre_weishu*100000/(Fc));
% F_WORD = 429219;
P_WORD = 0;%2^16


%% ROM表格的生成
A = 2^8;
ADC = A - 1; 

N = 18;

for j = 1 : 2^(N/2-1) 
        i = j - 1;
    sina(j) = round(A * sin(2*i*pi/(2^(N/2 + 1))));
%     cosa(j) = round(A * cos(2*i*pi/(2^(N/2 + 1))));
    cosa(j) = sina(j);
    sinb(j) = round(A/2 *  sin(2*i*pi/(2^N)));      %丢弃了一位符号位
    cosb(j) = round(A/2 * (1 - cos(2*i*pi/(2^N)))) ;
%    

%     sina(j) =   A * sin(2*i*pi/(2^(N/2 + 1)));
% %     cosa(j) =   A * cos(2*i*pi/(2^(N/2 + 1)));
%     cosa(j) = A*sina(j);
%     sinb(j) =   A * sin(2*i*pi/(2^N));
%     cosb(j) =   A * (1 - cos(2*i*pi/(2^N)));    
%     
end
%plot(sinb);


%% 绘图验证
%figure;
%plot(4096*t,car)


%% 相位累加器
for i = 1 : T
    
    if fre_add + F_WORD > 2^fre_weishu -1 %%累加判断是否溢出
        fre_add = fre_add + F_WORD - 2^fre_weishu + 1;
    else
        fre_add = fre_add + F_WORD ;
    end
        s1(i) = fre_add;
      
        
     TP1 = bitshift(fre_add, -P_jieduan);
     TP2 =  bitshift(fre_add, -jieduan);   
     TP3 = TP1 - TP2 * 2^(jieduan - P_jieduan);
     TP = 6*TP3/2^32;
     p(i) = TP3;   
      
        
   % 相位截断
    romaddr_reg = bitshift(fre_add, -jieduan)+ P_WORD;
    if romaddr_reg >= 2^N
        romaddr_reg = romaddr_reg  - 2^N;
    end  
        
    %符号判断
    fh = bitshift(romaddr_reg, - (N - 2)) ;
    romaddr_reg = romaddr_reg - fh * (2 ^ (N - 2));
    fuhao = fh + 1;
    s4(i) = fuhao;
    % 相位截断
    %romaddr_reg = bitshift(fre_add, -(fre_weishu - jieduan))+ P_WORD; 
    
    
    
    I = bitshift(romaddr_reg, -(N/2 - 1));
    s2(i) = romaddr_reg ;
    F = romaddr_reg - (2^(N - N/2 - 1))*I;

    if I == 2^(N / 2 -1)
        I = 0;

    end

    if F == 2^(N / 2 -1)
        F = 0;
    end
    aO = 0;
    if fuhao == 1
        dac_sin(i) = ( sina(I + 1) - sina(I + 1) * cosb(F + 1)/A) + cosa(2^(N / 2 -1) -I )*sinb(F + 1)/A ;    
        dac_cos(i) = cosa(2^(N / 2 -1) -I ) - cosa(2^(N / 2 -1) -I ) * cosb(F + 1)/A - sinb(F + 1)*sina(I + 1)/A ;
    elseif fuhao == 2
        I1 = 2^(N / 2 -1) - I;
        F1 = 2^(N / 2 -1) - F;
        % F1 = F + 1;
        dac_sin(i) = (sina(I1 ) - sina(I1 ) * cosb(F1)/A + cosa(2^(N / 2 -1) -I1 +  1)*sinb(F1)/A);    
        dac_cos(i) = aO*ADC - (cosa(2^(N / 2 -1) -I1 + 1 ) - cosa(2^(N / 2 -1) -I1 +1) * cosb(F1)/A - sinb(F1 )*sina(I1)/A);
    elseif fuhao == 3
        I1 = I + 1;
        F1 = F + 1;
        dac_sin(i) = aO*ADC - (sina(I1 ) - sina(I1 ) * cosb(F1)/A + cosa(2^(N / 2 -1) -I1 + 1 )*sinb(F1)/A);    
        dac_cos(i) = aO*ADC - (cosa(2^(N / 2 -1) -I1 + 1 ) - cosa(2^(N / 2 -1) -I1 + 1 ) * cosb(F1)/A- sinb(F1 )*sina(I1)/A);
    else
        I1 = 2^(N / 2 -1) - I;
       F1 = 2^(N / 2 -1) - F;        
%        F1 = F + 1;
        dac_sin(i) = aO*ADC - (sina(I1 ) - sina(I1) * cosb(F1)/A + cosa(2^(N / 2 -1) -I1 + 1)*sinb(F1)/A);    
        dac_cos(i) = cosa(2^(N / 2 -1) -I1 + 1 ) - cosa(2^(N / 2 -1) -I1 + 1 ) * cosb(F1)/A- sinb(F1 )*sina(I1)/A;
    end
        dac_sin(i) = dac_sin(i) + TP * dac_cos(i) - dac_sin(i)*(TP^2)/4 ;
        dac_cos(i) = dac_cos(i) - TP * dac_sin(i) + dac_cos(i)*(TP^2)/4;    
end
%% 结果进行验证
s3 = dac_sin;
figure;
subplot(3,1,1);

T1 = (1 : T)/(Fc);%Fc = 512Hz
plot(T1,s3);
title('时域波形');
xlabel('时间 (秒)');
ylabel('幅度');


Y = fft(s3);  % 计算离散傅里叶变换

%N = 10000;


P2 = abs(Y/T);  % 计算幅度谱


%取全部的频谱用下方注释掉的代码
% P1 = P2(1:T);  % 
% f = Fc*(0:(T-1))/T;  % 构建频率向量

%取一半的频谱用下方的代码
P1 = P2(1:T/2+1);  % 取一半的频谱（单侧频谱）
f = Fc*(0:(T/2))/T;  % 构建频率向量

Z = 20*log10(P1);

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
ylabel('幅度');

