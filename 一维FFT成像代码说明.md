# 一维FFT成像代码说明

## 1.设置基本参数

在这一部分中中需要首先对包括光速、频率、天线阵列数、天线阵列半径及涡旋波模式数等诸多参数进行设定。具体设定如下：<br>

```matlab
clc;clear;close all;
%% 基础参数
c = 299792458;              % 光速
f0 = 9e9;                   % 频率
lambda = c / f0 ;           % 波长
k = 2.0 * pi / lambda;      % 波数

arr_radius = 25 * lambda;   % 天线半径
arr_number = 16;            % 天线数量

moderange = 20;				%模式数变化范围
deltamode = 0.5;			%模式数梯度间隔
%% 反射点
target_num = 2;				%反射点数量
target_radius = [500,500] * lambda;		%各反射点距离值
target_elevation = [0.3,0.3] * pi;		%各反射点高度角
target_azimuth = [0.1,0.3] * pi;		%各反射点方位角
rcs = ones(1,target_num);				%反射点的电磁散射系数矩阵
```

# 2.回波信号生成

在这一部分中，需要按照论文所推导的多发多收（MIMO）和多发单收（MISO）两种收发模式下的回波信号公式，在程序中生成回波信号，以便后续操作。<br>

假设目标由多个理想散射点组成，以
$$
P_m(r_m,\theta_m,\phi_m)
$$
表示目标位置，对应的目标的电磁散射系数为
$$
\sigma_m(r_m,\theta_m,\phi_m),简写为\sigma_m
$$
则多发多收模式下的回波信号公式：<br>
$$
E_r(r)=N^2e^{il\pi}\sum^M_{m=1}\sigma_m\frac{e^{-i2kr_m}}{r^2_m}e^{il2\phi_m}J^2_l(ka\sin\theta_m)
$$
同理，将多发单收模式下的回波信号公式写为：
$$
E_r(r)=Ne^{il\frac{\pi}{2}}\sum^M_{m=1}\sigma_m\frac{e^{-i2kr_m}}{r^2_m}e^{il\phi_m}J_l(ka\sin\theta_m)
$$
根据以上两条公式，可以通过如下代码实现两种收发模式下回波信号的生成：<br>

```matlab
%%回波信号生成
modes = -moderange:deltamode:moderange;		%moderange设定为20，总变化范围为-20~20
echo_mimo = zeros(1,length(modes));			%声明mimo和miso的回波信号矩阵
echo_miso = zeros(1,length(modes));

for count = 1:length(modes)					%通过循环得到两种模式下的回波信号
    mode = modes(count);

    echo_mimo(count) = arr_number^2.0 * exp(1i * mode * pi) * rcs ...
                .* exp(-1i * 2.0 * k * target_radius)./(target_radius.^2.0) ...
                .* exp(1i * 2.0 * mode * target_azimuth) ...
                * besselj(mode, k*arr_radius' .*sin(target_elevation')).^2.0;

    echo_miso(count) = arr_number * exp(1i * mode * pi / 2.0) * rcs ...
                .* exp(-1i * 2.0 * k * target_radius)./(target_radius.^2.0) ...
                .* exp(1i * mode * target_azimuth) ...
                * besselj(mode, k*arr_radius' *sin(target_elevation'));
    
end
```

## 3.FFT和成像操作

由公式推导可以得出，回波信号中的模式域与方位角之间存在对偶关系，因此使用FFT方法对回波信号进行处理，来解出信号中的方位角信息，以实现一维FFT成像。<br>

在成像中需要注意的两点：首先由于MIMO回波函数中
$$
e^{il2\phi_m}
$$
的原因，经过FFT处理的MIMO回波信号中得到的方位角与正确值存在二倍关系。其次，当模式数梯度间隔（delta）由整数变为分数时，会使最终得到的FFT结果中有效的部分收缩，因此为了能够正确成像，也需要在成像时对坐标轴进行一定的变化。具体成像代码如下：<br>

```matlab
%%FFT与成像处理
N = 2048;			%FFT成像点数
echo_mimo_fft = fft(echo_mimo,N);							%对回波信号进行FFT变换
echo_mimo_y = abs(echo_mimo_fft)/max(abs(echo_mimo_fft));	%幅值归一化
echo_mimo_x = linspace(-90 , 90 , N * deltamode);			%生成对应的ｘ轴数据

%echo_miso_hilbert = hilbert((echo_miso));
echo_miso_fft = fft(echo_miso,N);							%MISO与MIMO类似
echo_miso_y = abs(echo_miso_fft)/max(abs(echo_miso_fft));
echo_miso_x = linspace(-180 , 180 , N * deltamode);


figure(1);
plot(echo_mimo_x,echo_mimo_y(1:N*deltamode));		%MIMO成像绘图
xlabel('方位角/°');
ylabel('归一化幅值');
title('MIMO成像');

figure(2);
plot(echo_miso_x,echo_miso_y(1:N*deltamode));		%MISO成像绘图
xlabel('方位角/°');
ylabel('归一化幅值');
title('MISO成像');

```

综合以上三段就可以正确实现一维FFT方法的涡旋波成像仿真。
