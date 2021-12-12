clc;clear;close all;
%% 基础参数
c = 299792458;              % 光速
f0 = 9e9;                   % 频率
lambda = c / f0 ;           % 波长
k = 2.0 * pi / lambda;      % 波数

arr_radius = 25 * lambda;   % 天线半径
arr_number = 16;            % 天线数量

moderange = 20;
deltamode = 0.5;
%% 反射点
target_num = 2;
target_radius = [500,500] * lambda;
target_elevation = [0.3,0.3] * pi;
target_azimuth = [0.1,0.3] * pi;
rcs = ones(1,target_num);

%%回波信号生成
modes = -moderange:deltamode:moderange;
echo_mimo = zeros(1,length(modes));
echo_miso = zeros(1,length(modes));

for count = 1:length(modes)
    mode = modes(count);

    echo_mimo(count) = arr_number^2.0 * exp(1i * mode * pi) * rcs ...
                .* exp(-1i * 2.0 * k * target_radius) ./ (target_radius .^ 2.0) ...
                .* exp(1i * 2.0 * mode * target_azimuth) ...
                * besselj(mode, k * arr_radius' .* sin(target_elevation')) .^ 2.0;

    echo_miso(count) = arr_number * exp(1i * mode * pi / 2.0) * rcs ...
                .* exp(-1i * 2.0 * k * target_radius) ./ (target_radius .^ 2.0) ...
                .* exp(1i * mode * target_azimuth) ...
                * besselj(mode, k * arr_radius' * sin(target_elevation'));
    
end

%%FFT与成像处理
N = 2048;
echo_mimo_fft = fft(echo_mimo,N);
echo_mimo_y = abs(echo_mimo_fft)/max(abs(echo_mimo_fft));
echo_mimo_x = linspace(-90 , 90 , N * deltamode);

%echo_miso_hilbert = hilbert((echo_miso));
echo_miso_fft = fft(echo_miso,N);
echo_miso_y = abs(echo_miso_fft)/max(abs(echo_miso_fft));
echo_miso_x = linspace(-180 , 180 , N * deltamode);


figure(1);
plot(echo_mimo_x,echo_mimo_y(1:N*deltamode));
xlabel('方位角/°');
ylabel('归一化幅值');
title('MIMO成像');

figure(2);
plot(echo_miso_x,echo_miso_y(1:N*deltamode));
xlabel('方位角/°');
ylabel('归一化幅值');
title('MISO成像');

