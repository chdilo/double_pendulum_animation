clear, clc, close all
%% 求解双摆的常微分方程
m = 1;       % 质量，kg
L = 1;       % 长度，m
a1 = 3*pi/4; % 角1，rad
a2 = 3*pi/8; % 角2，rad
p1 = 0;      % 动量1，kg*m/s
p2 = 0;      % 动量2，kg*m/s
frameRate = 30;        % 帧率，Hz，FPS
n = 8;                 % 每帧用n张图平均
steps = n * frameRate; % 每秒步数
dur = 4;               % 持续时间，second
nFrames = dur * frameRate; % 总帧数
nSteps = steps * dur;      % 总步数
time = 0:1/steps:dur-1/steps; % 全部时间点

% 指定绝对误差容限和相对误差容限
options = odeset('AbsTol', 1e-50, 'RelTol', 1e-13);
% 求解双摆的常微分方程
[T, Y] = ode15s(@(t, x) double_pendulum(t, x, m, L), ...
    time, [a1, a2, p1, p2], options);

% plot(t, Y(:,3:4)), legend('动量1', '动量2')
%% 生成动画帧序列
width = 640;
height = 360;
wCentre = (width+1) / 2;
hCentre = (height+1) / 2;
lineWidth = round(height/100);

x1 = L*sin(Y(:,1));
y1 = -L*cos(Y(:,1));
x2 = x1 + L*sin(Y(:,2));
y2 = y1 - L*cos(Y(:,2));

imageX1 = (height/4)*x1 + wCentre;
imageY1 = -(height/4)*y1 + hCentre;
imageX2 = (height/4)*x2 + wCentre;
imageY2 = -(height/4)*y2 + hCentre;

path = uigetdir('.\', '选择动画帧序列保存路径');
Fig = waitbar(0, '正在处理帧...');
for i = 1:nFrames
    image = zeros(height, width);
    for j = 1:n
        imageTemp = zeros(height, width);
        imageTemp = insertShape(imageTemp, 'Line', ...
            [wCentre hCentre imageX1((i-1)*n+j) imageY1((i-1)*n+j);
            imageX1((i-1)*n+j) imageY1((i-1)*n+j) imageX2((i-1)*n+j) imageY2((i-1)*n+j)],...
            'Opacity', 1, ...
            'LineWidth', lineWidth, ...
            'Color', [0.8500 0.3250 0.0980; 0 0.4470 0.7410]);
        imageTemp = insertShape(imageTemp, 'FilledCircle', ...
            [wCentre hCentre lineWidth;
            imageX1((i-1)*n+j) imageY1((i-1)*n+j) lineWidth;
            imageX2((i-1)*n+j) imageY2((i-1)*n+j) lineWidth], ...
            'Opacity', 1, ...
            'Color', [0.9290 0.6940 0.1250; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330]);
        image = image + imageTemp/n;
    end
    imwrite(image, fullfile(path, sprintf('%06u.png', i)))
    waitbar(i/nFrames, Fig, ...
        sprintf('正在处理帧...%.2f%% (%u/%u)', 100*i/nFrames, i, nFrames));
end
close(Fig)

%% 生成轨迹帧序列
image1 = zeros(height, width, 'gpuArray');
image2 = zeros(height, width, 'gpuArray');

colorMap1 = hot(256);
colorMap2 = myColorBlue(256);

radius = height/1024;
brightness = (height/radius/1024)^2;
halfLife = 4;
decayRate = 0.5^(1/(halfLife*steps));

xSequence = gpuArray(1:width);
ySequence = gpuArray(1:height);

path = uigetdir('.\', '选择轨迹帧序列保存路径');
Fig = waitbar(0, '正在处理帧...');
for i = 1:nFrames
    for j = 1:n
        [xArray1, yArray1] = meshgrid(xSequence-imageX1((i-1)*n+j), ySequence-imageY1((i-1)*n+j));
        [xArray2, yArray2] = meshgrid(xSequence-imageX2((i-1)*n+j), ySequence-imageY2((i-1)*n+j));
        image1 = decayRate*image1 + brightness./(1 + (xArray1/radius).^2 + (yArray1/radius).^2);
        image2 = decayRate*image2 + brightness./(1 + (xArray2/radius).^2 + (yArray2/radius).^2);
    end
    rgbImage = ind2rgb(uint8(256*image1), colorMap1) + ind2rgb(uint8(256*image2), colorMap2);
    imwrite(rgbImage, fullfile(path, sprintf('%06u.png', i)))
    waitbar(i/nFrames, Fig, ...
        sprintf('正在处理帧...%.2f%% (%u/%u)', 100*i/nFrames, i, nFrames));
end
close(Fig)

%% 生成音频
[file, path] = uiputfile({'*.flac';'*.wav'}, '保存音频文件', 'Audio');
fs = 48e3; % 音频采样率
t = 0:1/fs:dur-1/fs;
kf = 256; % 调频灵敏度
fc = 256; % 载波频率

v1 = (6 * (2*Y(:,3)-3*cos(Y(:,1)-Y(:,2)).*Y(:,4))) ./ ...
    ((m*L^2) * (16-9*cos(Y(:,1)-Y(:,2)).^2));

v2 = (6 * (8*Y(:,4)-3*cos(Y(:,1)-Y(:,2)).*Y(:,3))) ./ ...
    ((m*L^2) * (16-9*cos(Y(:,1)-Y(:,2)).^2));

vInterp1 = interp1(T, v1, t, 'spline');
vInterp2 = interp1(T, v2, t, 'spline');

includedAngle = interp1(T, Y(:,1)-Y(:,2), t, 'spline');

f1 = abs(vInterp1);
f2 = sqrt(vInterp1.^2 + vInterp2.^2 + 2*vInterp1.*vInterp2.*cos(includedAngle));
fmax = max([f1 f2]);
f1 = f1/fmax;
f2 = f2/fmax;

fIntegral1 = cumtrapz(t, f1);
fIntegral2 = cumtrapz(t, f2);

leftAmplitude1 = interp1(T, (1-x1)/2, t, 'spline');
leftAudio1 = f1.*leftAmplitude1.*cos(2*pi * (fc * t + kf * fIntegral1))/2;
leftAmplitude2 = interp1(T, (1-x2)/2, t, 'spline');
leftAudio2 = f2.*leftAmplitude2.*cos(2*pi * (fc * t + kf * fIntegral2))/2;

leftAudio = leftAudio1 + leftAudio2;

rightAmplitude1 = interp1(T, (1+x1)/2, t, 'spline');
rightAudio1 = f1.*rightAmplitude1.*cos(2*pi * (fc * t + kf * fIntegral1))/2;
rightAmplitude2 = interp1(T, (1+x2)/2, t, 'spline');
rightAudio2 = f2.*rightAmplitude2.*cos(2*pi * (fc * t + kf * fIntegral2))/2;

rightAudio = rightAudio1 + rightAudio2;

audiowrite([path, file], [leftAudio' rightAudio'], fs, 'BitsPerSample', 24)
