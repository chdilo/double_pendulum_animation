function dy = double_pendulum(~, y, m, L)
% 双摆的微分方程
% m：质量
% L：长度

g = 9.8;   % 重力加速度
a1 = y(1); % 角1
a2 = y(2); % 角2
p1 = y(3); % 动量1
p2 = y(4); % 动量2

% 导数
da1 = (6 * (2*p1-3*cos(a1-a2)*p2)) / ...
    ((m*L^2) * (16-9*cos(a1-a2)^2));
da2 = (6 * (8*p2-3*cos(a1-a2)*p1)) / ...
    ((m*L^2) * (16-9*cos(a1-a2)^2));
dp1 = -1/2*m*L^2 * (da1*da2*sin(a1-a2)+3*g/L*sin(a1));
dp2 = -1/2*m*L^2 * (-da1*da2*sin(a1-a2)+g/L*sin(a2));

dy(1) = da1; % 角1对时间的导数
dy(2) = da2; % 角2对时间的导数
dy(3) = dp1; % 动量1对时间的导数
dy(4) = dp2; % 动量2对时间的导数

dy = dy(:);
end