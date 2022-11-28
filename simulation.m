clear;clc

%% Set initial values
t = 1:100; % time
x = t; % x-coord of Tx
y = 9*t.^0.5; % y-coord of Tx
xMax = max(x, [], "all");
yMax = max(y, [], "all");
assert(length(t) == length(x))
assert(length(t) == length(y))

%% Set Rx values
xRx = [       0        0     xMax]; % x-coords of Rx
yRx = [       0     yMax        0]; % y-coords of Rx
assert(length(xRx) == length(yRx))

%% Solve the system
gidx = 1; % index for initial guess
rad = zeros(length(xRx)); % Distances between each Rx and the Tx. To be filled with values from Friis eqn
figure(1)
for i = 1:length(t)
    % Calculate and estimate object position
    rad(1) = pythag([xRx(1) yRx(1)], [x(i) y(i)]);
    rad(2) = pythag([xRx(2) yRx(2)], [x(i) y(i)]);
    rad(3) = pythag([xRx(3) yRx(3)], [x(i) y(i)]);
    sol = fsolve(@(pos) solvesys(pos, xRx, yRx, rad), [xRx(gidx) yRx(gidx)])
    plot(sol(1), sol(2), "b+", "markersize", 5)

    % Plot object position
    plot(x(i), y(i), "rx", "markersize", 5)

    % Format graph
    title("Animation of moving object")
    axis([0, xMax, 0, yMax])
    grid on
    hold off
    pause(0.05)
end

%% System solving function
%  @param pc - 2d-coords of Tx
%  @param x  - vector of x-coords for Rx
%  @param y  - vector of y-coords for Rx
%  @param d  - vector of Tx distances (one for each Rx)
%  @return posTx - estimated position of Tx
function posTx = solvesys(pc, x, y, d)
    posTx = [
        pc(1)^2 - 2*x(1)*pc(1) + x(1)^2 + pc(2)^2 - 2*y(1)*pc(2) + y(1)^2 - d(1)^2;...
        pc(1)^2 - 2*x(2)*pc(1) + x(2)^2 + pc(2)^2 - 2*y(2)*pc(2) + y(2)^2 - d(2)^2;...
        pc(1)^2 - 2*x(3)*pc(1) + x(3)^2 + pc(2)^2 - 2*y(3)*pc(2) + y(3)^2 - d(3)^2;...
    ];
end

%% Pythagorian theorem
function h = pythag(a, b)
    dx = a(1) - b(1);
    dy = a(2) - b(2);
    h = sqrt(dx^2 + dy^2);
end