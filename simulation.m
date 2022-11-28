clear;clc

%% Set Rx values
xRx = [       0        3       10]; % x coord
yRx = [       0        8        5]; % y coord
rad = [sqrt(32) sqrt(17) sqrt(37)]; % Distances between each Rx and the Tx (usually calculated by Friis transmission eqn)

%% Verify vector lengths
assert(length(xRx) == length(yRx))
assert(length(xRx) == length(rad))

%% Solve the system
gidx = 1; % index for initial guess
sol = fsolve(@(pos) solvesys(pos, xRx, yRx, rad), [xRx(gidx) yRx(gidx)])

%% System solving function
%  @param pc - 2d-coords of Tx
%  @param x  - vector of x-coords for Rx
%  @param y  - vector of y-coords for Rx
%  @param d  - vector of Tx distances (one for each Rx)
function posTx = solvesys(pc, x, y, d)
    posTx = [
        pc(1)^2 - 2*x(1)*pc(1) + x(1)^2 + pc(2)^2 - 2*y(1)*pc(2) + y(1)^2 - d(1)^2;...
        pc(1)^2 - 2*x(2)*pc(1) + x(2)^2 + pc(2)^2 - 2*y(2)*pc(2) + y(2)^2 - d(2)^2;...
        pc(1)^2 - 2*x(3)*pc(1) + x(3)^2 + pc(2)^2 - 2*y(3)*pc(2) + y(3)^2 - d(3)^2;...
    ];
end