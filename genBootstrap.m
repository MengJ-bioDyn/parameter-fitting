function  Yboot = genBootstrap(Y)

% purpose: overlapping block bootstrap for a matrix of time series
% Y=[time, variable 1, variable 2, ... variable n];

[tsn colN] = size(Y);

bootid = 1 + floor(tsn*rand(1, tsn));

sorted_bootid = sort(bootid);

Yboot = Y(sorted_bootid, :);

