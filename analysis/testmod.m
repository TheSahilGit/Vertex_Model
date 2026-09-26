clear; clc; close all;

ontime = 50;
offtime = 50;

for it = 1:7000
    if mod(it,ontime+offtime) < ontime
        s(it) = 1;
    else
        s(it) = 0;
    end

end