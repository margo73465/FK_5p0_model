figure, hold all;

for i=1:18
    plot(1:30000, rebuild_voltage_myStim((i - 1) * 30000 + 1:i*30000, 70))
end

%%

figure, hold all

for i=1:10
    plot(1:25000,mcurveshape((i - 1) * 25000 + 1:i * 25000))
end

%%

dVmax = zeros(size(down,1) - 1,3);
for i = 1:(size(down,1) - 1)
    upstroke = dV(up(i,1):down(i+1,1));
    [M,I] = max(upstroke);
    dVmax(i,1) = M;
    dVmax(i,2) = I + up(i,1);
    dVmax(i,3) = voltage(dVmax(i,2));
end

%% Add generation numbers to evaluated solutions file

shape(:,23) = 0;

for i=2:(size(shape,1))
    shape(i,23) = floor((i - 2) / 500) + 1;
end

%%
restitution(:,23) = 0;

for i=2:(size(restitution,1))
    restitution(i,23) = floor((i - 2) / 500) + 1;
end

