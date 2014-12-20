figure, hold all;

for i=1:18
    plot(1:30000, rebuild_voltage_myStim((i - 1) * 30000 + 1:i*30000, 70))
end

%%

figure, hold all

for i=1:10
    plot(1:25000,mcurveshape((i - 1) * 25000 + 1:i * 25000))
end