%% Creates APD restitution curve for stochastic pacing using APD_90

clear APD DI

%data = objective;
%data = objective_thresh73;
%data = objective_thresh73_expDI;
%data = objective_thresh70_newDI;
% data = objective_APD90_expDI;
%data = FK4V_cell_1p0_V_1p0_GAparams_thresh73_200x200(10001:30000,2);
%data = FK4V_cell_V_extraGAparams_200x200(10001:30000,2);
%data = FK4V_cell_1p0_V_1p0_GAparams_thresh70_200x200(10001:30000,2);
%data = FK4V_cell_1p0_V_1p0_GAparams_thresh73_expDI_200x200(10001:30000,2);
%data = FK4V_cell_1p0_V_1p0_thresh70_newDI_200x200(10001:30000,2);

% stim_times = APD90_expDI_stim_times;

len = size(data,1);
dVmax = 0;
V1 = data(1);
%V_low = min(data);
n = 1;
V_90 = -100.0;
V_low = 100.0;
DIstart = 0;
j = 1;
APD = zeros(60,3); 
DI = zeros(60,3);

for i=2:len 
    
    V1 = data(i);
    dV = data(i) - data(i-1);
    if dV >= dVmax
        disp(i);
        dVmax = dV;
        %disp(dVmax);
        j = 0;
    elseif dV < 0 && j==0
        disp('hello');
        APDstart = i;
        V_90 = data(i-1) - (data(i-1) - V_low)*.90;
        j = j+1;
        DI(n,1) = i;
        DI(n,2) = DIstart;
        DI(n,3) = (i-DIstart);
    end
    
    if data(i) < V_90 && data(i-1) >= V_90
        APD(n,1) = i;
        APD(n,2) = APDstart;
        APD(n,3) = (i - APDstart);
        dVmax = 0.0;
        n = n+1;
        DIstart = i;
    end
    
end

%len = size(APD,1);
%grandi_APD90_expDI_stochastic = [DI(2:len,3),APD(2:end,3)];

%figure, hold on;
%plot(DI(2:len,3),APD(2:end,3),'.b');

%% Creates APD restitution curve for stochastic pacing using APD_90

clear APD DI

%data = objective;
%data = objective_thresh73;
%data = objective_thresh73_expDI;
%data = objective_thresh70_newDI;
%data = objective_APD90_expDI;
%data = FK4V_cell_1p0_V_1p0_GAparams_thresh73_200x200(10001:30000,2);
%data = FK4V_cell_V_extraGAparams_200x200(10001:30000,2);
%data = FK4V_cell_1p0_V_1p0_GAparams_thresh70_200x200(10001:30000,2);
%data = FK4V_cell_1p0_V_1p0_GAparams_thresh73_expDI_200x200(10001:30000,2);
%data = FK4V_cell_1p0_V_1p0_thresh70_newDI_200x200(10001:30000,2);
data = FK4V_cell_1p0_V_1p0_APD90_expDI_200x200(10001:30000,2);

stim_times = APD90_expDI_stim_times/100;

len = size(data,1);
V_max = -100.0;
n = 2;
V_90 = -100.0;
V_low = 100.0;
APD = zeros(numel(stim_times),3);
DI = zeros(numel(stim_times),1);

for i=2:len 
    
    if n <= numel(stim_times)
        
        if i == floor(stim_times(n) - 1)
            V_low = data(i);
            V_max = -100.0;
        end
        
        if data(i) >= V_max
            V_max = data(i);
            V_90 = V_max - (V_max - V_low)*0.90;
            APDstart = i;
        end
        
        if data(i) < V_90 && data(i-1) >= V_90 && i - stim_times(n) > 50
            APD(n,1) = i;
            APD(n,2) = APDstart;
            APD(n,3) = (i - APDstart);
            n = n+1;
            disp(n);
            if n <= numel(stim_times)
                DI(n) = stim_times(n) - i;
            end
        end
    end
end

len = size(APD,1);
FK4V_APD90_expDI_stochastic = [DI(2:len),APD(2:end,3)];

figure, hold on;
plot(FK4V_APD90_expDI_stochastic(:,1),FK4V_APD90_expDI_stochastic(:,2),'.b');


%% Creates APD restitution curve for stochastic pacing using a threshold

data = objective;
%data = objective_thresh70;
%data = FK4V_cell_1p0_V_1p0_GAparams_thresh73_200x200(10001:30000,2);

threshold = -70.0;
len = size(data,1);
n = 1;
DIstart = 0;
clear DI2 APD2

for i=2:len 
    
    if data(i) > threshold && data(i-1) < threshold
        APDstart = i;
        DI2(n,1) = i;
        DI2(n,2) = DIstart;
        DI2(n,3) = (i-DIstart)/100;
        stim_time_TEST(n) = i/100;
    end
    
    if data(i) < threshold && data(i-1) > threshold
        APD2(n,1) = i;
        APD2(n,2) = APDstart;
        APD2(n,3) = (i - APDstart)/100;
        DIstart = i;
        n=n+1;
    end
end

figure
len = size(APD,1);
plot(DI2(2:len,3),APD2(2:end,3),'.b');


