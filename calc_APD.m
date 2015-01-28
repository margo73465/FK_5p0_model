function [up, down, V_90, dVmax, APD_DI] = calc_APD(voltage) 

dt = 0.1;
dV = diff(voltage) / dt;
sign_dV(dV >= 0) = 1;
sign_dV(dV < 0) = 0;
inflection_points = diff(sign_dV)';

y = abs(inflection_points) .* voltage(2:end - 1);

up_indices = find(inflection_points > 0);
down_indices = find(inflection_points < 0);

up = [up_indices + 1, y(up_indices)];
down = [down_indices + 1, y(down_indices)];

if size(up,1) > size(down,1)
    up = up(1:size(down,1),:);
elseif size(down,1) > size(up,1)
    down = down(1:size(up,1),:);
end

good_points = up(:,1) - down(:,1) > 30;
up = up(good_points, :);
down = down(good_points, :);

V_low = 0.04966;
V_90 = down(:,2) - (down(:,2) - V_low) * 0.9;
V_90_index = zeros(size(V_90,1),1);

dVmax = zeros(size(down,1) - 1, 3);
for i = 1:size(down,1)
    
    repol = voltage(down(i,1):up(i,1));
    repol(repol > V_90(i)) = 1;
    repol(repol < V_90(i)) = 0;
    
    if isempty(find(diff(repol),1))
        V_90_index(i) = up(i,1);
    else
        V_90_index(i) = find(diff(repol),1) + down(i,1);
    end
    
    if (i < size(down,1))
        upstroke = dV(up(i,1):down(i + 1, 1));
        [M,I] = max(upstroke);
        dVmax(i,1) = M;
        dVmax(i,2) = I + up(i,1);
        dVmax(i,3) = voltage(dVmax(i,2));
    end
    
end


V_90 = [V_90_index, V_90];

APD = (V_90(2:end,1) - down(2:end,1)) * dt;
% DI = (up(:,1) - V_90(:,1)) * dt;
DI = (dVmax(:,2) - V_90(1:end-1,1)) * dt;

APD_DI = [DI, APD];

end
