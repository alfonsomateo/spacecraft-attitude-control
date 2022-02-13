function [att_c,att_dot_c] = comanded_attitude(t,sim)
sz = size(sim.att_commands_Euler);
comp = sim.att_commands_t(2:end)-t;
for i = 1:sz(1)
    if comp(i) > 0
        if strcmp(sim.att_repr,'Euler')
            att_c = sim.att_commands_Euler(i,:)';
            att_dot_c = zeros(3,1);
        elseif strcmp(sim.att_repr,'quaternions')
            att_c = sim.att_commands_quaternions(i,:)';
            att_dot_c = zeros(4,1);
        end
        return;
    end
end
if strcmp(sim.att_repr,'Euler')
    att_c = sim.att_commands_Euler(end,:)';
    att_dot_c = zeros(3,1);
elseif strcmp(sim.att_repr,'quaternions')
    att_c = sim.att_commands_quaternions(end,:)';
    att_dot_c = zeros(4,1);
end
end