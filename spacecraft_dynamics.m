function X_dot = spacecraft_dynamics(X,Tc,sim)
if strcmp(sim.att_repr,'Euler')
    att = X(1:3);
    omega = X(4:6);
    N_theta = (1/cos(att(2)))*[cos(att(2)),sin(att(1))*sin(att(2)),cos(att(1))*sin(att(2));
        0,cos(att(1))*cos(att(2)),-sin(att(1))*cos(att(2));
        0,sin(att(1)),cos(att(1))];
    omega_NI = eul2rotm(att','ZYX')*[0,-sim.n,0]';
    att_dot = N_theta*(omega-omega_NI);
    omega_dot = sim.J\(-cross(omega,sim.J*omega)+sim.Td+Tc);
    X_dot = [att_dot;omega_dot];
elseif strcmp(sim.att_repr,'quaternions')
    att = X(1:4);
    omega = X(5:7);
    N_q = 0.5*[att(4),-att(3),att(2);
        att(3),att(4),-att(1);
        -att(2),att(1),att(4);
        -att(1),-att(2),-att(3)];
    quat = [att(4); att(1:3)]';
    omega_NI = quat2rotm(quat)*[0,-sim.n,0]';
    att_dot = N_q*(omega-omega_NI);
    omega_dot = sim.J\(-cross(omega,sim.J*omega)+sim.Td+Tc);
    X_dot = [att_dot;omega_dot];
end
end