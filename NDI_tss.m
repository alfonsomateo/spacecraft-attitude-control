function Tc = NDI_tss(t,X,sim)
[att_c,att_dot_c] = comanded_attitude(t,sim);
if strcmp(sim.att_repr,'Euler')
    att = X(1:3);
    omega = X(4:6);
    N_theta = (1/cos(att(2)))*[cos(att(2)),sin(att(1))*sin(att(2)),cos(att(1))*sin(att(2));
        0,cos(att(1))*cos(att(2)),-sin(att(1))*cos(att(2));
        0,sin(att(1)),cos(att(1))];
    omega_NI = eul2rotm(att','ZYX')*[0,-sim.n,0]';
    
    % Outer loop
    nu_2 = -sim.Kp*(att-att_c);
    omega_c = N_theta\nu_2+omega_NI;
    
    % Inner loop
    nu_1 = -sim.Kd*(omega-omega_c);
    Tc = sim.J*nu_1+cross(omega,sim.J*omega)-sim.Td;
elseif strcmp(sim.att_repr,'quaternions')
    att = X(1:4);
    omega = X(5:7);
    N_q = 0.5*[att(4),-att(3),att(2);
        att(3),att(4),-att(1);
        -att(2),att(1),att(4);
        -att(1),-att(2),-att(3)];
    quat = [att(4); att(1:3)]';
    omega_NI = quat2rotm(quat)*[0,-sim.n,0]';
  
    % Outer loop
    Rc = [att_c(4),att_c(3),-att_c(2),-att_c(1);
        -att_c(3),att_c(4),att_c(1),-att_c(2);
        att_c(2),-att_c(1),att_c(4),-att_c(3);
        att_c(1),att_c(2),att_c(3),att_c(4)];
    qe = Rc*att;
    nu_2 = -sim.Kp*qe(1:3);
    omega_c = N_q(1:3,1:3)\nu_2+omega_NI;
    
    % Inner loop
    nu_1 = -sim.Kd*(omega-omega_c);
    Tc = sim.J*nu_1+cross(omega,sim.J*omega)-sim.Td;
end
end