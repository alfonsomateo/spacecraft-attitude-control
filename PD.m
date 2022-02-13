function nu = PD(t,X,sim)
[att_c,att_dot_c] = comanded_attitude(t,sim);
if strcmp(sim.att_repr,'Euler')
    att = X(1:3);
    omega = X(4:6);
    
    % Obtain the derivatives from the angular velocity
    N_theta = (1/cos(att(2)))*[cos(att(2)),sin(att(1))*sin(att(2)),cos(att(1))*sin(att(2));
        0,cos(att(1))*cos(att(2)),-sin(att(1))*cos(att(2));
        0,sin(att(1)),cos(att(1))];
    omega_NI = eul2rotm(att','ZYX')*[0,-sim.n,0]';
    att_dot = N_theta*(omega-omega_NI);
    
    % PD control
    nu = -sim.Kp*wrapToPi(att-att_c)-sim.Kd*(att_dot-att_dot_c);
elseif strcmp(sim.att_repr,'quaternions')
    att = X(1:4);
    omega = X(5:7);
    % Obtain the derivatives from the angular velocity
    N_q = 0.5*[att(4),-att(3),att(2);
        att(3),att(4),-att(1);
        -att(2),att(1),att(4);
        -att(1),-att(2),-att(3)];
    quat = [att(4); att(1:3)]';
    omega_NI = quat2rotm(quat)*[0,-sim.n,0]';
    att_dot = N_q*(omega-omega_NI);
    
    % PD control
    Rc = [att_c(4),att_c(3),-att_c(2),-att_c(1);
        -att_c(3),att_c(4),att_c(1),-att_c(2);
        att_c(2),-att_c(1),att_c(4),-att_c(3);
        att_c(1),att_c(2),att_c(3),att_c(4)];
    Rc_dot = [att_dot_c(4),att_dot_c(3),-att_dot_c(2),-att_dot_c(1);
        -att_dot_c(3),att_dot_c(4),att_dot_c(1),-att_dot_c(2);
        att_dot_c(2),-att_dot_c(1),att_dot_c(4),-att_dot_c(3);
        att_dot_c(1),att_dot_c(2),att_dot_c(3),att_dot_c(4)];
    qe = Rc*att;
    qe_dot = Rc_dot*att+Rc*att_dot;
    nu = -sim.Kp*qe(1:3)-sim.Kd*qe_dot(1:3);
end
end