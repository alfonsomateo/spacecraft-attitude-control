function Tc = NDI(t,X,sim)
if strcmp(sim.att_repr,'Euler')
    att = X(1:3);
    omega = X(4:6);
    N_omega = [(cos(att(1))*omega(2)-sin(att(1))*omega(3))*tan(att(2)),(sin(att(1))*omega(2)+cos(att(1))*omega(3))/(cos(att(2))^2),0;
        -sin(att(1))*omega(2)-cos(att(1))*omega(3),0,0;
        (cos(att(1))*omega(2)-sin(att(1))*omega(3))/cos(att(2)),(sin(att(1))*omega(2)+cos(att(1))*omega(3))*tan(att(2))/cos(att(2)),0];
    N_theta = (1/cos(att(2)))*[cos(att(2)),sin(att(1))*sin(att(2)),cos(att(1))*sin(att(2));
        0,cos(att(1))*cos(att(2)),-sin(att(1))*cos(att(2));
        0,sin(att(1)),cos(att(1))];
    omega_NI = eul2rotm(att','ZYX')*[0,-sim.n,0]';
    dNdx = [N_omega, N_theta];
    A = dNdx*[zeros(3);sim.J^-1];
    b = dNdx*[N_theta*(omega-omega_NI);sim.J\(-cross(omega,sim.J*omega)+sim.Td)];
    nu = PD(t,X,sim);
    Tc = A\(nu-b);
elseif strcmp(sim.att_repr,'quaternions')
    att = X(1:4);
    omega = X(5:7);
    N_omega = 0.5*[0,omega(3),-omega(2),omega(1);
        -omega(3),0,omega(1),omega(2);
        omega(2),-omega(1),0,omega(3);
        -omega(1),-omega(2),-omega(3),0];
    N_q = 0.5*[att(4),-att(3),att(2);
        att(3),att(4),-att(1);
        -att(2),att(1),att(4);
        -att(1),-att(2),-att(3)];
    quat = [att(4); att(1:3)]';
    omega_NI = quat2rotm(quat)*[0,-sim.n,0]';
    dNdx = [N_omega, N_q];
    A = dNdx*[zeros(4,3);sim.J^-1];
    b = dNdx*[N_q*(omega-omega_NI);sim.J\(-cross(omega,sim.J*omega)+sim.Td)];
    nu = [PD(t,X,sim);0];
    Tc = A\(nu-b);
end
end