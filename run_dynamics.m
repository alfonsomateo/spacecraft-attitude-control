function ADCS_results = run_dynamics(sim)

% Initialise simulation loop

clock_cycle = 0;
control_cycle = 0;

final_cycle = sim.tf/sim.dt;
reset_counter = sim.ts/sim.dt;

if strcmp(sim.att_repr,'Euler')
    att = sim.initial_att';
    omega = zeros(3,1);
    raw_results = zeros(final_cycle+1,10);
elseif strcmp(sim.att_repr,'quaternions')
    quat = eul2quat(sim.initial_att,'ZYX');
    att = [quat(2:4) quat(1)]';
    omega = zeros(3,1);
    raw_results = zeros(final_cycle+1,11);
end
X = [att;omega];
X_sensor = [att;omega];
Tc = zeros(3,1);

% Select type of control

if strcmp(sim.control_type,'PD')
    control_law = @(t,X_sensor,Tc) PD(t,X_sensor,sim);
elseif strcmp(sim.control_type,'NDI')
    control_law = @(t,X_sensor,Tc) NDI(t,X_sensor,sim);
elseif strcmp(sim.control_type,'NDI_tss')
    control_law = @(t,X_sensor,Tc) NDI_tss(t,X_sensor,sim);
elseif strcmp(sim.control_type,'INDI')
    control_law = @(t,X_sensor,Tc) INDI(t,X_sensor,Tc,sim);
end

% Run simulation loop

while clock_cycle <= final_cycle
    % Store state
    t = clock_cycle*sim.dt;
    raw_results(clock_cycle+1,:) = [t,X',Tc'];
    
    % Compute control torques
    if control_cycle == reset_counter
        Tc = control_law(t,X_sensor,Tc);
        X_sensor = X;
        control_cycle = 0;
    end
    
    % Propagate dynamics
    fun = @(X) spacecraft_dynamics(X,Tc,sim);
    X = RK4(X,fun,sim.dt);
    if strcmp(sim.att_repr,'Euler')
        X(1:3) = wrapToPi(X(1:3));
    end
    
    % Increment cycle
    clock_cycle = clock_cycle+1;
    control_cycle = control_cycle+1;
end

% Process results

ADCS_results.t = raw_results(:,1);
if strcmp(sim.att_repr,'Euler')
    ADCS_results.Euler = raw_results(:,2:4);
    quat = eul2quat(raw_results(:,2:4),'ZYX');
    ADCS_results.quaternions = [quat(:,2:4),quat(:,1)];
    ADCS_results.omega = raw_results(:,5:7);
    ADCS_results.Tc = raw_results(:,8:10);
    ADCS_results.theta = zeros(final_cycle+1,1);
    for i = 1:final_cycle+1
        [att_c,att_dot_c] = comanded_attitude(ADCS_results.t(i),sim);
        R_att = eul2rotm(ADCS_results.Euler(i,:),'ZYX');
        R_att_c = eul2rotm(att_c','ZYX');
        R_e = R_att_c/R_att;
        angle_axis = rotm2axang(R_e);
        ADCS_results.theta(i) = angle_axis(4);
    end
elseif strcmp(sim.att_repr,'quaternions')
    ADCS_results.quaternions = raw_results(:,2:5);
    ADCS_results.Euler = quat2eul([raw_results(:,5),raw_results(:,2:4)],'ZYX');
    ADCS_results.omega = raw_results(:,6:8);
    ADCS_results.Tc = raw_results(:,9:11);
    for i = 1:final_cycle+1
        [att_c,att_dot_c] = comanded_attitude(ADCS_results.t(i),sim);
        R_att = quat2rotm([ADCS_results.quaternions(i,4),ADCS_results.quaternions(i,1:3)]);
        R_att_c = quat2rotm([att_c(4)',att_c(1:3)']);
        R_e = R_att_c/R_att;
        angle_axis = rotm2axang(R_e);
        ADCS_results.theta(i) = angle_axis(4);
    end
end
ADCS_results.omega_norm = vecnorm(ADCS_results.omega,2,2);
ADCS_results.Tc_norm = vecnorm(ADCS_results.Tc,2,2);
end

function X = RK4(X,fun,dt)
k1 = fun(X);
k2 = fun(X+0.5*k1*dt);
k3 = fun(X+0.5*k2*dt);
k4 = fun(X+k3*dt);
X = X+1/6*dt*(k1+2*k2+2*k3+k4);
end