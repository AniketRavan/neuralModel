function HH
    neurons = 2;
    g_Na = 120;  % mS/cm^2
    g_K = 36;  % mS/cm^2
    g_L = 0.3; % S/cm^2
    E_Na = 50 + 65;    % mV
    E_K = -77 + 65;    % mV
    E_L = -54.3 + 65;  % mV
    C = 1;    % mF/cm^2 
    firing_thresh = 20;
    dt = 0.01;
    TMAX = 70;
    t = 0:dt:TMAX;
    % Synapses
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ampa.tau1 = 1;
    ampa.tau2 = 6;
    ampa.E = 0;
    ampa.g_max = 1e-7;
    ampa.h = zeros(length(t), neurons);
    ampa.g = zeros(length(t), neurons);
    nmda.tau2 = 80;
    nmda.tau1 = 1;
    nmda.g_max = 6e-7;
    nmda.E = 60;
    nmda.h = zeros(length(t), neurons);
    nmda.g = zeros(length(t), neurons);
    glycine.tau1 = 1;
    glycine.tau2 = 2;
    glycine.E = -80;
    glycine.g_max = 1e-5;
    glycine.g = zeros(length(t), neurons);
    glycine.h = zeros(length(t), neurons);
    glycine.t_rise = glycine.tau1*glycine.tau2/(glycine.tau1-glycine.tau2);
    glycine.norm = ((glycine.tau2/glycine.tau1)^(glycine.t_rise/glycine.tau1) - (glycine.tau2/glycine.tau1)^(glycine.t_rise/glycine.tau2))^(-1);
    Vrest = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v = zeros(length(t), neurons);
    m = zeros(length(t), neurons);
    n = zeros(length(t), neurons);
    h = zeros(length(t), neurons);
    
    v(1, :) = Vrest;
    m(1, :) = a_m(Vrest)/(a_m(Vrest)+b_m(Vrest));
    h(1, :) = a_h(Vrest)/(a_h(Vrest)+b_h(Vrest));
    n(1, :) = a_n(Vrest)/(a_n(Vrest)+b_n(Vrest));
    connection(1) = 2;
    connection(2) = 1;
    for i = 2:length(t)
        for N = 1:neurons
            %I = zeros(1,length(t));
            %I(floor(length(t)/2):floor(length(t)/2) + 50) = 1;
            % update conductances
            if v(i - 1, connection(N)) > v(i, connection(N)) && v(i, connection(N)) > firing_thresh
                u = 1;
            else u = 0;
            end
            glycine.h(i, N) = glycine.h(i-1, N) + (- glycine.h(i-1, N)/glycine.tau1 + u)*dt;
            glycine.g(i, N) = glycine.g(i-1, N) + (- glycine.g(i-1)/glycine.tau2)*dt;
            I = 10*(t > 10 & t <= 60);
            gNa = g_Na * m(i - 1, N)^3 * h(i - 1, N);
            gK  = g_K  * n(i - 1, N)^4;
            mdot = a_m(v(i - 1, N))*(1-m(i - 1, N)) - b_m(v(i - 1, N))*m(i - 1, N);
            hdot = a_h(v(i - 1, N))*(1-h(i - 1, N)) - b_h(v(i - 1, N))*h(i - 1, N);
            ndot = a_n(v(i - 1, N))*(1-n(i - 1, N)) - b_n(v(i - 1, N))*n(i - 1, N);
            INa = gNa*(v(i - 1, N)-E_Na);
            IK = gK*(v(i - 1, N)-E_K);
            IL = g_L*(v(i - 1, N)-E_L);
            Iglycine = glycine.g_max*glycine.norm*glycine.g(i, N)*(v(i - 1, N)-glycine.E);
            vdot = (I(i-1) - INa - IK - IL - Iglycine)/C;
            m(i) = m(i-1) + mdot*dt;		
            h(i) = h(i-1) + hdot*dt;
            n(i) = n(i-1) + ndot*dt;
            v(i) = v(i-1) + vdot*dt;
        end
    
        figure, plot(t,v);
    end
end

function z = a_n(V)
    z = 0.01*(10 - V)./(exp((10 - V)/10) - 1);
    idx = find(10 - V == 0);
    z(idx) = 0.1;
end

function z = b_n(V)
    z = 0.125*exp(-V/80);
end

function z = a_m(V)
    z = 0.1*(25 - V)./(exp((25 - V)/10) - 1);
    idx = find(25 - V == 0);
    z(idx) = 1;
end

function z = b_m(V)
    z = 4*exp(-V/18);
end

function z = a_h(V)
    z = 0.07*exp(-V/20);
end

function z = b_h(V)
    z = 1./(exp((30 -V)/10) + 1);
end

