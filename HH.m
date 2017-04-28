function HH
    g_Na = 1200;  % S/cm^2
    g_K = 360;  % S/cm^2
    g_L = 3; % S/cm^2
    E_Na = 50 + 65;    % mV
    E_K = -77 + 65;    % mV
    E_L = -54.3 + 65;  % mV
    C = 1;    % mF/cm^2 
    
    dt = 0.001;
    TMAX = 70;
    t = 0:dt:TMAX;
    
    Vrest = 0;
    
    v = zeros(size(t));
    m = zeros(size(t));
    n = zeros(size(t));
    h = zeros(size(t));
    
    v(1) = Vrest;						
    m(1) = a_m(Vrest)/(a_m(Vrest)+b_m(Vrest));
    h(1) = a_h(Vrest)/(a_h(Vrest)+b_h(Vrest));
    n(1) = a_n(Vrest)/(a_n(Vrest)+b_n(Vrest));
    for i = 2:length(t)
        %I = zeros(1,length(t));
        %I(floor(length(t)/2):floor(length(t)/2) + 50) = 1;
        I = 10*(t > 10 & t <= 40);
        gNa = g_Na * m(i - 1)^3 * h(i - 1);
        gK  = g_K  * n(i - 1)^4;
        mdot = a_m(v(i - 1))*(1-m(i - 1)) - b_m(v(i - 1))*m(i - 1);
        hdot = a_h(v(i - 1))*(1-h(i - 1)) - b_h(v(i - 1))*h(i - 1);
        ndot = a_n(v(i - 1))*(1-n(i - 1)) - b_n(v(i - 1))*n(i - 1);
        vdot = (I(i-1) - gNa*(v(i - 1)-E_Na) - gK*(v(i - 1)-E_K) - g_L*(v(i - 1)-E_L))/C;
        m(i) = m(i-1) + mdot*dt;		
        h(i) = h(i-1) + hdot*dt;
        n(i) = n(i-1) + ndot*dt;
        v(i) = v(i-1) + vdot*dt;
    end

    figure, plot(t,v);
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

