function dm = valve_flow(valve_diameter, lift, P1, P2, T1, Cd, ...
        gamma, dm)
R = 0.287 + 0.02;

    A = lift * 2 * pi * valve_diameter;  % valve flow area

    rho = P1 / (R * T1);

    u = dm / (rho * A * Cd);

    P0 = P1 * (1 + ((gamma - 1) / 2) * (u / sqrt(gamma * R * T1)) ^ 2)^...
    (gamma / (gamma -1  ));   % stagnation pressure upstream of the valve

    critical_pressure_ratio = (2 / (gamma + 1) ) ^ (-gamma / (gamma - 1)); % pressure ratio at which the flow will choke

    pressure_ratio = P0 / P2;  % pressure ratio between upstream of the valve and downstream

    if pressure_ratio >= critical_pressure_ratio

        dm = Cd * A * P1 * sqrt((gamma / (R * T1)) * (2 / (gamma +1)) ^...
        ((gamma +1) / (gamma -1)));

    else

        dm = Cd * A * P1 * sqrt(((2 * gamma) / (gamma -1)) * (1 / (R * T1) ) * ...
            ((P2 / P1)^(2/gamma) - (P2 / P1)^((gamma + 1)/gamma)));

    end
end