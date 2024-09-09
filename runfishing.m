
function [ftot2, yieldtot2] = runfishing(slope)

    % Empty vectors
    yieldtot2 = zeros(999, 1);
    ftot2 = zeros(999, 1);

    for z = 1:999
        yhis = zeros(6, 100);
        yhis(:, 1) = [3*10^12; 3*10^9; 3*10^8; 3*10^8; 3*10^7; 3*10^7]; % Initial conditions

        for p = 1:99
            S2 = yhis(3, p) + yhis(4, p) + yhis(5, p) + yhis(6, p);
            T2 = @(t22) 15 + t22 * slope;
            T2_value = T2(p);
            Teff2 = @(T2) max(0, 1/5 * heaviside(-T2 + 17.219) .* (exp(0.2 * (T2 - 8)) - 1) + ...
                                1/5 * heaviside(T2 - 17.219) .* 4 .* log(-(T2 - 21)));
            Teff2_value = Teff2(T2_value);

            f22 = (6000 * exp(-0.0000000003 * S2)) * Teff2_value;
            f32 = (8000 * exp(-0.0000000003 * S2)) * Teff2_value;
            f42 = (10000 * exp(-0.0000000003 * S2)) * Teff2_value;
            f52 = (12000 * exp(-0.0000000003 * S2)) * Teff2_value;

            m232 = 0.6;
            m342 = 0.3;
            m452 = 0.4;
            m552 = 0.8;

            F232 = 0.0003 * z;
            F342 = 0.000325 * z;
            F452 = 0.00035 * z;
            F552 = 0.000375 * z;

            P232 = exp(-m232 - F232);
            P342 = exp(-m342 - F342);
            P452 = exp(-m452 - F452);
            P552 = exp(-m552 - F552);

            P012 = 0.001;
            P122 = 0.1;

            L2 = [0 0 f22 f32 f42 f52;
                  P012 0 0 0 0 0;
                  0 P122 0 0 0 0;
                  0 0 P232 0 0 0;
                  0 0 0 P342 0 0;
                  0 0 0 0 P452 P552];

            yhis(:, p+1) = L2 * yhis(:, p); % Update state
        end

        % Calculate yields
        yield12 = yhis(3, end) * (1 - P232) * F232 / (F232 + m232);
        yield22 = yhis(4, end) * (1 - P342) * F342 / (F342 + m342);
        yield32 = yhis(5, end) * (1 - P452) * F452 / (F452 + m452);
        yield42 = yhis(6, end) * (1 - P552) * F552 / (F552 + m552);
        yieldtot2(z) = yield12 + yield22 + yield32 + yield42;

        ftot2(z) = (F232 + F342 + F452 + F552) / 4;
    end

end



