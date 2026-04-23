function [picos_ddx_m, picos_x_m] = ResponseSpectrum(f_vector_accel, ref_accel, ref_disp, f_vector_disp)
    m    = 1;
    zeta = 0.05;

    accel_points = length(f_vector_accel);
    picos_ddx_m  = zeros(accel_points, 1);
    picos_x_m    = [];

    do_displacement = nargin == 5;  % both ref_disp and f_vector_disp provided

    if do_displacement
        picos_x_m = zeros(length(f_vector_disp), 1);
    end

    Ts=0.005;
    for i = 1:accel_points
        k = m * (2*pi*f_vector_accel(i))^2;
        c = zeta * 2 * m * 2*pi*f_vector_accel(i);

        sys_d  = c2d(tf([c k], [m c k]), Ts, 'zoh');
        [b, a] = tfdata(sys_d, 'v');

        picos_ddx_m(i) = max(abs(filter(b, a, ref_accel)));

        if do_displacement
            [is_in, idx] = ismember(f_vector_accel(i), f_vector_disp);
            if is_in
                picos_x_m(idx) = max(abs(filter(b, a, ref_disp)));
            end
        end
    end
end