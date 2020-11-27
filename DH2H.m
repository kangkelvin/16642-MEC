function H = DH2H(theta, d, a, alpha)
    Rz_theta = [cos(theta), -sin(theta) 0, 0;
                sin(theta), cos(theta), 0, 0;
                0, 0, 1, 0;
                0, 0, 0, 1];
    tz_d = [1, 0, 0, 0;
            0, 1, 0, 0;
            0, 0, 1, d;
            0, 0, 0, 1];
    tz_a = [1, 0, 0, a;
            0, 1, 0, 0;
            0, 0, 1, 0;
            0, 0, 0, 1];
    Rx_alpha = [1, 0, 0, 0;
                0, cos(alpha), -sin(alpha), 0;
                0, sin(alpha), cos(alpha), 0;
                0, 0, 0, 1];
    H = Rz_theta * tz_d * tz_a * Rx_alpha;
end