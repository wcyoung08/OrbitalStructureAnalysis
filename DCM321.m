function T=DCM321(ang)

ps = ang(1); th = ang(2); ph = ang(3);

T = [cos(th)*cos(ps)                            cos(th)*sin(ps)                             -sin(th);...
    sin(ph)*sin(th)*cos(ps)-cos(ph)*sin(ps)     sin(ph)*sin(th)*sin(ps)+cos(ph)*cos(ps)     sin(ph)*cos(th);...
    cos(ph)*sin(th)*cos(ps)+sin(ph)*sin(ps)     cos(ph)*sin(th)*sin(ps)-sin(ph)*cos(ps)     cos(ph)*cos(th)];

