function T=DCM123(ang)

ps = ang(1); th = ang(2); ph = ang(3);

T = [cos(th)*cos(ps)     sin(ph)*sin(th)*cos(ps)+cos(ph)*sin(ps)     -cos(ph)*sin(th)*cos(ps)+sin(ph)*sin(ps);...
    -cos(th)*sin(ps)     -sin(ph)*sin(th)*sin(ps)+cos(ph)*cos(ps)    cos(ph)*sin(th)*sin(ps)+sin(ps)*cos(ps);...
    sin(th)              -sin(ph)*cos(ph)     cos(ph)*cos(th)];

