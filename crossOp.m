function ans = crossOp(v);

v1 = v(1); v2 = v(2); v3 = v(3);

ans = [0 -v3 v2;...
    v3 0 -v1;...
    -v2 v1 0];