function e = step(e,w)
    A = [1 1;0 1];
    B = [0; 1];
    K = [-0.4 -1.2];
    e = (A+B*K)*e+w;
    e.plot()