


// void solve_sphere (double A, double b, double cc, double sc, double *cap, double *Bp)
// {
// 	double cb = cos(b), sb = sin(b);
// 	double sA, cA = cos(A);
// 	double x, y;
// 	double ca;
// 	double B;

// 	ca = cb*cc + sb*sc*cA;
// 	if (ca >  1.0) ca =  1.0;
// 	if (ca < -1.0) ca = -1.0;
// 	if (cap)
// 	    *cap = ca;

// 	if (!Bp)
// 	    return;

// 	if (sc < 1e-7){
// 	    B = cc < 0 ? A : PI-A;
// 	}
// 	else {
// 	    sA = sin(A);
// 	    y = sA*sb*sc;
// 	    x = cb - ca*cc;
// 	    B = y ? (x ? atan2(y,x) : (y>0 ? PI/2 : -PI/2)) : (x>=0 ? 0 : PI);
// 	}

// 	*Bp = B;
// 	range (Bp, 2*PI);
// }
