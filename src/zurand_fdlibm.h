/* exp() and log1p() for the ziggurat, from fdlibm 5.3 (netlib.org/fdlibm,
 * e_exp.c 1.6 04/04/22 and s_log1p.c 1.3 95/01/18).
 *
 * Why zurand carries its own: the platform's libm is not the same function
 * everywhere. Measured 2026-09-26 on the whole-stream digests: Apple's
 * libm and glibc/musl disagree often enough, over the ~4,300 calls in 3e6
 * normals, to change the xoshiro and philox streams, and neither is
 * correctly rounded (Apple's log1p is off by one ulp on 6% of those
 * calls). fdlibm is a fixed sequence of IEEE double operations, so it
 * gives the same bits on every platform that evaluates them as written --
 * which is also why Java's StrictMath is fdlibm.
 *
 * DERIVED, NOT VERBATIM. Unlike src/Random123 and src/numpyzig, this file
 * is edited, because the originals are not safe to compile as they are:
 *
 *  - every product that feeds an addition or subtraction goes through
 *    zurand_rounded(), because a fused multiply-add changes the result and
 *    GCC fuses by default on arm64 (see the note above zurand_rounded);
 *  - the __HI/__LO pointer casts, undefined behaviour under strict
 *    aliasing, are replaced by memcpy;
 *  - exp's `__HI(y) += k << 20` shifted a negative int, also undefined, and
 *    is done in unsigned arithmetic with the same result;
 *  - K&R declarations and the public names are gone; both are static.
 *
 * The algorithms, constants and branch structure are otherwise unchanged,
 * as are the comments below the notices. Accuracy: < 1 ulp (fdlibm's
 * error analysis), which is irrelevant statistically; what matters here is
 * that the result is the same everywhere.
 *
 * Include from zurand.c only, after zurand_rounded() is defined.
 */

/* ---- word access (replaces fdlibm.h's __HI / __LO) ---- */

static inline int32_t fd_hi(double x) {
    uint64_t b;
    memcpy(&b, &x, sizeof b);
    return (int32_t)(uint32_t)(b >> 32);
}

static inline uint32_t fd_lo(double x) {
    uint64_t b;
    memcpy(&b, &x, sizeof b);
    return (uint32_t)b;
}

static inline double fd_with_hi(double x, uint32_t hi) {
    uint64_t b;
    memcpy(&b, &x, sizeof b);
    b = ((uint64_t)hi << 32) | (b & UINT64_C(0xffffffff));
    memcpy(&x, &b, sizeof x);
    return x;
}

#define FD_R(p) zurand_rounded(p)

/* ====================================================
 * Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

/* __ieee754_exp(x)
 * Returns the exponential of x.
 *
 * Method
 *   1. Argument reduction:
 *      Reduce x to an r so that |r| <= 0.5*ln2 ~ 0.34658.
 *	Given x, find r and integer k such that
 *
 *               x = k*ln2 + r,  |r| <= 0.5*ln2.
 *
 *      Here r will be represented as r = hi-lo for better
 *	accuracy.
 *
 *   2. Approximation of exp(r) by a special rational function on
 *	the interval [0,0.34658]:
 *	Write
 *	    R(r**2) = r*(exp(r)+1)/(exp(r)-1) = 2 + r*r/6 - r**4/360 + ...
 *      We use a special Remes algorithm on [0,0.34658] to generate
 * 	a polynomial of degree 5 to approximate R. The maximum error
 *	of this polynomial approximation is bounded by 2**-59. In
 *	other words,
 *	    R(z) ~ 2.0 + P1*z + P2*z**2 + P3*z**3 + P4*z**4 + P5*z**5
 *  	(where z=r*r, and the values of P1 to P5 are listed below)
 *	and
 *	    |                  5          |     -59
 *	    | 2.0+P1*z+...+P5*z   -  R(z) | <= 2
 *	    |                             |
 *	The computation of exp(r) thus becomes
 *                             2*r
 *		exp(r) = 1 + -------
 *		              R - r
 *                                 r*R1(r)
 *		       = 1 + r + ----------- (for better accuracy)
 *		                  2 - R1(r)
 *	where
 *			         2       4             10
 *		R1(r) = r - (P1*r  + P2*r  + ... + P5*r   ).
 *
 *   3. Scale back to obtain exp(x):
 *	From step 1, we have
 *	   exp(x) = 2^k * exp(r)
 *
 * Special cases:
 *	exp(INF) is INF, exp(NaN) is NaN;
 *	exp(-INF) is 0, and
 *	for finite argument, only exp(0)=1 is exact.
 *
 * Accuracy:
 *	according to an error analysis, the error is always less than
 *	1 ulp (unit in the last place).
 *
 * Misc. info.
 *	For IEEE double
 *	    if x >  7.09782712893383973096e+02 then exp(x) overflow
 *	    if x < -7.45133219101941108420e+02 then exp(x) underflow
 *
 * Constants:
 * The hexadecimal values are the intended ones for the following
 * constants. The decimal values may be used, provided that the
 * compiler will convert from decimal to binary accurately enough
 * to produce the hexadecimal values shown.
 */

static const double
fd_one	= 1.0,
fd_halF[2]	= {0.5,-0.5,},
fd_huge	= 1.0e+300,
fd_twom1000= 9.33263618503218878990e-302,     /* 2**-1000=0x01700000,0*/
fd_o_threshold=  7.09782712893383973096e+02,  /* 0x40862E42, 0xFEFA39EF */
fd_u_threshold= -7.45133219101941108420e+02,  /* 0xc0874910, 0xD52D3051 */
fd_ln2HI[2]   ={ 6.93147180369123816490e-01,  /* 0x3fe62e42, 0xfee00000 */
	     -6.93147180369123816490e-01,},/* 0xbfe62e42, 0xfee00000 */
fd_ln2LO[2]   ={ 1.90821492927058770002e-10,  /* 0x3dea39ef, 0x35793c76 */
	     -1.90821492927058770002e-10,},/* 0xbdea39ef, 0x35793c76 */
fd_invln2 =  1.44269504088896338700e+00, /* 0x3ff71547, 0x652b82fe */
fd_P1   =  1.66666666666666019037e-01, /* 0x3FC55555, 0x5555553E */
fd_P2   = -2.77777777770155933842e-03, /* 0xBF66C16C, 0x16BEBD93 */
fd_P3   =  6.61375632143793436117e-05, /* 0x3F11566A, 0xAF25DE2C */
fd_P4   = -1.65339022054652515390e-06, /* 0xBEBBBD41, 0xC5D26BF1 */
fd_P5   =  4.13813679705723846039e-08; /* 0x3E663769, 0x72BEA4D0 */

static double zurand_exp(double x)
{
	double y,hi,lo,c,t;
	int k,xsb;
	uint32_t hx;

	hi = lo = 0.0;
	k = 0;
	hx  = (uint32_t)fd_hi(x);	/* high word of x */
	xsb = (int)((hx>>31)&1);	/* sign bit of x */
	hx &= 0x7fffffff;		/* high word of |x| */

    /* filter out non-finite argument */
	if(hx >= 0x40862E42) {			/* if |x|>=709.78... */
            if(hx>=0x7ff00000) {
		if(((hx&0xfffff)|fd_lo(x))!=0)
		     return x+x; 		/* NaN */
		else return (xsb==0)? x:0.0;	/* exp(+-inf)={inf,0} */
	    }
	    if(x > fd_o_threshold) return fd_huge*fd_huge; /* overflow */
	    if(x < fd_u_threshold) return fd_twom1000*fd_twom1000; /* underflow */
	}

    /* argument reduction */
	if(hx > 0x3fd62e42) {		/* if  |x| > 0.5 ln2 */
	    if(hx < 0x3FF0A2B2) {	/* and |x| < 1.5 ln2 */
		hi = x-fd_ln2HI[xsb]; lo=fd_ln2LO[xsb]; k = 1-xsb-xsb;
	    } else {
		k  = (int)(FD_R(fd_invln2*x)+fd_halF[xsb]);
		t  = k;
		hi = x - FD_R(t*fd_ln2HI[0]);	/* t*ln2HI is exact here */
		lo = FD_R(t*fd_ln2LO[0]);
	    }
	    x  = hi - lo;
	}
	else if(hx < 0x3e300000)  {	/* when |x|<2**-28 */
	    if(fd_huge+x>fd_one) return fd_one+x;/* trigger inexact */
	}
	else k = 0;

    /* x is now in primary range */
	t  = x*x;
	c  = x - FD_R(t*(fd_P1+FD_R(t*(fd_P2+FD_R(t*(fd_P3+FD_R(t*(fd_P4+FD_R(t*fd_P5)))))))));
	if(k==0) 	return fd_one-((x*c)/(c-2.0)-x);
	else 		y = fd_one-((lo-(x*c)/(2.0-c))-hi);
	if(k >= -1021) {
	    /* add k to y's exponent */
	    return fd_with_hi(y, (uint32_t)fd_hi(y) + ((uint32_t)k << 20));
	} else {
	    /* add k to y's exponent */
	    y = fd_with_hi(y, (uint32_t)fd_hi(y) + ((uint32_t)(k+1000) << 20));
	    return y*fd_twom1000;
	}
}

/* ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

/* double log1p(double x)
 *
 * Method :
 *   1. Argument Reduction: find k and f such that
 *			1+x = 2^k * (1+f),
 *	   where  sqrt(2)/2 < 1+f < sqrt(2) .
 *
 *      Note. If k=0, then f=x is exact. However, if k!=0, then f
 *	may not be representable exactly. In that case, a correction
 *	term is need. Let u=1+x rounded. Let c = (1+x)-u, then
 *	log(1+x) - log(u) ~ c/u. Thus, we proceed to compute log(u),
 *	and add back the correction term c/u.
 *	(Note: when x > 2**53, one can simply return log(x))
 *
 *   2. Approximation of log1p(f).
 *	Let s = f/(2+f) ; based on log(1+f) = log(1+s) - log(1-s)
 *		 = 2s + 2/3 s**3 + 2/5 s**5 + .....,
 *	     	 = 2s + s*R
 *      We use a special Reme algorithm on [0,0.1716] to generate
 * 	a polynomial of degree 14 to approximate R The maximum error
 *	of this polynomial approximation is bounded by 2**-58.45. In
 *	other words,
 *		        2      4      6      8      10      12      14
 *	    R(z) ~ Lp1*s +Lp2*s +Lp3*s +Lp4*s +Lp5*s  +Lp6*s  +Lp7*s
 *  	(the values of Lp1 to Lp7 are listed in the program)
 *	and
 *	    |      2          14          |     -58.45
 *	    | Lp1*s +...+Lp7*s    -  R(z) | <= 2
 *	    |                             |
 *	Note that 2s = f - s*f = f - hfsq + s*hfsq, where hfsq = f*f/2.
 *	In order to guarantee error in log below 1ulp, we compute log
 *	by
 *		log1p(f) = f - (hfsq - s*(hfsq+R)).
 *
 *	3. Finally, log1p(x) = k*ln2 + log1p(f).
 *		 	     = k*ln2_hi+(f-(hfsq-(s*(hfsq+R)+k*ln2_lo)))
 *	   Here ln2 is split into two floating point number:
 *			ln2_hi + ln2_lo,
 *	   where n*ln2_hi is always exact for |n| < 2000.
 *
 * Special cases:
 *	log1p(x) is NaN with signal if x < -1 (including -INF) ;
 *	log1p(+INF) is +INF; log1p(-1) is -INF with signal;
 *	log1p(NaN) is that NaN with no signal.
 *
 * Accuracy:
 *	according to an error analysis, the error is always less than
 *	1 ulp (unit in the last place).
 *
 * Constants:
 * The hexadecimal values are the intended ones for the following
 * constants. The decimal values may be used, provided that the
 * compiler will convert from decimal to binary accurately enough
 * to produce the hexadecimal values shown.
 */

static const double
fd_ln2_hi  =  6.93147180369123816490e-01,	/* 3fe62e42 fee00000 */
fd_ln2_lo  =  1.90821492927058770002e-10,	/* 3dea39ef 35793c76 */
fd_two54   =  1.80143985094819840000e+16,  /* 43500000 00000000 */
fd_Lp1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
fd_Lp2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
fd_Lp3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
fd_Lp4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
fd_Lp5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
fd_Lp6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
fd_Lp7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */

static double fd_zero = 0.0;

static double zurand_log1p(double x)
{
	double hfsq,f,c,s,z,R,u;
	int k,hx,hu,ax;

	f = c = 0.0;
	hu = 0;
	hx = fd_hi(x);		/* high word of x */
	ax = hx&0x7fffffff;

	k = 1;
	if (hx < 0x3FDA827A) {			/* x < 0.41422  */
	    if(ax>=0x3ff00000) {		/* x <= -1.0 */
		if(x==-1.0) return -fd_two54/fd_zero; /* log1p(-1)=+inf */
		else return (x-x)/(x-x);	/* log1p(x<-1)=NaN */
	    }
	    if(ax<0x3e200000) {			/* |x| < 2**-29 */
		if(fd_two54+x>fd_zero			/* raise inexact */
	            &&ax<0x3c900000) 		/* |x| < 2**-54 */
		    return x;
		else
		    return x - FD_R(x*x*0.5);
	    }
	    if(hx>0||hx<=((int)0xbfd2bec3)) {
		k=0;f=x;hu=1;}	/* -0.2929<x<0.41422 */
	}
	if (hx >= 0x7ff00000) return x+x;
	if(k!=0) {
	    if(hx<0x43400000) {
		u  = 1.0+x;
	        hu = fd_hi(u);		/* high word of u */
	        k  = (hu>>20)-1023;
	        c  = (k>0)? 1.0-(u-x):x-(u-1.0);/* correction term */
		c /= u;
	    } else {
		u  = x;
	        hu = fd_hi(u);		/* high word of u */
	        k  = (hu>>20)-1023;
		c  = 0;
	    }
	    hu &= 0x000fffff;
	    if(hu<0x6a09e) {
	        u = fd_with_hi(u, (uint32_t)(hu|0x3ff00000));	/* normalize u */
	    } else {
	        k += 1;
	        u = fd_with_hi(u, (uint32_t)(hu|0x3fe00000));	/* normalize u/2 */
	        hu = (0x00100000-hu)>>2;
	    }
	    f = u-1.0;
	}
	hfsq=FD_R(0.5*f*f);
	if(hu==0) {	/* |f| < 2**-20 */
	    if(f==fd_zero) { if(k==0) return fd_zero;
			else {c += FD_R(k*fd_ln2_lo); return FD_R(k*fd_ln2_hi)+c;} }
	    R = FD_R(hfsq*(1.0-FD_R(0.66666666666666666*f)));
	    if(k==0) return f-R; else
	    	     return FD_R(k*fd_ln2_hi)-((R-(FD_R(k*fd_ln2_lo)+c))-f);
	}
 	s = f/(2.0+f);
	z = s*s;
	/* R = z*(Lp1+z*(Lp2+z*(Lp3+z*(Lp4+z*(Lp5+z*(Lp6+z*Lp7)))))),
	 * written out step by step: the same operations in the same order. */
	R = fd_Lp6+FD_R(z*fd_Lp7);
	R = fd_Lp5+FD_R(z*R);
	R = fd_Lp4+FD_R(z*R);
	R = fd_Lp3+FD_R(z*R);
	R = fd_Lp2+FD_R(z*R);
	R = fd_Lp1+FD_R(z*R);
	R = FD_R(z*R);
	if(k==0) return f-(hfsq-FD_R(s*(hfsq+R))); else
		 return FD_R(k*fd_ln2_hi)-((hfsq-(FD_R(s*(hfsq+R))+(FD_R(k*fd_ln2_lo)+c)))-f);
}

#undef FD_R
