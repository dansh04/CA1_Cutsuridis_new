: from Dave et al (2021), Hippocampal CA1 Area

TITLE nmda synapse 

NEURON {
	POINT_PROCESS NMDA2
	NONSPECIFIC_CURRENT i
    RANGE gmax, enmda, eca, onset, tauF, tauS, mg, pnmda, pca, alpha, beta, T, v, gnmda, i, B
}

UNITS {
    (nS) = (microsiemens)
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (mM) = (milli/liter)
}

PARAMETER {
	gmax = 0.00065 (uS)
	enmda = 0      (mV)
	eca = 130      (mV)
	onset = 10     (ms)
   tauF = 9       (ms)
   tauS = 90      (ms)
   mg = 1         (mM)
   alpha = 0.072  <0,1e4>         
	beta = 0.0066
	T = 1
	util = 0.3
}

ASSIGNED {
	v 	    (mV)
	gnmda   (uS)
	i       (nA)
	B       : fraction of NMDA receptors not affected by magnesium block
}

STATE {
    r       : fraction of open channels
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit 
    B = mgblock(v)
    gnmda = gmax*((tauF*tauS)/(tauF-tauS))*duale((t-onset)/tauF,(t-onset)/tauS)
    i = gnmda*B*r*(v-enmda)
}

DERIVATIVE state {
    r' = alpha*T*(1-r)-beta*r
}

INITIAL {
    r = 0
}

FUNCTION mgblock(v(mV)) {
    mgblock = 1/(1+exp(-0.062(/mV)*v)*(mg/3.57(mM)))
}

FUNCTION duale(x,y) {
	if (x<0 || y<0) {
		duale = 0
	} 
	else {
		duale = exp(-x)-exp(-y)
	}
}

NET_RECEIVE(wgt,R,u,tlast (ms),nspike) {
    LOCAL x
    :printf("entry flag=%g t=%g\n", flag, tlast)
    COMMENT
    if (nspike==0) { 
        R=1  
        u=util 
    }
	else {
	     if (tauF>0) { u=util+(1-util)*u*exp(-(t-tlast)/tauF) }
	     if (tauD>0) { R=1+(R*(1-u)-1)*exp(-(t-tlast)/tauD) }
	}
	ENDCOMMENT
	x = wgt
	state_discontinuity(r,r+x)
    tlast = t
    nspike = nspike+1
}
