%%MATLABCONTROL - implicit solver (ode15s) for ODE system of Weiße et al.

%%INPUT
% init - initial conditions from bsim
% t0 - start point in time
% tf - end point in time

%%OUTPUT
% vector ynew
% because of java code issues the entries are given back instead of the
% whole vector

%%CHANGES
% changes to Weiße et al. ODE system: no dilution, scale parameters with
% M/Mref to get molecules per cell


%%driver
function [rmr,em,rmp,rmq,rmt,et,rmm,zmm,zmr,zmp,zmq,zmt,mt,mm,q,p,si,mq,mp,mr,r,a] = IMPLsolver(init, t0, tf)
    
% parameters
thetar= 426.8693338968694;
k_cm= 0.005990373118888;
gmax= 1260.0;
cl= 0;
nume0= 4.139172187824451;
s0= 1.0e4;
vm= 5800.0;
Km= 1.0e3;
numr0= 929.9678874564831;
nx= 300.0;
kq= 1.522190403737490e+05;
Kp= 180.1378030928276;
vt= 726.0;
nump0= 0.0;
numq0= 948.9349882947897;
Kt= 1.0e3;
nq= 4;
nr= 7549.0;
ns= 0.1;
thetax= 4.379733394834643;
parameters= [thetar k_cm gmax cl nume0 s0 vm Km numr0 nx kq Kp vt nump0 numq0 Kt nq nr ns thetax];

% define rate constants
b= 0;
dm= 0.1;
kb= 1;
ku= 1.0;
f= cl*k_cm;
rates= [b dm kb ku f];

% % define initial conditions
% rmr_0= 0;
% em_0= 0;
% rmp_0= 0;
% rmq_0= 0;
% rmt_0= 0;
% et_0= 0;
% rmm_0= 0;
% zmm_0= 0;
% zmr_0= 0;
% zmp_0= 0;
% zmq_0= 0;
% zmt_0= 0;
% mt_0= 0;
% mm_0= 0;
% q_0= 0;
% p_0= 0;
% si_0= 0;
% mq_0= 0;
% mp_0= 0;
% mr_0= 0;
% r_0= 10.0;
% a_0= 1000.0;
% init= [rmr_0 em_0 rmp_0 rmq_0 rmt_0 et_0 rmm_0 zmm_0 zmr_0 zmp_0 zmq_0 zmt_0 mt_0 mm_0 q_0 p_0 si_0 mq_0 mp_0 mr_0 r_0 a_0];
% 
% % call solver routine 
% %t0= 0;
% %tf= 2000;
[t,y]= ode15s(@(t,y) bsim_mainModel_odes(t, y, rates, parameters), [t0 tf], init);

    ynew = y(end,:);
    rmr= ynew(1);
    em= ynew(2);
    rmp= ynew(3);
    rmq= ynew(4);
    rmt= ynew(5);
    et= ynew(6);
    rmm= ynew(7);
    zmm= ynew(8);
    zmr= ynew(9);
    zmp= ynew(10);
    zmq= ynew(11);
    zmt= ynew(12);
    mt= ynew(13);
    mm= ynew(14);
    q= ynew(15);
    p= ynew(16);
    si= ynew(17);
    mq= ynew(18);
    mp= ynew(19);
    mr= ynew(20);
    r= ynew(21);
    a= ynew(22);
end


%%ODE-System
function dydt= bsim_mainModel_odes(t, y, rates, parameters)

	b= rates(1);
	dm= rates(2);
	kb= rates(3);
	ku= rates(4);
	f= rates(5);

	thetar= parameters(1);
	k_cm= parameters(2);
	gmax= parameters(3);
	cl= parameters(4);
	nume0= parameters(5);
	s0= parameters(6);
	vm= parameters(7);
	Km= parameters(8);
	numr0= parameters(9);
	nx= parameters(10);
	kq= parameters(11);
	Kp= parameters(12);
	vt= parameters(13);
	nump0= parameters(14);
	numq0= parameters(15);
	Kt= parameters(16);
	nq= parameters(17);
	nr= parameters(18);
	ns= parameters(19);
	thetax= parameters(20);

	rmr= y(1);
	em= y(2);
	rmp= y(3);
	rmq= y(4);
	rmt= y(5);
	et= y(6);
	rmm= y(7);
	zmm= y(8);
	zmr= y(9);
	zmp= y(10);
	zmq= y(11);
	zmt= y(12);
	mt= y(13);
	mm= y(14);
	q= y(15);
	p= y(16);
	si= y(17);
	mq= y(18);
	mp= y(19);
	mr= y(20);
	r= y(21);
	a= y(22);
    
    Mref = 1e8;
    M= nx*(q+et+em)+nr*(rmq+rmr+rmp+rmt+rmm);   %plus zm* if e.g. with GFP 
    
	Kg= gmax/Kp*M/Mref;
	gamma= gmax*a/(Kg+a);
	ttrate= (rmq+rmr+rmp+rmt+rmm)*gamma;
	lam= ttrate/M;
	nucat= em*vm*si/(Km*M/Mref+si);

	dydt(size(y,1),1)= 0;
	dydt(1)= +kb*Mref/M*r*mr+b*zmr-ku*rmr-gamma/nr*rmr-f*rmr;
	dydt(2)= +gamma/nx*rmm;
	dydt(3)= +kb*Mref/M*r*mp+b*zmp-ku*rmp-gamma/nx*rmp-f*rmp;
	dydt(4)= +kb*Mref/M*r*mq+b*zmq-ku*rmq-gamma/nx*rmq-f*rmq;
	dydt(5)= +kb*Mref/M*r*mt+b*zmt-ku*rmt-gamma/nx*rmt-f*rmt;
	dydt(6)= +gamma/nx*rmt;
	dydt(7)= +kb*Mref/M*r*mm+b*zmm-ku*rmm-gamma/nx*rmm-f*rmm;
	dydt(8)= +f*rmm-b*zmm;
	dydt(9)= +f*rmr-b*zmr;
	dydt(10)= +f*rmp-b*zmp;
	dydt(11)= +f*rmq-b*zmq;
	dydt(12)= +f*rmt-b*zmt;
	dydt(13)= +(nume0*M/Mref*a/(thetax*M/Mref+a))+ku*rmt+gamma/nx*rmt-kb*Mref/M*r*mt-dm*mt;
	dydt(14)= +(nume0*M/Mref*a/(thetax*M/Mref+a))+ku*rmm+gamma/nx*rmm-kb*Mref/M*r*mm-dm*mm;
	dydt(15)= +gamma/nx*rmq;
	dydt(16)= +gamma/nx*rmp;
	dydt(17)= +(et*vt*s0/(Kt+s0))-nucat; %why not Kt??
	dydt(18)= +(numq0*M/Mref*a/(thetax*M/Mref+a)/(1+(q/(kq*M/Mref))^nq))+ku*rmq+gamma/nx*rmq-kb*Mref/M*r*mq-dm*mq;
	dydt(19)= +(nump0*M/Mref*a/(thetax*M/Mref+a))+ku*rmp+gamma/nx*rmp-kb*Mref/M*r*mp-dm*mp;
	dydt(20)= +(numr0*M/Mref*a/(thetar*M/Mref+a))+ku*rmr+gamma/nr*rmr-kb*Mref/M*r*mr-dm*mr;
	dydt(21)= +ku*rmr+ku*rmt+ku*rmm+ku*rmp+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmp+gamma/nx*rmq-kb*Mref/M*r*mr-kb*Mref/M*r*mt-kb*Mref/M*r*mm-kb*Mref/M*r*mp-kb*Mref/M*r*mq;
	dydt(22)= +ns*nucat-ttrate;
end