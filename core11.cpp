#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "StdAfx.h"

// Konstants ****************************************************
    #define BurgersVec (0.286e-09)
	#define KBOLTZMANN (1.380658e-23)
	#define ELECTRONCHARGE (1.60217733e-19)
	#define _PI_ (3.1415926536)

// Inizialisation -- globale Veriablen **************************

	simData val;
	curveData exper;
	allRho rho;
	allRho rhoDot;
	fitParmSet fit;

// *************************

// *************************

	Vector n[24];
	Vector b[24];

	Vector n1[24];
	Vector b1[24];

	Vector hkl;
	Vector hkl_b;
	Vector Basis[3];
	Vector Z;
////////////////////////////////////////////////////////////
	Matrix g_operator;
	Matrix r_operator[24];

	Matrix StressTensor;
	Matrix StrainTensor;
////////////////////////////////////////////////////////////
Matrix ShearTensor[24];
Matrix RotationTensor[24];

double StressInSlipSystem[24];
double ShearInSlipSystem[24];


// **************************************************************************************************
	long i;
	long j;
	double dt;
	double SigmaShear;

	long k1;
	long k2;
	double tauwt, vG;
	double Norm;



// function deklaration - OPERATOR

void GenerateSterssTensor (void);	// procedure
void GenerateStrainTensor (void);	// procedure

void GenerateSetOfSlipSystems (void);	// procedure 
void GenerateG_operator (void);		  // procedure	

 double ScalarProduct(Vector a, Vector b); // Function
 Matrix  TensorProduct(Vector a, Vector b);
 Matrix _TensorProduct(Vector a, Vector b);

double ShearStress (Vector m, Vector n, Matrix Sigma);
Vector MatrixVector (Matrix mat, Vector vec);

Vector MatrixVectorReverse (Matrix mat, Vector vec);

// funktion deklaration - CORE 

	void prepareCurveConstants( void );
	void prepareIntegrationConstants( void );
	long diffEquations( void );
	long integrateRho( double timeStep );
	double v_SS(double StressSlipSystem);

// MAIN *****************************************************************

int main(void)
{

/////////////////// Inizialisation /////////////////////////
// ***************************************************** ///

// Inizilisation fit (fitParmSet fit) - model parameters
	
	fit.dAnn 	= 1.000000e+00 * BurgersVec; // [m]
	fit.dLock 	= 2.042741e+00 * BurgersVec; // [m]
	fit.dImmob 	= 1.131503e+00 * BurgersVec; // [m]
	fit.dClear 	= 3.211728e+00 * BurgersVec; // [m]
	fit.dClimb 	= 1.756252e+00 * BurgersVec; // [m]

	fit.fWall 	= 1.043873e-01;

	fit.betaLeffIntern = 1.254158e+02;
	fit.betaLeffWall   = 9.784541e+01;
	
	fit._betaLeffIntern = 1 /1.254158e+02;
	fit._betaLeffWall   = 1 /9.784541e+01;
	
	fit.alphaTaylor 	= 5.000000e-01;

	fit.d0solute 	= 1.000000e+00 * BurgersVec; // [m]
	fit.Vcross 	= 9.942005e+01 * BurgersVec * BurgersVec * BurgersVec; // [m^3]

	fit.Qcross 	= 1.363390e+00 * ELECTRONCHARGE;  // [J]
	fit.Qclimb 	= 9.992049e+00 * ELECTRONCHARGE;  // [J]
	
	fit.rhoMobil0 	= 8.282318e+10;  //8.282318e+10; [m^-2]
	fit.rhoIntern0 	= 2.455309e+09;  //2.455309e+09; [m^-2]
	fit.rhoWall0 	= 5.778842e+13;  //5.778842e+13; [m^-2]
	
	fit.thermSolStrength1 = 1.778076e+02 * 1.0e+06; // [MPa]
	fit.thermSolStrength2 = 9.999999e-01 * 1.0e+06; // [MPa]
	
	fit.athermSolStrength1 = 2.524088e+00 * 1.0e+06;// [MPa] 
	fit.athermSolStrength2 = 9.999999e-01 * 1.0e+06;// [MPa]

	fit.particleStrength = 1.000000e+01 * 1.0e+06;  // [MPa]
	fit.numActivePlanes  = 1.000000e+00;

// Inizialisation val (simData val) - work Variablen *** *** *** *** 

	val.burgersVec = 0.286000e-09;  // [m]	// real constants
	val._burgersVec = 0;             // [1/m]
	val.shearMod0  = 27.00000e+09;  // [Pa]
	val.T0G_Kelvin = 300.0;         // [K]
	val.dG_dT      = -0.05e+09;     // [Pa/K]
	val.Lgrain0    = 0.001;	       // [m]
	val._Lgrain0   = 0;
	val.diffCoeff0 = 1.3e-04;       // [m^2/s]
	
	val.shearMod  = 0;	// T-, epsp-, M- dependent
	val.kT        = 0;
	val._kT       = 0;
	val.Gb_4p     = 0;
	val.ttFac     = 0;
	val.MepsDot_b = 0;

	val.fWall  =  0;		// fitParm-dependent
	val._fWall  = 0;
	val.fIntern = 0;
	
	val.tauSolTherm0 = 0;	// chemistry-dependent
	val.tauSolAtherm = 0;

	val.lambdaSol    = 0;
	val._lambdaSol   = 0;
	val.tauParticle  = 0;
	
	val._DGsol0  = 0;
	val.annFac   = 0;
	val.lockFac  = 0;
	val.dipolFac = 0;
	val.immobFac = 0;
	val.clearFac = 0;
	val.climbFac = 0;

//	val.numActivePlanes = 4.000000e+00;
//	val._numActivePlanes = 0;
	val.simTau = 0;

// Inizialisation Variablen rho (rho_)

	rho.mobil = 0;
	rho.intern = 0;
	rho.wall = 0;

// Inizialisaion Variablen rhoDot (rhoDot_)

	rhoDot.mobil = 0;
	rhoDot.intern = 0;
	rhoDot.wall = 0;


// Inizialisation exper (experiment exper)*************

	exper.expSigma = 0;
	exper.simSigma = 0;

//	exper.MTaylor  = 3.0;
	exper.T_Kelvin = 300;

	exper.strain = 0;
	exper.time   = 0;
	exper.strainRate = 1;

	exper.rhoM = 0;
	exper.rhoW = 0;
	exper.rhoI = 0;


	exper.Lgrain0 = 0.001;
	exper.atConc1 = 0.01;
	exper.atConc2 = 0.000000001;
	exper.particleVolumeFraction = 0.00000000001;
	exper.particleRadius_b = 0.00000000001;


// ***************************************************** ///
////////////////////////////////////////////////////////////
//////////  Reinizialisation  //////////////////////////////

	exper.T_Kelvin = 300;
	double MachineStrainRate = 1;

	double Epsilon = 0;
	double EpsilonEl = 0; 
	double YoungModulus = 70.00000e+09;

// Geberate Basis //////////////////////////////////////////
	Basis[0].i1 = 1; Basis[1].i1 = 0; Basis[2].i1 = 0;
	Basis[0].i2	= 0; Basis[1].i2 = 1; Basis[2].i2 = 0;
	Basis[0].i3	= 0; Basis[1].i3 = 0; Basis[2].i3 = 1;

// Generate 
	Z.i1 = 0;
	Z.i2 = 0;
	Z.i3 = 1;

	hkl.i1 = 2; 
	hkl.i2 = 1;
	hkl.i3 = 3;
	
	printf ("hkl = %e %e %e\n",hkl.i1, hkl.i2, hkl.i3);
	
	Norm = sqrt(hkl.i1 * hkl.i1 + hkl.i2 * hkl.i2 + hkl.i3 * hkl.i3);
	hkl.i1 = hkl.i1 /Norm;
	hkl.i2 = hkl.i2 /Norm;
	hkl.i3 = hkl.i3 /Norm;
	
	printf ("hkl = %e %e %e norm\n",hkl.i1, hkl.i2, hkl.i3);

		FILE * pFile;
		pFile = fopen ("C:\dat.txt","w");

		FILE * dFile;
		dFile = fopen ("C:\debug1.txt","w");

////////////////////////////////////////////////////////////////
	
	

	prepareCurveConstants();

// Work Space ********************* OPERATORS

	GenerateSterssTensor (); 
	GenerateSetOfSlipSystems ();
	GenerateG_operator ();

// Print G_operator 

	printf ( "** G_operator **\n");

	printf ( " %e  %e  %e  \n", g_operator.i11, g_operator.i12, g_operator.i13 );
	printf ( " %e  %e  %e  \n", g_operator.i21, g_operator.i22, g_operator.i23 );
	printf ( " %e  %e  %e  \n", g_operator.i31, g_operator.i32, g_operator.i33 );

// Cristall in work position (Slip System + Basis)

	for (i=0; i<24; i++) 
	{
	n[i] = MatrixVector (g_operator, n1[i]);
	b[i] = MatrixVector (g_operator, b1[i]);
	}


	for (i=0; i<3; i++) { Basis[i] = MatrixVector (g_operator, Basis[i]); }

	printf ( "** Basis **\n");

	printf ( " 1x %e  2x %e  3x %e  \n", Basis[0].i1, Basis[1].i1, Basis[2].i1 );
	printf ( " 1y %e  2y %e  3y %e  \n", Basis[0].i2, Basis[1].i2, Basis[2].i2 );
	printf ( " 1z %e  2z %e  3z %e  \n", Basis[0].i3, Basis[1].i3, Basis[2].i3 );

	printf ("Information\n" );

	printf ("MachineStrainRate  = %e\n", MachineStrainRate );
	printf ("exper.T_Kelvin = %e\n", exper.T_Kelvin);
	printf ("tauSolTherm0   = %e\n", val.tauSolTherm0);

    fprintf (pFile,"%e %e\n", 0.0, 0.0 ); // zero poin in dat.txt file

//   *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
//   *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** 

// dt calculatinon ( dt - time step)

	if (MachineStrainRate >= 0 )
	{dt= 1/MachineStrainRate * 1e-06;}
	else 
	{dt=-1/MachineStrainRate * 1e-06;}
	
// main ziklus

	for (k1=0; k1<100000; k1++)
	for (k2=0; k2<10; k2++)
	
	{
	{

// Generation and calculation of Stresstensor

	GenerateSterssTensor ();

	StressTensor.i33 = YoungModulus * EpsilonEl + YoungModulus * MachineStrainRate * dt;

// calculate Stress In Slip Systems (in All SS)

	for (i=0; i<24; i++) { StressInSlipSystem[i] = ShearStress (n[i], b[i], StressTensor);}	

// calculate Shear In Slip Systems


	for (i=0; i<24; i++) 
	{

	if (StressInSlipSystem[i] > 0) 
	{ShearInSlipSystem[i] = val.burgersVec * 0.083 * rho.mobil * v_SS(StressInSlipSystem[i]) * dt;}
	else 
	{ShearInSlipSystem[i] = 0;};

	}


// calculate Sigma Shear

	SigmaShear = 0;
	for (i=0; i<24; i++) { SigmaShear = SigmaShear + ShearInSlipSystem[i];}
	

	// key for 3IVM

	if (SigmaShear > 0) { 


// calculate ShearTensor

	for (i=0; i<24; i++) { ShearTensor[i] = TensorProduct(n[i], b[i]);}	

// calculate StrainTensor

	GenerateStrainTensor();

	for (i=0; i<24; i++)
	{
	StrainTensor.i11 = StrainTensor.i11 + ShearInSlipSystem[i] * ShearTensor[i].i11;
	StrainTensor.i12 = StrainTensor.i12 + ShearInSlipSystem[i] * ShearTensor[i].i12;
	StrainTensor.i13 = StrainTensor.i13 + ShearInSlipSystem[i] * ShearTensor[i].i13;

	StrainTensor.i21 = StrainTensor.i21 + ShearInSlipSystem[i] * ShearTensor[i].i21;
	StrainTensor.i22 = StrainTensor.i22 + ShearInSlipSystem[i] * ShearTensor[i].i22;
	StrainTensor.i23 = StrainTensor.i23 + ShearInSlipSystem[i] * ShearTensor[i].i23;

	StrainTensor.i31 = StrainTensor.i31 + ShearInSlipSystem[i] * ShearTensor[i].i31;
	StrainTensor.i32 = StrainTensor.i32 + ShearInSlipSystem[i] * ShearTensor[i].i32;
	StrainTensor.i33 = StrainTensor.i33 + ShearInSlipSystem[i] * ShearTensor[i].i33;
	}

// calculate num Active Planes


	double maxShear = ShearInSlipSystem[0];
	
	for (i=0; i<24; i++)
	{
		if (ShearInSlipSystem[i] > maxShear){ maxShear = ShearInSlipSystem[i];}
	}

	fit.numActivePlanes = SigmaShear/maxShear;

// ***********************************

	EpsilonEl = EpsilonEl + MachineStrainRate * dt - StrainTensor.i33;

//	exper.strainRate = StrainTensor.i33/dt;
	
	Epsilon = Epsilon + MachineStrainRate * dt; // exper_strainRate;

	exper.simSigma = YoungModulus * EpsilonEl;

	if (k2 == 0) 
	{
		fprintf (pFile,"%e %e ", Epsilon, exper.simSigma );
		fprintf (pFile,"\n");
	}	

	
// the 3IVM main function
	
		prepareIntegrationConstants();
		long nonsense = integrateRho (dt);

// Calculation of new cristallposition

// calculate RotateTensor 
	
	for (i=0; i<24; i++)
	{
	r_operator[i].i11 = 0; r_operator[i].i12 = 0; r_operator[i].i13 = 0;
	r_operator[i].i21 = 0; r_operator[i].i22 = 0; r_operator[i].i23 = 0;
	r_operator[i].i31 = 0; r_operator[i].i32 = 0; r_operator[i].i33 = 0;
	}
	
		
	for (i=0; i<24; i++) { RotationTensor[i] = _TensorProduct( b[i], n[i]);}

	
	for (i=0; i<24; i++)
	{
	r_operator[i].i11 = ShearInSlipSystem[i] * RotationTensor[i].i11 + 1;
	r_operator[i].i12 = ShearInSlipSystem[i] * RotationTensor[i].i12;
	r_operator[i].i13 = ShearInSlipSystem[i] * RotationTensor[i].i13;

	r_operator[i].i21 = ShearInSlipSystem[i] * RotationTensor[i].i21;
	r_operator[i].i22 = ShearInSlipSystem[i] * RotationTensor[i].i22 + 1;
	r_operator[i].i23 = ShearInSlipSystem[i] * RotationTensor[i].i23;

	r_operator[i].i31 = ShearInSlipSystem[i] * RotationTensor[i].i31;
	r_operator[i].i32 = ShearInSlipSystem[i] * RotationTensor[i].i32;
	r_operator[i].i33 = ShearInSlipSystem[i] * RotationTensor[i].i33 + 1;
	}

//
	for (j=0; j<24; j++)
	{
	for (i=0; i<24; i++) 
	{
	n[i] = MatrixVector (r_operator[j], n[i]);
	b[i] = MatrixVector (r_operator[j], b[i]);


	}
	}

	for (j=0; j<24; j++)
	{
	for (i=0; i<3; i++) 
	{
	Basis[i] = MatrixVector (r_operator[j], Basis[i]);
	}
	}

// Calculation of orintation for triangle presantation

	hkl.i1 = ScalarProduct (Basis[0], Z);
	hkl.i2 = ScalarProduct (Basis[1], Z);
	hkl.i3 = ScalarProduct (Basis[2], Z);

// Output of orintation in 'debug1.txt'
	if (k2 == 0) 
	{fprintf (dFile,"%e %e\n", hkl.i1, hkl.i2 );}


	} // end if (key for 3IVM)
	else
	{

	EpsilonEl = EpsilonEl + MachineStrainRate * dt;

	Epsilon = Epsilon + MachineStrainRate * dt; // otnositel'noe udlinnenie;

	exper.simSigma = YoungModulus * EpsilonEl;

// Output in filt 'dat.txt'
	if (k2 == 0) {fprintf (pFile,"%e %e\n", Epsilon, exper.simSigma );}
	

	} // end if (key for 3IVM)

	}// *****************
	}// *****************

		fclose (dFile);	
		fclose (pFile);
		printf ("The End\n" );

		return 0;
} // end main(void)	


// ********************************************************
//			FUNCTIONS
// ********************************************************

//  funktion prepareCurveConstants

void prepareCurveConstants(void)
{
	val._burgersVec = 1 / val.burgersVec;
	val._Lgrain0 = 1 / exper.Lgrain0;
	val._fWall = 1 / fit.fWall;
	val.fIntern = 1 - fit.fWall;

		rho.mobil  = fit.rhoMobil0;
		rho.intern = fit.rhoIntern0;
		rho.wall   = fit.rhoWall0;

		double tsf1q = fit.thermSolStrength1 * fit.thermSolStrength1;
		double tsf2q = fit.thermSolStrength2 * fit.thermSolStrength2;
		double asf1q = fit.athermSolStrength1 * fit.athermSolStrength1;
		double asf2q = fit.athermSolStrength2 * fit.athermSolStrength2;
		
		double c1 = exper.atConc1;
		double c2 = exper.atConc2;
	
	val.tauSolTherm0 = sqrt( tsf1q * c1 + tsf2q * c2 );
	val.tauSolAtherm = sqrt( asf1q * c1 + asf2q * c2 );
	val._lambdaSol = 1 / (val.burgersVec / sqrt( c1 + c2 ));

		double b = val.burgersVec;
		double r_b = exper.particleRadius_b;
		double gamma_b = fit.particleStrength;
		double f = exper.particleVolumeFraction; 
		double wf = sqrt((4/_PI_) * f);
		double tauUnit = wf * gamma_b;
		double fi = 2*r_b*gamma_b / val.shearMod0;
		double tauUnder = (0.9-3*f)*sqrt(fi) + 0.4*wf;
		double tauOver = 0.8/(fi*(1-wf));
		double tau = 0.8 + f;						// peak
	
		if(tauUnder < tau)		tau = tauUnder;
		if(tauOver < tau)		tau = tauOver;
	
	val.tauParticle = tau * tauUnit;
	val._numActivePlanes = 1 / fit.numActivePlanes;

		
		double T_Kelvin = exper.T_Kelvin;

	val.kT = KBOLTZMANN * T_Kelvin;
	val._kT = 1 / val.kT;
	val.shearMod = val.shearMod0 + (T_Kelvin - val.T0G_Kelvin) * val.dG_dT;

	val.Gb_4p = val.shearMod * val.burgersVec * (0.25/_PI_);
	val.ttFac = val.shearMod * val.burgersVec * fit.alphaTaylor;

		double clf = val.burgersVec * val.burgersVec * val.diffCoeff0 * val._kT;
	
	val.climbFac = 2 * fit.dClimb * val._numActivePlanes * clf * exp(- fit.Qclimb * val._kT);

		double db = val.burgersVec * fit.d0solute;
	val._DGsol0 = 1 / sqrt( db*db*db * val.shearMod * val.tauSolTherm0 );

}


	
void prepareIntegrationConstants()	
{
	val.MepsDot_b = SigmaShear/dt * val._burgersVec;
	
		double MepsDot2_bn = val.MepsDot_b * 2 * val._numActivePlanes;
	
	val.annFac = MepsDot2_bn * fit.dAnn;
	val.lockFac = 2 * MepsDot2_bn * fit.dLock * (fit.numActivePlanes - 1);
	val.dipolFac = MepsDot2_bn;
	val.immobFac = MepsDot2_bn * fit.dImmob;
	val.clearFac = MepsDot2_bn * fit.dClear;

}

// funktion diffEquations


long diffEquations( void )
{
	// hier sind T, fw, G, b, ... konstant
	// davon abhngiges ist extern in val vorgekaut

	if(rho.mobil < 0)			return 1;
	if(rho.intern < 0)			return 1;
	if(rho.wall < 0)			return 1;
	if(rho.intern > 1e15)		return 1;
	if(rho.wall > 1e16)			return 1;

		const double  nyy0 = 3e+10;			// Anlauffrequenz 10^10...10^11 Hz
		const double _nyy0 = 1 / nyy0;

		double rhoWallTotal   = rho.wall   + rho.mobil * fit.fWall;
		double rhoInternTotal = rho.intern + rho.mobil * val.fIntern;

		double _lambdaWall   =  sqrt(rhoWallTotal);
		double _lambdaIntern =  sqrt(rhoInternTotal);

		double tauWallTaylor   = val.ttFac * _lambdaWall;
		double tauInternTaylor = val.ttFac * _lambdaIntern;

		double vGlide = val.MepsDot_b / rho.mobil;		// Orowan-Gleichung
		vG = val.MepsDot_b / rho.mobil;
		double DGsol = - val.kT * log( vGlide * val._lambdaSol * _nyy0 );
	
		double tauSolRel = 2 / ( 2 + ( DGsol * val._DGsol0 ) );
	
	if(vGlide < 1e-30)	tauSolRel = 0;

		double tauSolTherm =  val.tauSolTherm0 * tauSolRel;
		double tauChem = tauSolTherm + val.tauSolAtherm + val.tauParticle;
	
		double tauWall   = 	tauWallTaylor 	+ tauChem;
		double tauIntern = 	tauInternTaylor + tauChem;

//	val.simTau = fit.fWall * tauWall + val.fIntern * tauIntern;

	//	double _Leff;
		double _Leff = val._Lgrain0 * 2;

		

	_Leff += fit._betaLeffIntern * _lambdaIntern * val.fIntern;
	_Leff += fit._betaLeffWall   * _lambdaWall   * fit.fWall;
	
		double dDipol = val.Gb_4p / tauWallTaylor - fit.dAnn;

		double rdmOrowan = val.MepsDot_b * _Leff;

		double rdmAnn = val.annFac * rho.mobil;	
	
		double rdwImmob = val.immobFac * rho.wall;
		double rdiImmob = val.immobFac * rho.intern;
		double rdmImmob = fit.fWall * rdwImmob + val.fIntern * rdiImmob;

		double rdmLock = val.lockFac * rho.mobil;	
		double rdwLock = rdmLock;
		double rdiLock = rdmLock;
	
		double rdmDipol = val.dipolFac * (dDipol - fit.dAnn) * rho.mobil;	
		double rdwDipol = rdmDipol * val._fWall;

		double rdwClear = val.clearFac * rho.wall;
		double rdiClear = val.clearFac * rho.intern;

		double gi = fit.Qcross - fit.Vcross * tauIntern;
		double gw = fit.Qcross - fit.Vcross * tauWall;

		double rdiCross = rho.intern * nyy0 * exp( - gi * val._kT );
		double rdwCross = rho.wall * nyy0 * exp( - gw * val._kT );

double clf = val.climbFac * (val.MepsDot_b * val.burgersVec) * (val.MepsDot_b * val.burgersVec) * 10 * 1e12;
	
	double rdiClimb = val.climbFac * tauInternTaylor * rho.intern * rho.intern;
	double rdwClimb = val.climbFac * tauWallTaylor   * rho.wall   * rho.wall;

	rhoDot.mobil  = rdmOrowan - rdmAnn - rdmLock - rdmDipol - rdmImmob;
	rhoDot.intern = rdiLock + rdiImmob - rdiClear - rdiCross - rdiClimb;
	rhoDot.wall   = rdwLock + rdwImmob + rdwDipol - rdwClear - rdwCross - rdwClimb;

	return 0;
}


// funktion integrateRho 

long integrateRho( double timeStep )
{
	long nonsense;
	
	nonsense = diffEquations();
	if( nonsense )	return 1;
	
	double drm	= timeStep * rhoDot.mobil;
	double dri	= timeStep * rhoDot.intern;
	double drw	= timeStep * rhoDot.wall;

	rho.mobil	+= drm;
	rho.intern	+= dri;
	rho.wall	+= drw;

	return 0;
}

// funktion v_SS 

	double v_SS(double tauSlipSystem)
	{
	
		double v_SS;
		double vWall;
		double vIntern;

		
		const double nyy0 = 3e+10;			// Anlauffrequenz 10^10...10^11 Hz

		double rhoWallTotal   = rho.wall + rho.mobil * fit.fWall;
		double rhoInternTotal = rho.intern + rho.mobil * val.fIntern;

		double _lambdaWall   =  sqrt(rhoWallTotal);
		double _lambdaIntern =  sqrt(rhoInternTotal);

		double tauWallTaylor   = val.ttFac * _lambdaWall;
		double tauInternTaylor = val.ttFac * _lambdaIntern;

		double _Leff;

_Leff = fit._betaLeffIntern * _lambdaIntern * val.fIntern + fit._betaLeffWall * _lambdaWall * fit.fWall;
/////////////////////////////////////////////////////////////////////////////
	
//		val.lambdaSol = val.burgersVec / sqrt(0.01);
		
		val.lambdaSol = 1/_Leff;

/////////////////////////////////////////////////////////////////////////////

		if ((tauSlipSystem - tauWallTaylor - val.tauSolAtherm) <= 0 ) vWall = 0;

		if ((tauSlipSystem - tauWallTaylor - val.tauSolAtherm) >  0 )
		
		{

		double tauRelativ = (tauSlipSystem - tauWallTaylor - val.tauSolAtherm) / val.tauSolTherm0;

		double gWall = 2 * 1/tauRelativ * exp(1.25 * log(1 - tauRelativ));

		double GWall = val.burgersVec * val.burgersVec * val.burgersVec * sqrt(val.shearMod0/2 * val.tauSolTherm0) * gWall;
		
		vWall   = val.lambdaSol * nyy0 * exp ( - GWall * val._kT);

		}
///////////////////////////////////////////////////////////////////////////////
		
		if ((tauSlipSystem - tauInternTaylor - val.tauSolAtherm) <= 0 ) vIntern = 0;

		if ((tauSlipSystem - tauInternTaylor - val.tauSolAtherm) >  0 )
		
		{

		double tauRelativ = (tauSlipSystem - tauInternTaylor - val.tauSolAtherm) / val.tauSolTherm0;

		double gIntern = 2 * 1/tauRelativ * exp(1.25 * log(1 - tauRelativ));

		double GIntern = val.burgersVec * val.burgersVec * val.burgersVec * sqrt(val.shearMod0/2 * val.tauSolTherm0) * gIntern;

		vIntern = val.lambdaSol * nyy0 * exp ( - GIntern * val._kT);
		
		}
		
		v_SS = fit.fWall * vWall + (1 - fit.fWall) * vIntern;
		
		return v_SS;
	}
//	Funktion ///////////////////////////////////////////////////////////////////

	Matrix TensorProduct(Vector n, Vector b)
	{
		Matrix c;

		c.i11 = 0.5 * (b.i1 * n.i1 + n.i1 * b.i1);
		c.i12 = 0.5 * (b.i1 * n.i2 + n.i1 * b.i2);	
		c.i13 = 0.5 * (b.i1 * n.i3 + n.i1 * b.i3);

		c.i21 = 0.5 * (b.i2 * n.i1 + n.i2 * b.i3);
		c.i22 = 0.5 * (b.i2 * n.i2 + n.i2 * b.i2);
		c.i23 = 0.5 * (b.i2 * n.i3 + n.i2 * b.i3);

		c.i31 = 0.5 * (b.i3 * n.i1 + n.i3 * b.i1);
		c.i32 = 0.5 * (b.i3 * n.i2 + n.i3 * b.i2);
		c.i33 = 0.5 * (b.i3 * n.i3 + n.i3 * b.i3);
		
		return c;
	}

///////////////////////////////////////////////////////////////////////////////
	Matrix _TensorProduct( Vector b, Vector n)
	{
		Matrix c;

		c.i11 = 0.5*(b.i1 * n.i1 - n.i1 * b.i1); 
		c.i12 = 0.5*(b.i2 * n.i1 - n.i2 * b.i1);
		c.i13 = 0.5*(b.i3 * n.i1 - n.i3 * b.i1);

		c.i21 = 0.5*(b.i1 * n.i2 - n.i1 * b.i2);
		c.i22 = 0.5*(b.i2 * n.i2 - n.i2 * b.i2);
		c.i23 = 0.5*(b.i3 * n.i2 - n.i3 * b.i2);

		c.i31 = 0.5*(b.i1 * n.i3 - n.i1 * b.i3);
		c.i32 = 0.5*(b.i2 * n.i3 - n.i2 * b.i3);
		c.i33 = 0.5*(b.i3 * n.i3 - n.i3 * b.i3);

		return c;
	}
///////////////////////////////////////////////////////////////////////////////
	double ScalarProduct(Vector a, Vector b)
	{
		double c;
	
		c = a.i1 * b.i1 + a.i2 * b.i2 + a.i3 * b.i3;
		
		return c;
	}


/////////////////////////////////////////////////////////////////////////////
	double ShearStress (Vector n, Vector b, Matrix Sigma)
	{
		double SS,SS1,SS2,SS3;

		SS1 = b.i1*(Sigma.i11*n.i1 + Sigma.i12*n.i2 + Sigma.i13*n.i3);
		SS2 = b.i2*(Sigma.i21*n.i1 + Sigma.i22*n.i2 + Sigma.i23*n.i3);
		SS3 = b.i3*(Sigma.i31*n.i1 + Sigma.i32*n.i2 + Sigma.i33*n.i3);
		
		SS = SS1 + SS2 +SS3;
		return SS;
	}

//////////////////////////////////////////////////////////////////////////////
	Vector MatrixVector (Matrix mat, Vector vec)
	{
		Vector Prod;

		Prod.i1 = mat.i11 * vec.i1 + mat.i12 * vec.i2 + mat.i13 * vec.i3;
		Prod.i2 = mat.i21 * vec.i1 + mat.i22 * vec.i2 + mat.i23 * vec.i3;
		Prod.i3 = mat.i31 * vec.i1 + mat.i32 * vec.i2 + mat.i33 * vec.i3;
		return Prod;
	}

//////////////////////////////////////////////////////////////////////////////
	Vector MatrixVectorReverse (Matrix mat, Vector vec)
	{
		Vector Prod;

		Prod.i1 = mat.i11 * vec.i1 + mat.i21 * vec.i2 + mat.i31 * vec.i3;
		Prod.i2 = mat.i12 * vec.i1 + mat.i22 * vec.i2 + mat.i32 * vec.i3;
		Prod.i3 = mat.i13 * vec.i1 + mat.i23 * vec.i2 + mat.i33 * vec.i3;
		return Prod;
	}

//////////////////////////////////////////////////////////////////////////////
void GenerateSetOfSlipSystems (void)
	{
	
	n1[0].i1 =  1;  b1[0].i1 =  0;
	n1[0].i2 =  1;  b1[0].i2 =  1;
	n1[0].i3 =  1;  b1[0].i3 = -1;

	n1[1].i1 =  1;  b1[1].i1 =  0;
	n1[1].i2 =  1;  b1[1].i2 = -1;
	n1[1].i3 =  1;  b1[1].i3 =  1;
//
	n1[2].i1 =  1;  b1[2].i1 =  1;
	n1[2].i2 =  1;  b1[2].i2 =  0;
	n1[2].i3 =  1;  b1[2].i3 = -1;

	n1[3].i1 =  1;  b1[3].i1 = -1;
	n1[3].i2 =  1;  b1[3].i2 =  0;
	n1[3].i3 =  1;  b1[3].i3 =  1;
//
	n1[4].i1 =  1;  b1[4].i1 =  1;
	n1[4].i2 =  1;  b1[4].i2 = -1;
	n1[4].i3 =  1;  b1[4].i3 =  0;

	n1[5].i1 =  1;  b1[5].i1 = -1;
	n1[5].i2 =  1;  b1[5].i2 =  1;
	n1[5].i3 =  1;  b1[5].i3 =  0;

	/////////////////////////////////////////////////////////

	n1[6].i1 = -1;  b1[6].i1 =  0;
	n1[6].i2 =  1;  b1[6].i2 =  1;
	n1[6].i3 =  1;  b1[6].i3 = -1;

	n1[7].i1 = -1;  b1[7].i1 =  0;
	n1[7].i2 =  1;  b1[7].i2 = -1;
	n1[7].i3 =  1;  b1[7].i3 =  1;
//
	n1[8].i1 = -1;  b1[8].i1 =  1;
	n1[8].i2 =  1;  b1[8].i2 =  0;
	n1[8].i3 =  1;  b1[8].i3 =  1;

	n1[9].i1 = -1;  b1[9].i1 = -1;
	n1[9].i2 =  1;  b1[9].i2 =  0;
	n1[9].i3 =  1;  b1[9].i3 = -1;
//
	n1[10].i1 = -1;  b1[10].i1 =  1;
	n1[10].i2 =  1;  b1[10].i2 =  1;
	n1[10].i3 =  1;  b1[10].i3 =  0;

	n1[11].i1 = -1;  b1[11].i1 = -1;
	n1[11].i2 =  1;  b1[11].i2 = -1;
	n1[11].i3 =  1;  b1[11].i3 =  0;

	//////////////////////////////////////////////////////////
	
	n1[12].i1 =  1;  b1[12].i1 =  0;
	n1[12].i2 = -1;  b1[12].i2 =  1;
	n1[12].i3 =  1;  b1[12].i3 =  1;

	n1[13].i1 =  1;  b1[13].i1 =  0;
	n1[13].i2 = -1;  b1[13].i2 = -1;
	n1[13].i3 =  1;  b1[13].i3 = -1;
//
	n1[14].i1 =  1;  b1[14].i1 =  1;
	n1[14].i2 = -1;  b1[14].i2 =  0;
	n1[14].i3 =  1;  b1[14].i3 = -1;

	n1[15].i1 =  1;  b1[15].i1 = -1;
	n1[15].i2 = -1;  b1[15].i2 =  0;
	n1[15].i3 =  1;  b1[15].i3 =  1;
//
	n1[16].i1 =  1;  b1[16].i1 =  1;
	n1[16].i2 = -1;  b1[16].i2 =  1;
	n1[16].i3 =  1;  b1[16].i3 =  0;

	n1[17].i1 =  1;  b1[17].i1 = -1;
	n1[17].i2 = -1;  b1[17].i2 = -1;
	n1[17].i3 =  1;  b1[17].i3 =  0;

	//////////////////////////////////////////////////////////

	n1[18].i1 = -1;  b1[18].i1 =  0;
	n1[18].i2 = -1;  b1[18].i2 =  1;
	n1[18].i3 =  1;  b1[18].i3 =  1;

	n1[19].i1 = -1;  b1[19].i1 =  0;
	n1[19].i2 = -1;  b1[19].i2 = -1;
	n1[19].i3 =  1;  b1[19].i3 = -1;
//
	n1[20].i1 = -1;  b1[20].i1 =  1;
	n1[20].i2 = -1;  b1[20].i2 =  0;
	n1[20].i3 =  1;  b1[20].i3 =  1;

	n1[21].i1 = -1;  b1[21].i1 = -1;
	n1[21].i2 = -1;  b1[21].i2 =  0;
	n1[21].i3 =  1;  b1[21].i3 = -1;
//
	n1[22].i1 = -1;  b1[22].i1 =  1;
	n1[22].i2 = -1;  b1[22].i2 = -1;
	n1[22].i3 =  1;  b1[22].i3 =  0;

	n1[23].i1 = -1;  b1[23].i1 = -1;
	n1[23].i2 = -1;  b1[23].i2 =  1;
	n1[23].i3 =  1;  b1[23].i3 =  0;

///////////////////////////////////////////////////////

		for (i=0; i<24; i++)
		{
		n1[i].i1 = n1[i].i1/1.732050;
		n1[i].i2 = n1[i].i2/1.732050;
		n1[i].i3 = n1[i].i3/1.732050;
   		
		b1[i].i1 = b1[i].i1/1.414213;
		b1[i].i2 = b1[i].i2/1.414213;
		b1[i].i3 = b1[i].i3/1.414213;
		}

	return;
	}

/////////////////////////////////////////////////////////////////////////

	void GenerateSterssTensor (void)
	{
	StressTensor.i11 = 0;
	StressTensor.i12 = 0;
	StressTensor.i13 = 0;
	StressTensor.i21 = 0; 
	StressTensor.i22 = 0;
	StressTensor.i23 = 0;
	StressTensor.i31 = 0;
	StressTensor.i32 = 0;
	StressTensor.i33 = 0;

	return;
	}
/////////////////////////////////////////////////////////////////////////

	void GenerateStrainTensor (void)
	{
	StrainTensor.i11 = 0;
	StrainTensor.i12 = 0;
	StrainTensor.i13 = 0;
	StrainTensor.i21 = 0; 
	StrainTensor.i22 = 0;
	StrainTensor.i23 = 0;
	StrainTensor.i31 = 0;
	StrainTensor.i32 = 0;
	StrainTensor.i33 = 0;

	return;
	}
//////////////////////////////////////////////////////////////////////////
	void GenerateG_operator (void)
	{
	double h_index = hkl.i1;    
	double k_index = hkl.i2;    
	double l_index = hkl.i3; 

	double u_index;
	double v_index;
	double w_index;


	if ( (h_index == 0) && (k_index == 0) && (l_index == 1) ) 
	{
		u_index = 1;
		v_index = 0;
		w_index = 0;
	}
	else
	{
		u_index = - k_index / h_index;
		v_index =  1;
		w_index =  0;
	}

	double mod_hkl = sqrt(h_index * h_index + k_index * k_index + l_index * l_index);
	double mod_uvw = sqrt(u_index * u_index + v_index * v_index + w_index * w_index);

	h_index = h_index/mod_hkl;    
	k_index = k_index/mod_hkl;    
	l_index = l_index/mod_hkl; 

	u_index = 1/mod_uvw * u_index;
	v_index = 1/mod_uvw * v_index;
	w_index = 1/mod_uvw * w_index;

	double q_index =  k_index * w_index - l_index * v_index;
	double r_index = -h_index * w_index + l_index * u_index;
	double s_index =  h_index * v_index - k_index * u_index;
	
	double mod_qrs = sqrt(q_index * q_index + r_index * r_index + s_index * s_index);

	q_index = q_index/mod_qrs;
	r_index = r_index/mod_qrs;
	s_index = s_index/mod_qrs;

		g_operator.i11 =  u_index;  g_operator.i12 =  v_index;  g_operator.i13 =  w_index;  
		g_operator.i21 =  q_index;  g_operator.i22 =  r_index;  g_operator.i23 =  s_index; 
		g_operator.i31 =  h_index;  g_operator.i32 =  k_index;  g_operator.i33 =  l_index;


	return;
	}