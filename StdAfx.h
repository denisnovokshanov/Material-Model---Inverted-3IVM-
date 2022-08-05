// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#if !defined(AFX_STDAFX_H__F2C15E41_5DE5_4020_B5D2_E09739CB9888__INCLUDED_)
#define AFX_STDAFX_H__F2C15E41_5DE5_4020_B5D2_E09739CB9888__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

#include <stdio.h>

// TODO: reference additional headers your program requires here

struct Vector
{
double i1;
double i2;
double i3;
};

struct Matrix
{
double i11; double i12; double i13;
double i21; double i22; double i23;
double i31; double i32; double i33;
};

// Inizialisation exper 

struct curveData
{

	float expSigma;
	float simSigma;

	float MTaylor;
	float T_Kelvin;

	float strain;
	float time;
	float strainRate;

	float rhoM;
	float rhoW;
	float rhoI;


	double Lgrain0;
	double atConc1;
	double atConc2;
	double particleVolumeFraction;
	double particleRadius_b;

};

struct fitParmSet
{
	double dAnn;   // [m]
	double dLock;  // [m]
	double dImmob; // [m]
	double dClear; // [m]
	double dClimb; // [m]

	double fWall;

	double betaLeffIntern;
	double betaLeffWall;
	
	double _betaLeffIntern;
	double _betaLeffWall;
	
	double alphaTaylor;

	double d0solute; // [m]
	double Vcross  ; // [m^3]

	double Qcross;   // [J]
	double Qclimb;   // [J]
	
	double rhoMobil0;   //8.282318e+10; [m^-2]
	double rhoIntern0;  //2.455309e+09; [m^-2]
	double rhoWall0;    //5.778842e+13; [m^-2]
	
	double thermSolStrength1; // [MPa]
	double thermSolStrength2; // [MPa]
	
	double athermSolStrength1;// [MPa]
	double athermSolStrength2;// [MPa]

	double particleStrength;  // [MPa]
	double numActivePlanes;
};

struct simData
{
	double burgersVec;  // [m]	// real constants
	double _burgersVec; // [1/m]
	double shearMod0;   // [Pa]
	double T0G_Kelvin;  // [K]
	double dG_dT;       // [Pa/K]
	double Lgrain0;     // [m]
	double _Lgrain0;    // [m^-1]

	double diffCoeff0;  // [m^2/s]
	
	double shearMod;	// T-, epsp-, M- dependent
	double kT;
	double _kT;
	double Gb_4p;
	double ttFac;
	double MepsDot_b;

	double fWall;		// fitParm-dependent
	double _fWall;
	double fIntern;
	
	double tauSolTherm0;	// chemistry-dependent
	double tauSolAtherm;

	double lambdaSol;
	double _lambdaSol;
	double tauParticle;
	
	double _DGsol0;
	double annFac;
	double lockFac;
	double dipolFac;
	double immobFac;
	double clearFac;
	double climbFac;

	double numActivePlanes;
	double _numActivePlanes;
	double simTau;
};

struct  allRho
{
	double mobil;
	double intern;
	double wall;
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__F2C15E41_5DE5_4020_B5D2_E09739CB9888__INCLUDED_)
