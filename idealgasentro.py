#-*- coding: utf-8 -*-
#  AMDG
# JMJ-V!
# 
# -------------------------------------------------------------------
#  Entropie als Funktion von Y, rho, T für ein ideales Gas
#  
#   30.10.2019 (c) Gabriel-Dominique Marleau, Uni Tübingen
#   gabriel.marleau@uni-tuebingen.de
# 
# 
#  
#  * Gilt im Moment nur für T >> thetarot ~ 85 K ! *
#    vgl. Abb. 2 von Vaytet et al. (2014) -> wohl ca. wie Gleichgewichtsverhältnis
# -------------------------------------------------------------------

import numpy as np
from scipy.special import gamma

# 
from scipy.optimize import brentq

# nur in (i)python sinnvoll:
#%load_ext autoreload
#%autoreload 2
# this will reload automatically when the library (the source file idealgasentro.py) is changed

# Naming convention:
#  - Species:
#     Htot: all hydrogen species
#           H2, HI, HII, eH: H_2, neutral H, ionised H, "hydrogen electrons"
#     HeI, HeII, HeIII, eHe: similar
# 
#  - Functions without the suffix _rhoT are functions of the X_i's (provided by SCvH tables
#    or to be computed separately)
# 
#  - Uppercase S for entropy per molecule, lowercase s for entropy per gram
# 
# Notes:
#  Following SCvH, the hydrogen and the helium parts can be computed separately, as though
#    there were only species of that element, and be combined afterwards.
# 
#  All entropies are either in kB/molecule or kB/mH (i.e., the kB is always dropped
#  
#  There are different layers of functions:
#    1.) x,y,z1,z2 as functions of rho and T
#    2.) X_i's (from SCvH) as a function of x,y,z1,z2 (from D'AB13)
#    3.) thermodynamic functions as functions of the X_i's
#    4.) thermodynamic functions as functions of rho and T
#   
#   Usually only functions of level 4 are called by the user but the other ones
#   can be called too without problem.
#  
#  The hydrogen and helium functions could be perfectly symmetric since the species
#  are but for historical reasons there are slight differences. This is only a minor
#  cosmetico-programmatical issue, though, which does not affect the answers!
# 

# ===============================================================================
#  The dissociation and ionisation fractions from D'Angelo & Bodenheimer (2013)
# ===============================================================================
# 
## Eqs. (15--19) of D'Angelo & Bodenheimer (2013), http://adsabs.harvard.edu/abs/2013ApJ...778...77D
##   The quantities are defined above their Eq. (15) and are:
#
# x    Degree of ionisation of atomic hydrogen = rho_{H+} / ( rho_{H}  + rho_{H+} ),
# y    Degree of dissociation of molecular hydrogen = rho_{H}  / ( rho_{H2} + rho_{H}  )
# z1   Degree of single ionisation of helium = rho_{He+}  / ( rho_{He}   + rho_{He+}  )
# z2   Degree of double ionisation of helium = rho_{He2+} / ( rho_{He+}  + rho_{He2+} )
#
#    If there are actually metals (Z!=0), they should be counted in the helium part (Y->Yeff = Y+Z, Z->0)
# 
###  Note: Putting a*a in the sqrt instead of a**2. avoids the need to truncate at low a values.
###        Still need to truncate at high values by putting "?:" to avoid numerical problems
###        Probably because of the exponential (truncates at y=1-1e-8 at high T for rho=1e-13 and Y=0.25)
###        Similar for z1 and z2 (maybe, not quite tested)
# 

def yDAB13(Y,rho,T,Z=0):
	DeltaH2HI = 4.48
	
	#if (Y+Z == 1.):
		#return 1.0
	
	# with a*a instead of a**2.: can also get low abundances! Same for z1 and z2 below
	a = amu/2/(1.-Y-Z)/rho * (pi*amu*kB*T)**1.5/hPlanck**3 *np.exp(-DeltaH2HI*eV/(kB*T))
	
	# needed:
	return 1.0 if a>1e8 else 0.5*(-a+np.sqrt(a*a+4*a))

def xDAB13(Y,rho,T,Z=0):
	DeltaHIHII = 13.6
	
	#if (Y+Z == 1.):
		#return 1.0
	
	a = amu/(1.-Y-Z)/rho * (2.*pi*me*kB*T)**1.5/hPlanck**3 *np.exp(-DeltaHIHII*eV/(kB*T))
	
	return 1.0 if a>1e8 else 0.5*(-a+np.sqrt(a*a+4*a))

def z1DAB13(Y,rho,T,Z=0):
	DeltaHeIHeII = 24.59
	
	X = 1.-Y-Z
	
	#if (X == 0.):
		#return 1.0
		
	b = 4*amu/rho * (2*pi*me*kB*T)**1.5/hPlanck**3 *np.exp(-DeltaHeIHeII*eV/(kB*T))
	
	return 1.0 if b>1e12*Y else 1./(0.5*Y)*(-(b+X) + np.sqrt((b+X)*(b+X) + Y*b))

def z2DAB13(Y,rho,T,Z=0):
	DeltaHeIIHeIII = 54.42
	
	X = 1.-Y-Z
	
	#if (X == 0.):
		#return 1.0
		
	b = amu/rho * (2*pi*me*kB*T)**1.5/hPlanck**3 *np.exp(-DeltaHeIIHeIII*eV/(kB*T))
	
	return 1.0 if b>1e12*Y else 1./(0.5*Y)*(-(b+X+Y/4.)+np.sqrt((b+X+Y/4.)*(b+X+Y/4.)+Y*b))


# =================================================================
#  X_i's (SCvH) ab x,y,z1,z2 || X_i's (DAB13) from x,y,z1,z2
# =================================================================

def XH2(x,y):
	if (y==1.):
		return 0.
	else:
		return 1./( 4*x*y/( (1-x)*(1-y)) + 2*y/(1-y) + 1 )

def XHI(x,y):
	if (y==0.):
		return 0.
	elif (x==1.):
		return 0.
	else:
		return 1./( 2*x/(1-x) + 1 + 0.5*(1-y)/y )

# -----------------------------------------------------------------

def XHeI  (z1,z2):
	if (z1==1. or z2==1.):
		return 0.
	elif (z1==0.):
		return 1.
	else:
		return 1./( 3.*z2/(1.-z2)*z1/(1.-z1) + 2.*z1/(1.-z1) + 1. )

def XHeII (z1,z2):
	if (z1==0.):
		return 0.
	elif (z2==1.):
		return 0.
	else:
		return 1./( 3.*z2/(1.-z2) + 2. + (1.-z1)/z1 )

def XHeIII(z1,z2):
	if (z2==0. or z1==0.):
		return 0.
	else:
		return ( 1./( 3. + 2*(1.-z2)/z2 + (1.-z1)/z1*(1.-z2)/z2 ) )

def XeHe  (z1,z2):
	return XHeI(z1,z2) + 2*XHeIII(z1,z2)


# =================================================================
#  X_i's und mu ab den anderen X_i's || X_i's and mu from the other X_i's
# =================================================================

def XeH_SCvH(XH2,XHI):
	x = 0.5*(1-XH2-XHI)
	minXe = 1e-5  # bei 1e-6 wird es rauschig für Wasserstoff, zumindest für logP = 4
	return 0. if x < minXe else x

def mu_Htot(XH2,XHI):
	'''Mittlere Molekularmasse für Wasserstoff || Mean molecular weight for hydrogen species'''
	return 0.5*(1+3*XH2+XHI)

# -----------------------------------------------------------------

def XHeIII_SCvH(XHe,XHeII):
	return 1/3.*(1-XHe-2*XHeII)

def XeHe_SCvH(XHeI,XHeII):
	x = 1/3.*(2-2*XHeI-XHeII)
	# man braucht kein minXe für Wasserstoff?
	return x

def mu_Hetot(XHeI,XHeII):
	'''Mittlere Molekularmasse für Helium || Mean molecular weight for helium species'''
	
	return 4./3.*(1+2*XHeI+XHeII)

# =================================================================
#  X_i's ab rho und T || X_i's from rho and T
# =================================================================
# 
# For the following functions the dependence in Y cancels out,
#  so it is arbitrary. The X_i's are defined for "only hydrogen"
#  or "only helium".
# Problems only for Y=1-eps, unsurprisingly
# TODO ZUTUN: Could one simplify the formulae?
# 

def XH2_rhoT (rho,T):
	Y = 0.5
	return XH2( xDAB13(Y,rho,T), yDAB13(Y,rho,T) )

def XHI_rhoT (rho,T):
	Y = 0.5
	return XHI( xDAB13(Y,rho,T), yDAB13(Y,rho,T) )

def XHII_rhoT(rho,T):
	Y = 0.5
	return XeH_rhoT(rho,T)

def XeH_rhoT (rho,T):
	return XeH_SCvH(XH2_rhoT(rho,T), XHI_rhoT(rho,T) )

# -----------------------------------------------------------------

def XHeI_rhoT  (rho,T):
	Y = 0.5
	return XHeI  ( z1DAB13(Y,rho,T), z2DAB13(Y,rho,T) )

def XHeII_rhoT (rho,T):
	Y = 0.5
	return XHeII ( z1DAB13(Y,rho,T), z2DAB13(Y,rho,T) )

def XHeIII_rhoT(rho,T):
	Y = 0.5
	return XHeIII( z1DAB13(Y,rho,T), z2DAB13(Y,rho,T) )

def XeHe_rhoT  (rho,T):
	return XHeII_rhoT(rho,T) + 2.*XHeIII_rhoT(rho,T)

# =================================================================
#  Entropien pro Teilchen für Wasserstoff oder Helium || Entropies per particle for hydrogen or helium
# =================================================================

# Entropie pro monatomisches Teilchen
def Smono_proTeilchen(rho,T,mu):
	'''Entropie pro Teilchen für monatomische Teilchen: Sackur--Tetrode-Gleichung'''
	
	return 2.5 + np.log( (2*pi* mu*mH* kB*T/hPlanck**2.)**1.5 / (rho/(mu*mH)) )

# Entropie pro diatomisches Teilchen
def Sdia_proTeilchen(rho,T,OP):
	'''Entropie pro diatomisches Teilchen (spezifisch: H2)
	     OP = "T.gg.Trot":  Ortho-zu-Para wenn T >> theta_rot
	     [einzige Möglichkeit für jetzt]
	'''
	if (OP == "T.gg.Trot"):
		Srot = np.log(T/(2.*thetarot_DAB13))
	else:
		Srot = 0.
	
	#  ZUTUN: Vibrationsteil nie sichtbar wichtig?
	Svib = np.log(1-np.exp(-thetavib_DAB13*1./T))
	
	return 1. + Smono_proTeilchen(rho,T,2.) + Srot - Svib


# =================================================================
#  Teilentropien pro Masseneinheit für Wasserstoff und Helium
#  || Partial entropies per unit mass for hydrogen and helium
# =================================================================

# -----------------------------------------------------------------
#   Wasserstoff || hydrogen
# -----------------------------------------------------------------

def sH2_proMasse_rhoT(rho,T,OP):
	'''Entropie pro Masseneinheit wenn nur H2 || Entropy per unit mass of only H2
	     OP = "T.gg.Trot":  Ortho-zu-Para wenn T >> theta_rot
	          [nur diese Modus ist zurzeit verfügbar]
	'''
	XH2 = XH2_rhoT(rho,T)
	XHI = XHI_rhoT(rho,T)
	return 1/mu_Htot(XH2,XHI)* XH2  * Sdia_proTeilchen(rho,T,OP)

def sHI_proMasse_rhoT(rho,T):
	XH2 = XH2_rhoT(rho,T)
	XHI = XHI_rhoT(rho,T)
	return 1/mu_Htot(XH2,XHI)* XHI  * Smono_proTeilchen(rho,T,1)

def sHII_proMasse_rhoT(rho,T):
	XH2 = XH2_rhoT(rho,T)
	XHI = XHI_rhoT(rho,T)
	XHII= XHII_rhoT(rho,T)
	return 1/mu_Htot(XH2,XHI)* XHII * Smono_proTeilchen(rho,T,1)

def seH_proMasse_rhoT (rho,T,K=0):
	'''Elektronenentropie pro Masseneinheit Wasserstoff.
	  K: integration constant for the electron entropy; best match to SCvH for K=0
	 -> entropy per mass in units of kB/mH (dimensionless)
	'''
	XH2 = XH2_rhoT(rho,T)
	XHI = XHI_rhoT(rho,T)
	XeH = XeH_rhoT(rho,T)
	if (XeH < 1e-5):
		SEl = 0.
	else:
		SEl = (Se_ideal_proTeilchen("nurH",rho,T) + K)
	
	return 1/mu_Htot(XH2,XHI)* XeH * SEl

# -----------------------------------------------------------------
#   Helium || helium
# -----------------------------------------------------------------

# Helium
def sHeI_proMasse_rhoT  (rho,T):
	XHeI  = XHeI_rhoT(rho,T)
	XHeII = XHeII_rhoT(rho,T)
	
	return 1/mu_Hetot(XHeI,XHeII)* XHeI   * Smono_proTeilchen(rho,T,4)

def sHeII_proMasse_rhoT (rho,T):
	XHeI  = XHeI_rhoT(rho,T)
	XHeII = XHeII_rhoT(rho,T)
	
	return 1/mu_Hetot(XHeI,XHeII)* XHeII  * Smono_proTeilchen(rho,T,4)

def sHeIII_proMasse_rhoT(rho,T):
	XHeI   = XHeI_rhoT(rho,T)
	XHeII  = XHeII_rhoT(rho,T)
	XHeIII = XHeIII_rhoT(rho,T)
	
	return 1/mu_Hetot(XHeI,XHeII)* XHeIII * Smono_proTeilchen(rho,T,4)

def seHe_proMasse_rhoT (rho,T,K=0):
	'''Elektronenentropie pro Masseneinheit Helium.
	  K: integration constant for the electron entropy; best match to SCvH for K=0
	 -> entropy per mass in units of kB/mH (dimensionless)
	'''
	XHeI  = XHeI_rhoT(rho,T)
	XHeII = XHeII_rhoT(rho,T)
	XeHe  = XeHe_rhoT(rho,T)
	
	return 1/mu_Hetot(XHeI,XHeII)* XeHe * (Se_ideal_proTeilchen("nurHe",rho,T) + K)



# =================================================================
#  Hilfsfunktionen und Entropie pro Elektron für Elektronen
#  || Auxiliary functions and entropies per electron for electrons
# =================================================================

def F_einhalb(ne,T):
	## ZUTUN hälfe hier etwas Umschreiben? -> Probleme bei kleinen T
	return ne / (4*pi/hPlanck**3.*(2*me*kB*T)**1.5)

def alphaIdeal(F_einhalb):
	# Fukushima (2015), App. Math. & Comp., 259, 698--707
	return np.log(F_einhalb/gamma(1.5))

# -----------------------------------------------------------------

# TODO ZUTUN man könnte vielleicht ne_H und ne_He mit einer Kombination
#      aus anderen Funktionen in dieser Datei ersetzen
def ne_H (rho,T):
	'''Elektronenanzahldichte für „Wasserstoffelektronen“ || Electron number density of "hydrogen electrons"
	'''
	XH2 = XH2_rhoT(rho,T)
	XHI = XHI_rhoT(rho,T)
	
	return rho/mH * ( 1 - XH2 - XHI )/( 1 + 3*XH2 + XHI )

def ne_He(rho,T):
	'''Elektronenanzahldichte für „Heliumelektronen“ || Electron number density of "helium electrons"
	'''
	XHeI  = XHeI_rhoT(rho,T)
	XHeII = XHeII_rhoT(rho,T)
	
	return rho/mH * 0.25 *( 2 - 2*XHeI - XHeII )/( 1 + 2*XHeI + XHeII )

# -----------------------------------------------------------------

def Se_ideal_proTeilchen(Modus,rho,T):
	'''Ideale Entropie (d.h. im Nichtentartetenlimes) pro Elektron (kB/Elektron)
	|| Ideal entropy (i.e. in the non-degenerate limit) per electron (kB/electron)
	Fukushima (2015), App. Math. & Comp., 259, 698--707
	'''
	if (Modus == "nurH"):
		ne = ne_H(rho,T)
	else:
		ne = ne_He(rho,T)
	
	return  0. if ne < 1e-5 else 15./4. - alphaIdeal(F_einhalb( ne, T ))



# =================================================================
# =================================================================
#  Hauptentropieformeln:
#    Gesamtentropien für Wasserstoff
#    Gesamtentropien für Helium
#    Mischungsentropie
#    Gesamtentropien für eine Mischung
# =================================================================
# =================================================================

def sHtot_proMasse_rhoT(rho,T,OP="T.gg.Trot", K=0,Protonspin=True):
	'''
	Total entropy per mass of hydrogen (in kB/mH)
	
	rho: density (g/cm^3)
	  T: temperature (K)
	 OP: mode for ortho-to-para ratio; only "T.gg.Trot" available for now
	  K: integration constant for the electron entropy; best match to SCvH for K=0
	 -> entropy per mass in units of kB/mH (dimensionless)
	
	Note: The published SCvH tables have a proton spin term
	but some versions of the code circulating do not include it.
	See Footnote 2 of Mordasini, Marleau & Mollière (2017)
	for the description of who has it and who not. It leads to
	an offset of Delta S = (1-Y)*ln(2) when combining with helium.
	'''
	return sH2_proMasse_rhoT(rho,T,OP) \
		+ sHI_proMasse_rhoT(rho,T) \
		+ sHII_proMasse_rhoT(rho,T) \
		+ seH_proMasse_rhoT (rho,T,K) \
		+ (np.log(2) if Protonspin else 0.)


def sHetot_proMasse_rhoT(rho,T,K=0):
	'''
	Total entropy per mass of helium (in kB/mH)
	  
	rho: density (g/cm^3)
	  T: temperature (K)
	  K: integration constant for the electron entropy; best match to SCvH for K=0
	 -> entropy per mass in units of kB/mH (dimensionless)
	'''
	return sHeI_proMasse_rhoT(rho,T) \
		+ sHeII_proMasse_rhoT(rho,T) \
		+ sHeIII_proMasse_rhoT(rho,T) \
		+ seHe_proMasse_rhoT (rho,T,K)

# -----------------------------------------------------------------

# siehe Korrektur in Fußnote 4 im Anhang von Baraffe et al. (2008)
def sMisch_proMasse(Y,XH2,XHI,XeH, XHeI,XHeII,XeHe):
	'''
	Ideal mixing entropy from Eqs. (53--56) in SCvH but with correction of typo
	   according to Footnote 4 (Appendix) of Baraffe et al. (2008).
	   Output in kB/mH'''
	
	betagamma = 1./4.*Y/(1-Y) * 1.5*(1+XHI+3*XH2)/(1.+2*XHeI+XHeII)
	delta     = 2./3.*(2-2*XHeI-XHeII)/(1.-XH2-XHI) * betagamma
	
	logdTerm = 0. if XeH==0 else np.log(1+delta)
	loginvdTerm = 0. if XeHe==0 else np.log(1+1./delta)
	
	return (1.-Y)*2./(1+XHI+3*XH2) * ( np.log(1+betagamma) - XeH * logdTerm \
		+ betagamma* ( np.log(1+1./betagamma) - XeHe * loginvdTerm ) )


def sMisch_proMasse_rhoT(Y,rho,T):
	'''Ideal mixing entropy from Eqs. (53--56) in SCvH but with correction of typo
	   according to Footnote 4 (Appendix) of Baraffe et al. (2008).
	   Output in kB/mH'''
	
	XH2   = XH2_rhoT(rho,T)
	XHI   = XHI_rhoT(rho,T)
	XeH   = XeH_rhoT(rho,T)
	XHeI  = XHeI_rhoT(rho,T)
	XHeII = XHeII_rhoT(rho,T)
	XeHe  = XeHe_rhoT(rho,T)
	
	return sMisch_proMasse(Y,XH2,XHI,XeH, XHeI,XHeII,XeHe)

# -----------------------------------------------------------------

# 
# die Gesamtentropie der Mischung pro Masseneinheit, in kB/mH
#  als Funktion von (P,T)
# 
def stot_proMasse_PT(Y, P, T, OP="T.gg.Trot", K=0,Protonspin=True):
	'''Gesamtentropie pro Gramm einer H+He-Mischung || Total entropy per unit mass of an H+He mixture
	Changes (P,T) to (rho,T) and calls stot_proMasse_rhoT()
	-> see its documentation
	'''
	rho = Dichte_PT(Y,P,T)
	return stot_proMasse_rhoT(Y, rho, T, OP=OP, K=K,Protonspin=Protonspin)

def stot_proMasse_PT_bystro(Y, P, T, OP="T.gg.Trot", K=0,Protonspin=True):
	'''Gesamtentropie pro Gramm einer H+He-Mischung || Total entropy per unit mass of an H+He mixture
	Changes (P,T) to (rho,T) and calls stot_proMasse_rhoT_bystro()
	-> see its documentation
	'''
	rho = Dichte_PT(Y,P,T)
	return stot_proMasse_rhoT_bystro(Y, rho, T, OP=OP, K=K,Protonspin=Protonspin)


# 
# die Gesamtentropie der Mischung pro Masseneinheit, in kB/mH
#  als Funktion von (rho,T)
#
def stot_proMasse_rhoT(Y, rho, T, OP="T.gg.Trot", K=0,Protonspin=True):
	'''Gesamtentropie pro Gramm einer H+He-Mischung || Total entropy per unit mass of an H+He mixture
	  Y: helium mass fraction
	rho: density (g/cm^3)
	  T: temperature (K)
	 OP: mode for ortho-to-para ratio; only "T.gg.Trot" available for now
	  K: integration constant for the electron entropy; best match to SCvH for K=0
	 -> entropy per mass in units of kB/mH (dimensionless)
	'''
	return (1.-Y)*sHtot_proMasse_rhoT(rho,T,OP,K,Protonspin=True) \
		+ Y*sHetot_proMasse_rhoT(rho,T,K) \
		+ sMisch_proMasse_rhoT(Y,rho,T)



# wie stot_proMasse_rhoT aber schneller -- Formeln sonst gleich
def stot_proMasse_rhoT_bystro(Y, rho, T, OP="T.gg.Trot", K=0,Protonspin=True):
	'''Gesamtentropie pro Gramm einer H+He-Mischung || Total entropy per unit mass of an H+He mixture
	   * Schnellere Ausrechnung als mit stot_proMasse_rhoT weil Unterfunktionen hineinkopiert
	     || Faster execution than with stot_proMasse_rhoT because subroutines copied into this function
	  
	  Y: helium mass fraction
	rho: density (g/cm^3)
	  T: temperature (K)
	 OP: mode for ortho-to-para ratio; only "T.gg.Trot" available for now
	  K: integration constant for the electron entropy; best match to SCvH for K=0
	 -> entropy per mass in units of kB/mH (dimensionless)
	'''
	
	XH2    = XH2_rhoT(rho,T)
	XHI    = XHI_rhoT(rho,T)
	XHII   = XHII_rhoT(rho,T)
	XeH    = XeH_rhoT(rho,T)
	XHeI   = XHeI_rhoT(rho,T)
	XHeII  = XHeII_rhoT(rho,T)
	XHeIII = XHeIII_rhoT(rho,T)
	XeHe   = XeHe_rhoT(rho,T)
	
	# -------------------------------------------------------------
	# -------------------------------------------------------------
	
	Smono_proTeilchen_1 = Smono_proTeilchen(rho,T,1)
	
	mu_Htot_ = mu_Htot(XH2,XHI)
	
	sH2_proMasse_rhoT  = 1/mu_Htot_* XH2  * Sdia_proTeilchen(rho,T,OP)
	
	sHI_proMasse_rhoT  = 1/mu_Htot_* XHI  * Smono_proTeilchen_1
	
	sHII_proMasse_rhoT = 1/mu_Htot_* XHII * Smono_proTeilchen_1
	
	if (XeH < 1e-5):
		SeH = 0.
	else:
		SeH = (Se_ideal_proTeilchen("nurH",rho,T) + K)
	
	seH_proMasse_rhoT = 1/mu_Htot_ * XeH * SeH
	
	sHtot_proMasse_rhoT = sH2_proMasse_rhoT \
		+ sHI_proMasse_rhoT \
		+ sHII_proMasse_rhoT \
		+ seH_proMasse_rhoT \
		+ (np.log(2) if Protonspin else 0.)
	
	# -------------------------------------------------------------
	
	Smono_proTeilchen_4 = Smono_proTeilchen(rho,T,4)
	
	mu_Hetot_ = mu_Hetot(XHeI,XHeII)
	
	sHeI_proMasse_rhoT   = 1/mu_Hetot_ * XHeI   * Smono_proTeilchen_4
	
	sHeII_proMasse_rhoT  = 1/mu_Hetot_ * XHeII  * Smono_proTeilchen_4
	
	sHeIII_proMasse_rhoT = 1/mu_Hetot_ * XHeIII * Smono_proTeilchen_4
	
	seHe_proMasse_rhoT   = 1/mu_Hetot_ * XeHe * (Se_ideal_proTeilchen("nurHe",rho,T) + K)
	
	sHetot_proMasse_rhoT = sHeI_proMasse_rhoT \
		+ sHeII_proMasse_rhoT \
		+ sHeIII_proMasse_rhoT \
		+ seHe_proMasse_rhoT
	
	# -------------------------------------------------------------
	# -------------------------------------------------------------
	
	betagamma = 1./4.*Y/(1-Y) * 1.5*(1+XHI+3*XH2)/(1.+2*XHeI+XHeII)
	delta     = 2./3.*(2-2*XHeI-XHeII)/(1.-XH2-XHI) * betagamma
	
	logdTerm = 0. if XeH==0 else np.log(1+delta)
	loginvdTerm = 0. if XeHe==0 else np.log(1+1./delta)
	
	sMisch_proMasse_rhoT = (1.-Y)*2./(1+XHI+3*XH2) * ( np.log(1+betagamma) - XeH * logdTerm \
		+ betagamma* ( np.log(1+1./betagamma) - XeHe * loginvdTerm ) )
	
	# -------------------------------------------------------------
	# -------------------------------------------------------------
	
	return (1.-Y)*sHtot_proMasse_rhoT \
		+ Y*sHetot_proMasse_rhoT \
		+ sMisch_proMasse_rhoT



# =================================================================
#  Zusätzliche Funktionen
# =================================================================

def muDAB13(Y,rho,T,Z=0):
	'''Mittlere Molekularmasse für eine H+He-Mischung mit den D'Angelo & Bodenheimer (2013)-Formeln
	|| Mean molecular weight with the D'Angelo & Bodenheimer (2013) formulae'''
	if (Y<0 or Y>1):
		raise Exception('[muDAB13] Y<0 oder >0: {}'.format(Y))
	x  = xDAB13(Y,rho,T,Z=Z)
	y  = yDAB13(Y,rho,T,Z=Z)
	z1 = z1DAB13(Y,rho,T,Z=Z)
	z2 = z2DAB13(Y,rho,T,Z=Z)
	
	return 4.0/(2*(1.-Y)*(1.+y+2*x*y)+Y*(1.+z1+z1*z2))

def Druck_rhoT(Y,rho,T,Z=0):
	return rho/muDAB13(Y,rho,T,Z=Z) /mH *kB*T

# Achtung: Reihenfolge anders weil rho das 1. sein muß
def DeltaDruck(rho, Y,T,PZiel):
	P = Druck_rhoT(Y,rho,T)
	
	#lg = np.log(P/PZiel)	
	#x = (lg**2. + (P - PZiel)**2.)
	#print '  > ',PZiel, P, rho, x, lg
	##return x if lg > 0 else -x
	##return lg
	
	# ist doch am besten:
	return P - PZiel

# von https://github.com/andrewcumming/gasgiant (David Berardo) übernommen und (sehr) angepaßt
def Dichte_PT(Y,P,T):
	'''Dichte aus (P,T) mit den DAB13-Funktionen invertieren || Invert density from (P,T) with the DAB13 functions
	Wenn keine Lösung: rho = 0 || If no solution is found, return 0
	'''
	# großzügige Grenzen || generous boundaries
	rhomin = 1e-40
	rhomax = 1e2
	eps = 1e-8
	# WICHTIG!:
	#  der Nullstellenfinder benutzt
	#    np.allclose(x, x0, atol=xtol, rtol=rtol)
	#  als Kriterium aber der Standardwert für xtol ist 1e-5...
	#    http://lagrange.univ-lyon1.fr/docs/numpy/1.11.0/reference/generated/numpy.allclose.html
	# -> das geht also unter Umständen nicht, wenn rho klein sein soll!
	#   Da atol=0 nicht erlaubt ist, benutzen wir 1e-308
	atol = 1e-308
	
	f1 = DeltaDruck(rhomin,Y,T,P)
	f2 = DeltaDruck(rhomax,Y,T,P)
	
	if (f1*f2 > 0.0):
		print("[Dichte_PT] keine Lösung -> rho = 0 || no solution") # T,P,f1,f2)
		return 0.0
	rho = brentq(DeltaDruck,rhomin,rhomax,xtol=atol,rtol=eps,args=(Y,T,P))
	
	## überprüfen, ob es gut übereinstimmt:
	#Pcheck = Druck_rhoT(Y,rho,T)   #rho*kB*T/muDAB13(Y,rho,T)/mH
	#print '*  ', P, Pcheck, rho, P/Pcheck
	
	return rho


def nHDAB13(Y,rho,T,Z=0):
	'''number density of atomic hydrogen'''
	
	XH2 = XH2_rhoT(rho,T)
	XHI = XHI_rhoT(rho,T)
	
	# ACHTUNG: in Gl. 36 und 37 von Saumon, Chabrier & van Horn (1995)
	#          wird so getan, als ob das rho nur von H bzw. He komme
	#          -> Faktor X=1-Y-Z
	return 2*rho*(1.-Y-Z)/amu / ( 1 + 3.*XH2 + XHI ) * XHI


# =================================================================
#  delad berechnen um effektiv die Entropie mit DAB13 vergleichen zu können || compute delad to compare with DAB13
# =================================================================

def deladDAB13_rhoT(Y,rho,T, OP="T.gg.Trot", K=0,Protonspin=True):
	'''Berechnet delad = (dlnT/dlnP)_{konst. S} mit den DAB13-Funktionen
	ohne Invertierung obwohl S eine Funktion von rho und T ist (aber dafür mit vielen Ableitungen)
	|| computes delad from DAB13 without derivatives'''
	
	eps = 1e-2  # viel kleiner wird rauschig und ist nicht nötig
	
	# Schreibweise: QX_Y \equiv (dQ/dX)_{konst. Y}
	
	S0 = stot_proMasse_rhoT_bystro(Y, rho,         T,         OP=OP,K=K,Protonspin=Protonspin)
	S1 = stot_proMasse_rhoT_bystro(Y, rho*(1+eps), T,         OP=OP,K=K,Protonspin=Protonspin)
	S2 = stot_proMasse_rhoT_bystro(Y, rho,         T*(1+eps), OP=OP,K=K,Protonspin=Protonspin)
	
	Srho_T = (S1-S0)/(eps*rho)
	ST_rho = (S2-S0)/(eps*T)
	
	P0 = Druck_rhoT(Y, rho,         T        )
	P1 = Druck_rhoT(Y, rho*(1+eps), T        )
	P2 = Druck_rhoT(Y, rho,         T*(1+eps))
	
	Prho_T = (P1-P0)/(eps*rho)
	PT_rho = (P2-P0)/(eps*T)
	
	TP_S = - Srho_T / (ST_rho*Prho_T - Srho_T*PT_rho)
	
	return P0/T * TP_S


# =================================================================
#  Daten ausgeben
# =================================================================

def ZGTabellespeichern(Y, OP="T.gg.Trot", K=0,Protonspin=True):
	
	#// Achtung: wir schreiben Y0.750 usw. aber meinen eigentlich X0.750...
	#//   der Rückwätskompatibilität wegen so lassen
	#sprintf(Name, "Modelldaten/ZG_Y%.3f%s%s%s%s_230x800.dat",
		#H_MASS_FRAC,
		#(HELIUM_IONIZATION == YES? "": "_kHeIonis"),
		#(ORTHO_PARA_MODE == 0? "_nurParaH": (ORTHO_PARA_MODE == 1? "_OPGg": "_OP31")),
		#(PV_TEMPERATURE_TABLE == YES? "_TempTabJA":""),
		#(TV_ENERGY_TABLE == YES? "_EnTabJA":"")
	#);
	Dateiname = 'ZGTabelle_DAB13_mitEntropie.dat'
	
	nrho=231
	nt=801
	
	logTmin = 1.
	logTmax = 5.
	logrhomin = -23.
	logrhomax = 0.
	
	# ACHTUNG: viele Größen werden nicht ausgerechnet und sind deswegen Null
	with open(Dateiname,'w') as fd:
		# Kopfzeile wie in /home/gabriel/Dokumente/Rolf_2013-09-17/src/ModifiedPluto/main.c
		fd.write("#-------------1:T       2:rho         3:P          4:mu     5:Gamm1        6:cv(kB/mH)7:eint(erg/g)  8:gammaeff                 9:x            10:y           11:z1           12:z2      13:delad   14:cp(kB/mH)  15:dlnrhodlnT_P_rhoT   16:chiT   17:Gamm3   18:Entropie(kB/mH)\n")
		for j in range(nrho):
			
			rho = 10.**(logrhomin + j*(logrhomax-logrhomin)/1./(nrho-1))
			
			for i in range(nt):
				
				T = 10.**(logTmin + i*(logTmax-logTmin)/1./nt)
				P = Druck_rhoT(Y,rho,T)

				mu = muDAB13(Y,rho,T)
				Gamm1 = 0
				cv = 0
				eint = 0
				gammaeff = 0
				x = xDAB13(Y,rho,T)
				y = yDAB13(Y,rho,T)
				z1 = z1DAB13(Y,rho,T)
				z2 = z2DAB13(Y,rho,T)
				del_ad = deladDAB13_rhoT(Y,rho,T, OP=OP, K=K,Protonspin=Protonspin)
				cp = 0
				dlnrhodlnT_P = 0  # könnte man schon ausgeben aber wird nicht benötigt...
				chiT = 0
				Gamm3 = 0
				
				stot = stot_proMasse_rhoT(Y,rho,T,OP=OP,K=K,Protonspin=Protonspin)
				
				fd.write(" %15.4e  %15.4e  %10.4e    %10.4f  %10.4f  %10.4f  %15.4e  %10.4f      %14.3e  %14.3e  %14.3e  %14.3e    %10.4f  %15.3e   %15.4e   %10.4g   %10.4f  %10.4f  \n" \
				  %(T, rho, P, mu,Gamm1,cv, eint, gammaeff, x,y,z1,z2, del_ad, cp, dlnrhodlnT_P, chiT, Gamm3, stot))
				
				if (i%100==0 and j%10==0): print(i, rho, T, stot)
			
			fd.write("\n")
	
	print(" > "+Dateiname+" gemacht")



# =================================================================
#  Konstanten
# =================================================================

pi = 3.14159265358979    # Torte
me = 9.11e-28            # g
amu = 1.661e-24          # g; atomic mass unit
mH = 1.673e-24           # g
kB = 1.381e-16           # erg/K
eV = 1.602e-12           # erg
c = 3e10                 # cm/s
e = 4.80e-10             # e.s.u. = dyn^0.5 cm
G = 6.67e-8              # cgs
RJ = 7.15e9              # cm
MJ = 1.898e30            # g
ME = 5.9736e27           # g; PLUTO (CONST_Mearth)
annum = 0.9461e18/ 2.99792458e10  # s; PLUTO (CONST_ly/CONST_c)
nmtocm = 1e-7            # conversion from nm to cm
hPlanck = 6.626e-27      # Planck's constant
thetarot_DAB13 = 85.5    # Rotationstemperatur von Wasserstoff in D'Angelo & Bodenheimer (2013)
thetavib_DAB13 = 6140.   # Vibrationstemperatur " " " "

