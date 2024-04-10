import math as m
import numpy as np

def FDeltaGV(TempK):
    """" Calculates the driving force for nucleation in J/m^3
    Based on a Thermocalc derived curve at different temperatures (verify against Porter and Easterling)
    This is different than what Makiewicz used and is the true nucleation driving force from ThermoCalc (shared tangent after Winterbottom)
    Fits are 6th order polynomials, used two different T ranges for a better fit
    Beta transus is +/- 1263 K for straight Ti64, but this is affected by other elements
    Currently for only 1 composition Ti6Al4V, no contaminants """

    if TempK < 650:

        a6 = -9.129099800E-08
        a5 = 2.190555350E-04
        a4 = -2.121865270E-01
        a3 = 1.055251740E+02
        a2 = -2.822220780E+04
        a1 = 4.339870240E+06
        a0 = -7.609064690E+08
        return a6 * (TempK ** 6) + a5 * (TempK ** 5) + a4 * (TempK ** 4) + a3 * (TempK ** 3) + a2 * (
                TempK ** 2) + a1 * (TempK ** 1) + a0

    elif 650 <= TempK <= 1400:

        a6 = 3.944502580E-11
        a5 = -2.498349840E-07
        a4 = 6.663610880E-04
        a3 = -9.520011940E-01
        a2 = 5.581605910E+02
        a1 = 4.445091590E+05
        a0 = -5.864387420E+08
        return a6 * (TempK ** 6) + a5 * (TempK ** 5) + a4 * (TempK ** 4) + a3 * (TempK ** 3) + a2 * (
                    TempK ** 2) + a1 * (TempK ** 1) + a0

    else:

        return 0

def NucleationRate(TempK, NZero, gammaAB, deltaGS, theta):
    """ Calculate heterogeneous nucleation rate per m^3, or simplifies to homogeneous at theta=180
    TempK is the temperature at which nucleation is occurring in this time step in K
    NZero is the number of nucleation sites available per unit volume and time in 1/(m^3*s)
    gammaAB is surface/interfacial energy between the Alpha and Beta phases per unit area in J/m^2
    deltaGV is the free energy change due to the phase transformation per unit volume in J/m^3
    deltaGS is the strain energy change due to the phase transformation per unit volume in J/m^3
    theta is the phase contact angle between Alpha and Beta ("droplet" of Alpha on Beta) in degrees"""

    #Constants:
    k = 1.38064852E-23                  # Boltzman's constant in J/K
    R = 8.3144598                        # Gas Constant  J/(K*Mole)

    #Calculate effect of contact angle on activation energy - from Porter and Easterling 2nd Edition pp. 272
    STheta = ((2 + m.cos(m.radians(theta)))*(1 - m.cos(m.radians(theta)))**2)/2

    #Calculate Gibbs Free Energy change due to phase change
    deltaGV = FDeltaGV(TempK)

    #Calculate Activation Energy
    if deltaGV + deltaGS < 0:

        QHetero = 16*m.pi*gammaAB**3/(3*(deltaGV + deltaGS)**2)*STheta
        #Calculate Nucleation Rate
        #Account for Machine Epsilon
        if NZero * m.exp(-QHetero / (k * (TempK))) == 0 :

            return 1E-10

        else:

            return NZero * m.exp(-QHetero / (k * (TempK)))

    else:

        return 1E-10


def tEquivalent(Vol,NucRate,GrowthRate,JMAKexp):
    """Calculate equivalent time as defined by SM Kelly, 2004
    Vol is volume from the previous time step
    This calculates the time it would've taken to from an input extent using the kinetic function defined
    at the temperature implicit to the nucleation and growth rate input"""

    if (NucRate*8*m.pi/15*GrowthRate**3) <= 0:

        return 0

    else:

        return (Vol/(NucRate*8*m.pi/15*GrowthRate**3))**(1/JMAKexp)

def DVinBeta(TempK):
    """Calculates the diffusivity of V in Beta phase as a function of temperature in m^2/s
    TempK is the temperature of the process in Kelvin
    This is based on Semiatin et al. "Diffusion Coefficients for Modeling the Heat Treatment of Ti-6Al-4V" 2004"""

    D0V = 1E5        #Diffusion constant
    QDVk = 17460     #Activation Energy Divided by Boltzmann's constant

    return D0V * m.exp(-QDVk/TempK)*1E-12            #initially in micrometers^2/2, converted by 1E-12 to m^2/s

def TZero(Van):
    """TZero is Beta to Alpha spontaneous transition temperature at a V content of Van in Kelvin
    This is where the free energies are equal - need to reproduce this in ThermoCalc
    Van is Vanadium content in mole fraction
    This is not yet used"""

    if Van <= 0.28601:

        return 1274.7 - 4548.4*Van + 8174*Van**2 - 11316*Van**3

def VEqConc(TempK,Phase):
    """Calculates fractional molar Equilibrium Vanadium Composition at TempK in the phase
    TempK is the temperature in Kelvin
    Functions are from Thermocalc - didn't bother modeling Beta below 400 K, previous function is about 0.95 at that T
    Significant differences from Makiewicz results, particularly in Alpha
    Phase is which phase to calculate the concentration in, Alpha, or Beta"""

    if Phase == "Alpha":

        a6 = 7.979842E-19
        a5 = -3.416356E-15
        a4 = 5.528480E-12
        a3 = -4.337816E-09
        a2 = 1.809364E-06
        a1 = -3.867381E-04
        a0 = 3.314268E-02
        return a6 * (TempK ** 6) + a5 * (TempK ** 5) + a4 * (TempK ** 4) + a3 * (TempK ** 3) + a2 * (
                TempK ** 2) + a1 * (TempK ** 1) + a0

    elif Phase == "Beta":

        if TempK >= 1262.4734:

            return 0.0360096

        elif 1262.4734 >= TempK >= 400:

            a6 = -3.602880E-17
            a5 = 1.882117E-13
            a4 = -3.965737E-10
            a3 = 4.306646E-07
            a2 = -2.538350E-04
            a1 = 7.577572E-02
            a0 = -7.943998
            return a6 * (TempK ** 6) + a5 * (TempK ** 5) + a4 * (TempK ** 4) + a3 * (TempK ** 3) + a2 * (
                    TempK ** 2) + a1 * (TempK ** 1) + a0

        elif TempK < 400:

            return 0.95

        else:

            return 0

    else:

        return 0


def GrowthRate(TempK):
    """Calculates growth rate of Alpha at TempK based on Zener approach (Porter and Easterling 2nd Ed. pp 280)
    Resulting units are m/s^(1/2)
    X0 is alloy average composition of Vanadium (i.e. 0.04 for Ti64)
    teq is equivalent time for location on the sigmoidal curve"""

    XVMolar = 0.0360096  # Average V alloy composition in molar concentration
    XAlpha = VEqConc(TempK, "Alpha")  # X_Beta in Zener equation, here it is V Conc. in Alpha
    XBeta = VEqConc(TempK, "Beta")  # X_e in Zener Equation, here it is V Conc. in Beta
    Diff = DVinBeta(TempK) #Diffusivity in the Beta phase, assumed to be kinetically limiting

    if (XVMolar - XBeta) / ((XAlpha - XBeta)) * m.sqrt(Diff) < 0:

        return 1E-30

    else:

        return (XVMolar - XBeta) / ((XAlpha - XBeta)) * m.sqrt(Diff)

def EquilibriumAlpha(TempK):
    """Calculate Equilibrium Alpha volume fraction based on ThermoCalc for Ti64
    TempK is the temperature in Kelvin
    Curve is adjusted to reflect the ~91% alpha fraction typically seen. TC yields more like 97%,
    with the Beta phase composed of nearly 100% V, not seen in reality. This is done by a simple scaling factor
    outside of this module"""

    if TempK >= 1262.4734:      #Above the Beta Transus for straight Ti64

        return 0

    elif 1002 <= TempK <= 1262.4734:     #First curve fit for volumetric Alpha fraction

        a6 = 2.721442922E-15
        a5 = -2.077090893E-11
        a4 = 6.450547205E-8
        a3 = -1.049698524E-4
        a2 = 9.477852654E-2
        a1 = -45.14394796
        a0 = 8.880180232E3
        return a6 * (TempK ** 6) + a5 * (TempK ** 5) + a4 * (TempK ** 4) + a3 * (TempK ** 3) + a2 * (
                    TempK ** 2) + a1 * (TempK ** 1) + a0

    else:                       #Second curve fit for volumetric Alpha fraction

        a6 = 0
        a5 = 0
        a4 = 0
        a3 = 0
        a2 = 0
        a1 = 0
        a0 = 0.91
        return a6 * (TempK ** 6) + a5 * (TempK ** 5) + a4 * (TempK ** 4) + a3 * (TempK ** 3) + a2 * (
                    TempK ** 2) + a1 * (TempK ** 1) + a0

def MartensiteFraction(Ms,TempK,MartBPara,BetaFrac):
    """Based on Magee theory / Koistenin-Marburger - 1959 - "A general equation prescribing the extent of
    austenite-martnesite transformation"
    TempK is current temperature in Kelvin
    Ms in Martensite Start Temperature in Kelvin
    BetaFrac is the remaining fraction of Beta at the timestep"""

    if 1-m.exp(-MartBPara*(Ms-TempK)) < 0:

        return 0

    else:

        return BetaFrac*(1-m.exp(-MartBPara*(Ms-TempK)))


def CalculateExtentAlpha(IZero,GrowthRate,time,JMAKexp):
    """Calculates the extent of beta to alpha tranformation in a time step
    IZero is the nucleation rate per m^3
    GrowthRate is the growth rate in m/(s^(1/2))
    EqAlpha is the amount of Alpha at equilibrium in the time step (at isotherm)
    time is the equivalent time in seconds"""


    #Calculate Real Alpha Volume
    VReal = IZero*(8*m.pi/15)*(GrowthRate**3)*(time**(JMAKexp))

    #Calculate Extent this represents
    if 1-m.exp(-VReal) >= 1:

        return 1

    else:

        return 1-m.exp(-VReal)

def LathThickness(tnlath,FracWid,deltafracwid,TempK):
    """Calculates the lath thickness formed at a given temperature and averages it with previous
    calculated average lath thickness (from the previous steps).
    Originally established by Charles (Charles - 2008 - Dissertation - Modeling Microstructure
    Evolution of Weld Deposited Ti-6Al-4V), she explained more later (Charles-Murgau - 2016-PhD_thesis-
    Lulea-Murgau-microstructure-model-Ti64-AM)
    This is designed to work only for either the colony or the basketweave morphology separately.
    This is a result of their different initiation physics.
    tnlath is the average calculated lath thickness from the previous steps
    FracWid is the fraction of the morphology already existing
    deltafracwid is the fraction formed in the current timestep
    TempK is the temperature in Kelvin of the timestep"""

    klath = 5                   # Arrhenius Prefactor
    Rlath = 2300                  # Activation Temperature
    teqlath = klath * m.exp(-Rlath / TempK)
    """ Equilibrium calculated lath thickness - 
    based on theory that higher undercooling leads to more individual nucleation events / lamellar direction growth """


    return (tnlath*FracWid + teqlath*deltafracwid) / (FracWid + deltafracwid)
