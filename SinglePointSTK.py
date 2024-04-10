import numpy as np
from scipy.signal import find_peaks
import NGFunctions as NG



# In[3]:

def ThermalStep(Timesa, TempCa):

    TempKa = [x + 273.15 for x in TempCa]

    imax = len(Timesa)

    MartStart = 848             # Martensite start temperature in Kelvin
    TBetaTransus = 1262.4734  # Beta Transus Temperature in Kelvin
    TDiss = 1252.4734  # Defines lower temperature peak bound for layer band formation
    MartBPara = 3E-3  # Murgau-Charles 2016 Thesis pp. 49 from Elmer et al. 2005b
    NZeroGB = 5E16  # sites / m^3
    NZeroBW = 3E20  # sites / m^3
    NZeroCOL = 2E17  # sites / m^3
    NZeroCOLSpecial = 4E19  # sites / m^3 for peaks below TBetaTransus and above TDiss
    gammaAB = 0.01          # J / m^2
    deltaGSGB = 0           # J / m^3
    deltaGSBW = 45E6         # J / m^3
    deltaGSCOL = 0.5E6     # J / m^3
    thetaGB = 5
    thetaBW = 5
    thetaCOL = 5
    JMAKexpBW = 5/2
    JMAKexpCOL = 3/2

    fractionGB = [0.01]        #initialize fraction arrays
    fractionBW = [0.8]
    fractionCOL = [0.1]
    fractionMART = [0]
    fractionMASS = [0]
    alpha = fractionGB[0] + fractionBW[0] + fractionCOL[0] + fractionMART[0] + fractionMASS[0]
    fractionalpha = [alpha]
    tlath = [0.5]

    x = np.asarray(TempKa)
    x = np.ravel(x)
    peaks, value = find_peaks(x, height=1, threshold=None, distance=5)
    heights = value['peak_heights']
    peaks_special = np.logical_and(heights >= TDiss, heights < TBetaTransus)
    if np.any(peaks_special):

        peaks_max = np.argmax(heights[peaks_special])
        original_idx = (peaks[peaks_special])[peaks_max]

    else:

        original_idx = None

    for i in range(0,imax-1):                   # Loop with thermal data and times to calculate microstructural development

        if (TempKa[i] + TempKa[i+1]) / 2  <= 0:

            TempK = 0.1

        else:
            TempK = (TempKa[i] + TempKa[i+1]) / 2
        TempStepK = TempKa[i+1] - TempKa[i]         # Calculate deltaT
        timestep = Timesa[i+1] - Timesa[i]          # Start with 2nd point for timestep calculation

        EqAlpha = NG.EquilibriumAlpha(TempK)
        beta = 1 - fractionalpha[i]
        fractionbeta = [beta]

        if EqAlpha > fractionalpha[i]:          #Alpha nucleaction and growth calculations

            if fractionGB[i] < EqAlpha and NG.NucleationRate(TempK, NZeroGB, gammaAB, deltaGSGB, thetaGB) > 0:     #GB Alpha calculations

                NucRateGB = NG.NucleationRate(TempK, NZeroGB, gammaAB, deltaGSGB, thetaGB)
                GRateGB = NG.GrowthRate(TempK)
                VoleGB = -1*np.log(1 - fractionGB[i]/EqAlpha)
                EquivTimeGB = NG.tEquivalent(VoleGB,NucRateGB,GRateGB,JMAKexpCOL)
                TimeGB = EquivTimeGB + timestep
                fractionGBE = EqAlpha*NG.CalculateExtentAlpha(NucRateGB,GRateGB,TimeGB,JMAKexpCOL)
                deltafracGBE = fractionGBE - fractionGB[i]

            else:

                deltafracGBE = 0

            if fractionBW[i] < EqAlpha and NG.NucleationRate(TempK, NZeroBW, gammaAB, deltaGSBW, thetaBW) > 0:    #BW Alpha calculations

                NucRateBW = NG.NucleationRate(TempK, NZeroBW, gammaAB, deltaGSBW, thetaBW)
                GRateBW = NG.GrowthRate(TempK)
                VoleBW = -1*np.log(1 - fractionBW[i]/EqAlpha)
                EquivTimeBW = NG.tEquivalent(VoleBW,NucRateBW,GRateBW,JMAKexpBW)
                TimeBW = EquivTimeBW + timestep
                fractionBWE = EqAlpha*NG.CalculateExtentAlpha(NucRateBW,GRateBW,TimeBW,JMAKexpBW)
                deltafracBWE = fractionBWE - fractionBW[i]

            else:

                deltafracBWE = 0

            if fractionCOL[i] < EqAlpha and NG.NucleationRate(TempK, NZeroCOL, gammaAB, deltaGSCOL,
                                                           thetaCOL) > 0:  # COL Alpha calculations

                if original_idx == None:

                    NucRateCOL = NG.NucleationRate(TempK, NZeroCOL, gammaAB, deltaGSCOL, thetaCOL)

                else:

                    if i < original_idx:

                        NucRateCOL = NG.NucleationRate(TempK, NZeroCOL, gammaAB, deltaGSCOL, thetaCOL)

                    else:

                        NucRateCOL = NG.NucleationRate(TempK, NZeroCOLSpecial, gammaAB, deltaGSCOL, thetaCOL)

                GRateCOL = NG.GrowthRate(TempK)
                VoleCOL = -1 * np.log(1 - fractionCOL[i] / EqAlpha)
                EquivTimeCOL = NG.tEquivalent(VoleCOL, NucRateCOL, GRateCOL, JMAKexpCOL)
                TimeCOL = EquivTimeCOL + timestep
                fractionCOLE = EqAlpha * NG.CalculateExtentAlpha(NucRateCOL, GRateCOL, TimeCOL, JMAKexpCOL)
                deltafracCOLE = fractionCOLE - fractionCOL[i]

            else:

                deltafracCOLE = 0

            if TempK >= MartStart:  # Need to be below Martensite start temperature for diffusionless transformation

                deltafracMARTE = 0
                deltafracMASSE = 0

            else:

                if NG.MartensiteFraction(MartStart,TempK,MartBPara,beta) <= 1:

                    deltafracMARTE = NG.MartensiteFraction(MartStart, TempK, MartBPara, beta)

                else:

                    deltafracMARTE = 1

                fractionMASSE = 0  # Need to look at Charles
                deltafracMASSE = 0

            """Use JMAK correction adapted for multiple competing morphologies as detailed in Jones and Bhadesia, 1997
            Note: Instead of (1-X_transformed), it makes more sense to use (EqAlpha-X_Transformed) since this transformation
            is a precipitation one and occurs into a two phase field, unlike the transformation from FCC to BCC described 
            by Jones and Bhadesia, which is eutectic, with the total equilibrium transformed fraction being 1 once the 
            temperature falls below the transformation temperature."""

            deltafracGB = (1 - (fractionalpha[i]) / EqAlpha) * deltafracGBE
            deltafracBW = (1 - (fractionalpha[i]) / EqAlpha) * deltafracBWE
            deltafracCOL = (1 - (fractionalpha[i]) / EqAlpha) * deltafracCOLE
            deltafracMART = (1 - (fractionalpha[i]) / EqAlpha) * deltafracMARTE
            deltafracMASS = (1 - (fractionalpha[i]) / EqAlpha) * deltafracMASSE

            # Append the new fraction to the list for each morphology
            fractionGB.append(fractionGB[i] + deltafracGB)
            fractionBW.append(fractionBW[i] + deltafracBW)
            fractionCOL.append(fractionCOL[i] + deltafracCOL)
            fractionMASS.append(fractionMASS[i] + deltafracMASS)
            fractionMART.append(fractionMART[i] + deltafracMART)


        else:  # If Alpha exceed equilibrium value, Alpha Dissolution mode - hierarchical, dissolving highest Gibbs energy phases first

            if EqAlpha == 0:

                deltafracMART = -fractionMART[i]
                deltafracMASS = -fractionMASS[i]
                deltafracGB = -fractionGB[i]
                deltafracBW = -fractionBW[i]
                deltafracCOL = -fractionCOL[i]

            elif fractionalpha[i] - EqAlpha < fractionMART[i] >= 0:

                deltafracMART = EqAlpha - fractionalpha[i]
                deltafracMASS = 0
                deltafracGB = 0
                deltafracBW = 0
                deltafracCOL = 0

            elif fractionMART[i] + fractionMASS[i] >= fractionalpha[i] - EqAlpha:

                deltafracMART = -fractionMART[i]
                deltafracMASS = EqAlpha - fractionalpha[i] + fractionMART[i]
                deltafracGB = 0
                deltafracBW = 0
                deltafracCOL = 0

            elif fractionMART[i] + fractionMASS[i] + fractionBW[i] >= fractionalpha[i] - EqAlpha:

                deltafracMART = -fractionMART[i]
                deltafracMASS = -fractionMASS[i]
                deltafracGB = 0
                deltafracBW = EqAlpha - fractionalpha[i] + fractionMART[i] + fractionMASS[i]
                deltafracCOL = 0

            elif fractionMART[i] + fractionMASS[i] + fractionBW[i] + fractionCOL[i] >= fractionalpha[i] - EqAlpha:

                deltafracMART = -fractionMART[i]
                deltafracMASS = -fractionMASS[i]
                deltafracGB = 0
                deltafracBW = -fractionBW[i]
                deltafracCOL = EqAlpha - fractionalpha[i] + fractionMART[i] + fractionMASS[i] + fractionBW[i]

            else:

                deltafracMART = -fractionMART[i]
                deltafracMASS = -fractionMASS[i]
                deltafracGB = EqAlpha - fractionalpha[i] + fractionMART[i] + fractionMASS[i] + fractionBW[i] + \
                              fractionCOL[i]
                deltafracBW = -fractionBW[i]
                deltafracCOL = -fractionCOL[i]

            fractionGB.append(fractionGB[i] + deltafracGB)
            fractionBW.append(fractionBW[i] + deltafracBW)
            fractionCOL.append(fractionCOL[i] + deltafracCOL)
            fractionMART.append(fractionMART[i] + deltafracMART)
            fractionMASS.append(fractionMASS[i] + deltafracMASS)

        fractionalpha.append(
            (fractionGB[i + 1] + fractionBW[i + 1] + fractionCOL[i + 1] + fractionMART[i + 1] + fractionMASS[i + 1]))

        FracWid = fractionBW[i]
        deltafracwid = deltafracBW
        fractionGBend = fractionGB[i + 1]
        fractionBWend = fractionBW[i + 1]
        fractionCOLend = fractionCOL[i + 1]
        fractionMARTend = fractionMART[i + 1]
        fractionMASSend = fractionMASS[i + 1]

        if FracWid == 0:

            tlath.append(0.3)

        elif deltafracwid > 0 and fractionalpha[i+1] > 0:

            if tlath[i] == 0:
                tlath[i] == 0.3

            LathT = NG.LathThickness(tlath[i], FracWid, deltafracwid, TempK)
            tlath.append(LathT)

        else:

            tlath.append(tlath[i])

        fraca = fractionalpha[i + 1]
        tlathend = tlath[i+1]

    return [fraca, fractionGBend, fractionBWend,fractionCOLend, fractionMASSend, fractionMARTend, tlathend]