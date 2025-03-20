Below is a final, combined reference for the variables in the `invts` branch (which you renamed `dis_*`), merging the code snippet from the Monte Carlo generator with earlier explanations. First is a list‐style summary, then a table that shows the definitions in plain text (no LaTeX).

---

## 1) List‐Style Summary

**A. Electron–Ion Initial State Invariants**

1. **dis_twopdotk**
    - Code definition: `TwoPdotk = 2.*(PIncident_Vertex.Dot(kIncident_Vertex))`
    - Meaning: 2 × (dot product of the ion’s 4-momentum and the electron’s 4-momentum).
    - Used to build the total center-of-mass energy of the electron + ion system.

2. **dis_s_e**
    - Code definition: `s_e = MIon*MIon + mElectron*mElectron + TwoPdotk`
    - Meaning: The Mandelstam “s” for the electron–ion system (their total energy squared in the center-of-mass frame).

**B. Virtual Photon Variables**

1. **dis_twopdotq**
    - Code definition: `TwoPdotq = 2.*(PIncident_Vertex.Dot(qVirtual_Vertex))`
    - Meaning: 2 × (dot product of the ion’s 4-momentum and the virtual photon’s 4-momentum).

2. **dis_s_q**
    - Code definition: `s_q = MIon*MIon + TwoPdotq`
    - Meaning: A version of center-of-mass energy squared where the electron is replaced by the virtual photon.

**C. Core DIS Kinematics**

1. **dis_q2 (Q2)**
    - Code definition: `Q2 = Q2Max*uu + Q2Min*(1.-uu)`
    - Meaning: The negative four-momentum transfer squared of the virtual photon.
    - Represents how “hard” the electron scatters off the ion.

2. **dis_xbj (Bjorken x)**
    - Code definition: `xBj = pow(xMin,1.-uv)*pow(xMax,uv)`
    - Meaning: Bjorken x, the fraction of the proton or nucleon momentum carried by the struck quark.

3. **dis_x_d**
    - Code definition: `x_d = xBj*(MProton/MIon)`
    - Meaning: A rescaled xBj that accounts for the ion mass (if the ion is heavier than a proton).

4. **dis_y_d (inelasticity)**
    - Code definition: `y_d = Q2/(x_d*TwoPdotk)`
    - Meaning: Fraction of the electron’s energy transferred to the target, adapted for the ion’s mass.

5. **dis_yplus (Yplus)**
    - Code definition: `Yplus = 1 + ((1-invts.y_D)*(1-invts.y_D))`
    - Meaning: A factor often appearing in cross‐section formulas, 1 + (1 - y)^2.

**D. Proton (Ion) Momentum in Rest Frame**

1. **dis_pdrest**
    - Code definition:
      ```
      pDrest = sqrt(PIncident_Rest(0)*PIncident_Rest(0)
                     + PIncident_Rest(1)*PIncident_Rest(1)
                     + PIncident_Rest(2)*PIncident_Rest(2));
      ```
    - Meaning: Magnitude of the 3-momentum of the incoming proton/ion in its own rest frame.

**E. Missing Mass**

1. **dis_mx2**
    - Code definition: `MX2 = PX_Vertex.M2();`
    - Meaning: Squared invariant mass of the remaining hadronic system.
    - Can help identify exclusive vs. inclusive events or hidden final-state particles.

**F. Spectator Kinematics**

1. **dis_alphas (alphaS)**
    - Code definition:
      ```
      alphaS = ABeam*(pS_rest*csThRecoil + pSpectator_Rest.E()) / MIon
      ```
    - Meaning: The light-cone momentum fraction for the spectator nucleon or recoil baryon.
    - Tells how much of the beam momentum the spectator carries along the beam direction.

2. **dis_pPerpS**
    - Code definition:
      ```
      pPerpS = pS_rest*sqrt(1.-csThRecoil*csThRecoil);
      ```
    - Meaning: The transverse momentum of the spectator in its rest frame, relative to the beam axis.

**G. Other Variables in Code Snippet**

- **dis_nu** (energy transfer), **dis_tspectator** (t), **dis_tprime** (t minus t_min), and **dis_tempvar** (dummy) were mentioned in prior discussions. They may appear in the final file but are not explicitly defined in the snippet above.

---

## 2) Table of Definitions (Plain Text)

Below is a concise table grouping these variables by category, showing the code snippet or formula, their meaning, and any extra notes. (Equations are given as text for clarity.)

```
┌───────────────────────┬──────────────────────────────────────────────────────────────────┬───────────────────────────────────────────────────────────────────────────────────────────┐
│ Electron–Ion Invariants                                                                                                             │
├───────────────────────┼──────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
│ dis_twopdotk          │ TwoPdotk = 2. * PIncident_Vertex.Dot(kIncident_Vertex)         │ Dot product (times 2) of the ion 4-momentum and electron 4-momentum.                     │
│ dis_s_e               │ s_e = MIon*MIon + mElectron*mElectron + TwoPdotk               │ Mandelstam "s" for e + ion system (total CM energy squared).                              │
├───────────────────────┼──────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
│ Virtual Photon                                                                                                                        │
├───────────────────────┼──────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
│ dis_twopdotq          │ TwoPdotq = 2. * PIncident_Vertex.Dot(qVirtual_Vertex)          │ Dot product (times 2) of ion 4-momentum and virtual photon 4-momentum.                    │
│ dis_s_q               │ s_q = MIon*MIon + TwoPdotq                                     │ CM energy squared but replacing the electron with the virtual photon.                     │
├───────────────────────┼──────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
│ Core DIS Kinematics                                                                                                                   │
├───────────────────────┼──────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
│ dis_q2                │ Q2 = Q2Max*uu + Q2Min*(1.-uu)                                   │ Negative four-momentum transfer of the virtual photon.                                    │
│ dis_xbj               │ xBj = pow(xMin,1.-uv)*pow(xMax,uv)                             │ Bjorken x, fraction of momentum carried by the struck quark.                               │
│ dis_x_d               │ x_d = xBj*(MProton / MIon)                                      │ Rescaled xBj for an ion heavier than a proton.                                            │
│ dis_y_d               │ y_d = Q2 / (x_d * TwoPdotk)                                     │ Inelasticity in that ion context.                                                         │
│ dis_yplus             │ yplus = 1 + (1 - y_d)*(1 - y_d)                                 │ Factor in DIS cross sections: 1 + (1 - y)^2.                                              │
├───────────────────────┼──────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
│ Proton (Ion) Momentum in Rest Frame                                                                                                   │
├───────────────────────┼──────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
│ dis_pdrest            │ pDrest = sqrt( sum of squares of PIncident_Rest(...) )         │ Magnitude of the ion's 3-momentum in its own rest frame.                                  │
├───────────────────────┼──────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
│ Missing Mass                                                                                                                           │
├───────────────────────┼──────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
│ dis_mx2               │ MX2 = PX_Vertex.M2()                                           │ Squared invariant mass of the unobserved hadronic system.                                 │
├───────────────────────┼──────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
│ Spectator Kinematics                                                                                                                  │
├───────────────────────┼──────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
│ dis_alphas            │ alphaS = ABeam*(pS_rest*csThRecoil + pSpectator_Rest.E())      │ Light-cone momentum fraction of the spectator (or recoil baryon).                          │
│                       │              / MIon                                             │                                                                                           │
│ dis_pPerpS            │ pPerpS = pS_rest* sqrt(1 - csThRecoil*csThRecoil)              │ Transverse momentum of the spectator in its rest frame.                                    │
└───────────────────────┴──────────────────────────────────────────────────────────────────┴───────────────────────────────────────────────────────────────────────────────────────────┘
```

**Extra Notes**  
• Units: Typically, all energies and momenta are in GeV or GeV^2 (for squared quantities).  
• Variables like xBj, y_d, alphaS are dimensionless fractions.  
• Some additional variables (dis_nu, dis_tspectator, dis_tprime, dis_tempvar) might also appear in your files but were not fully defined in the specific code snippet.

This wraps up the consolidated reference, showing how each `dis_*` variable is defined in the MC generator code and describing what it means for the final analysis.