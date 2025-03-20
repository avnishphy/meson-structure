**Below is a single, combined markdown table** for the `invts` (now `dis_*`) variables.
- Three columns:
    1. **Variable (dis_*)**
    2. **Definition** (code snippet or formula in plain text)
    3. **Physical Meaning**
- The **Comments** column is omitted per your request.

| **Variable (dis_*)**        | **Definition**                                                          | **Physical Meaning**                                                             |
|-----------------------------|-------------------------------------------------------------------------|----------------------------------------------------------------------------------|
| **Electron–Ion Invariants** |                                                                         |                                                                                  |
| twopdotk                    | 2 * PIncident_Vertex.Dot(kIncident_Vertex)                              | Scalar product (times 2) of ion 4-momentum and incident electron 4-momentum      |
| s_e                         | MIon^2 + mElectron^2 + twopdotk                                         | Total squared energy (Mandelstam s) of the electron–ion system                   |
| **Virtual Photon**          |                                                                         |                                                                                  |
| twopdotq                    | 2 * PIncident_Vertex.Dot(qVirtual_Vertex)                               | Scalar product (times 2) of ion 4-momentum and virtual photon 4-momentum         |
| s_q                         | MIon^2 + twopdotq                                                       | Modified CM energy squared, replacing electron with virtual photon               |
| **Core DIS Kinematics**     |                                                                         |                                                                                  |
| q2                          | Q2 = Q2Max * uu + Q2Min * (1 - uu)                                      | Negative four-momentum transfer squared (Q^2) of the virtual photon              |
| xBj                         | xBj = (xMin)^(1 - uv) * (xMax)^(uv)                                     | Bjorken-x, fraction of the target nucleon’s momentum carried by the struck quark |
| x_d                         | x_d = xBj * (MProton / MIon)                                            | Rescaled x to account for the ion’s mass                                         |
| y_d                         | y_d = q2 / (x_d * twopdotk)                                             | Inelasticity parameter in ion-level kinematics                                   |
| Yplus                       | 1 + (1 - y_d) * (1 - y_d)                                               | Common DIS factor 1 + (1-y)^2 in the generator’s notation                        |
| **Proton (Ion) Momentum**   |                                                                         |                                                                                  |
| pDrest                      | sqrt( PIncident_Rest(0)^2 + PIncident_Rest(1)^2 + PIncident_Rest(2)^2 ) | Magnitude of the proton/ion 3-momentum in its own rest frame                     |
| **Missing Mass**            |                                                                         |                                                                                  |
| MX2                         | PX_Vertex.M2()                                                          | Squared invariant mass of the unobserved (remaining) hadronic system             |
| **Spectator Kinematics**    |                                                                         |                                                                                  |
| alphaS                      | ABeam * ( pS_rest * csThRecoil + pSpectator_Rest.E() ) / MIon           | Light-cone momentum fraction of the spectator/recoil nucleon or baryon           |
| pPerpS                      | pS_rest * sqrt( 1 - csThRecoil * csThRecoil )                           | Transverse momentum of the spectator in the rest frame                           |
| **Miscellaneous**           |                                                                         |                                                                                  |
| nu                          | E_incident_e - E_scattered_e (or similar)                               | Energy transfer from the electron to the hadronic system                         |
| t_spectator, tprime         | (p_in - p_recoil)^2, and tprime = t - t_min (not explicitly in snippet) | Momentum transfer to the spectator or “excess” above minimal t                   |
| tempvar                     | (placeholder variable)                                                  | Debug or intermediate quantity, often zero                                       |

- **Units:**
    - Most energy and momentum quantities are in GeV, and squared forms (like s_e, s_q, q2, MX2) are in GeV^2.
    - Dimensionless variables include xBj, x_d, alphaS, y_d, Yplus.

- **Notes:**
    - Variable naming in code was `invts.TwoPdotk`, etc. You have them as `dis_twopdotk`, etc.
    - The definitions come directly from the MC generator snippet you provided.