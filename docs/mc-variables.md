
## 1. List-Style Description

1. **Electron–Ion Initial State Invariants**  
   - **`dis_twopdotk`**: Defined as

   $\text{TwoPdotk} = 2 \times P_{Incident} \cdot k_{Incident}$

   where \(P_{\text{Incident}}\) is the 4-momentum of the ion (proton or nucleus), and \(k_{\text{Incident}}\) is the 4-momentum of the incoming electron.

   - **`dis_s_e`**: Defined as

   $$
   \texttt{s\_e} = M_{\text{Ion}}^2 + m_{\text{Electron}}^2 + \texttt{TwoPdotk}.
   $$

   This is essentially the Mandelstam \(s\) for the **electron–ion** system, i.e., the total energy squared in their CM frame.

2. **Virtual Photon Variables**  
   - **`dis_twopdotq`**:

   $$
   \texttt{TwoPdotq} = 2 \times P_{\text{Incident}} \cdot q_{\text{Virtual}},
   $$

   where \(q_{\text{Virtual}}\) is the 4-momentum of the virtual photon.

   - **`dis_s_q`**:

   $$
   \texttt{s\_q} = M_{\text{Ion}}^2 + \texttt{TwoPdotq}.
   $$

   Similar to \(\texttt{s\_e}\), but now the electron is replaced by the virtual photon in computing an effective CM energy squared.

3. **Core DIS Kinematics**  
   - **`dis_q2`** \((Q^2)\): The negative four-momentum transfer squared of the virtual photon. Computed or interpolated between `Q2Min` and `Q2Max` in the generator.  
   - **`dis_xbj`** \((x_{\mathrm{Bj}})\): The Bjorken \(x\) scaling variable, also interpolated between `xMin` and `xMax`.  
   - **`dis_x_d`** \((x_D)\): A rescaled version of \(x_{\mathrm{Bj}}\) that accounts for the ion’s mass, i.e. \(\texttt{x\_Bj} \times (M_{\text{Proton}} / M_{\text{Ion}})\).  
   - **`dis_y_d`** \((y_D)\): An inelasticity parameter for the ion context, given by

     $$
     y_D = \frac{Q^2}{x_D \times \texttt{TwoPdotk}}.
     $$

   - **`dis_yplus`** \((\mathrm{Yplus})\): A factor commonly used in DIS cross sections, \(\displaystyle 1 + (1 - y_D)^2\).

4. **Proton (Ion) Momentum in Rest Frame**  
   - **`dis_pdrest`**: The magnitude of the 3-momentum

   $$
   \sqrt{\,P_{\text{Incident\_Rest}}(0)^2 + P_{\text{Incident\_Rest}}(1)^2 + P_{\text{Incident\_Rest}}(2)^2}\,
   $$

   of the incoming proton (or ion) in its own rest frame.

5. **Missing Mass**  
   - **`dis_mx2`** \((M_X^2)\): The squared invariant mass of the **hadronic final state** not explicitly accounted for by the measured or tagged particles. Computed as \(\texttt{PX\_Vertex.M2()}\) in the code, where `PX_Vertex` is the 4‐momentum of the final hadronic system.

6. **Spectator Kinematics**  
   - **`dis_alphas`** \((\alpha_S)\): The **light-cone momentum fraction** carried by the spectator. Defined as

     $$
     \texttt{alphaS} \;=\; A_{\text{Beam}} \times 
         \frac{\,pS_{\text{rest}} \cdot \cos(\theta_{\text{Recoil}})\;+\;E_{\text{Spectator\_Rest}}\,}{\,M_{\text{Ion}}\,},
     $$

     (exact expression depends on the code’s naming of `pS_rest`, `csThRecoil`, etc.).

   - **`dis_pPerpS`** \((p_{\perp S})\): The spectator’s transverse momentum in the rest frame, typically

     $$
     \texttt{pPerpS} \;=\; pS_{\text{rest}} \times \sqrt{\,1 - \cos^2(\theta_{\text{Recoil}})\,}.
     $$

7. **Other Variables**  
   - **`dis_nu`**: Often the energy transfer \(\nu = E_{\text{beam}} - E_{\text{scattered}}\).  
   - **`dis_tspectator`**, **`dis_tprime`**: The generator may compute the momentum transfer \(t\) to the spectator (\(\texttt{tspectator}\)) and \(t^\prime = t - t_{\min}\).  
   - **`dis_tempvar`**: A placeholder or debug variable, not used for final physics.

---

## 2. Table of Definitions

Below is a **condensed table** listing each variable, the **equation or code snippet**, and a **short description**. Variables are grouped for clarity:

| **Variable (dis_*)** | **Definition (Code or Formula)**                                                                                                  | **Physical Meaning**                                                                                                                 | **Comments**                                                                                                                                                                                             |
|----------------------|-----------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Electron–Ion Invariants** |                                                                                                                                   |                                                                                                                                       |                                                                                                                                                                                                           |
| `twopdotk`          | $$2 \times P_{\text{Incident}} \cdot k_{\text{Incident}}$$                                                                         | Scalar product of ion & electron 4‐momenta (with factor 2).                                                                           | Used to build total CM energy of \((e + \text{Ion})\).                                                                                                                                                  |
| `s_e`               | $$M_{\text{Ion}}^2 + m_{\text{Electron}}^2 + \texttt{twopdotk}$$                                                                   | Mandelstam \(s\) for the electron–ion system (CM energy squared).                                                                     | Shows total energy scale for the initial \((e + \text{Ion})\).                                                                                                                                          |
| **Virtual Photon**  |                                                                                                                                   |                                                                                                                                       |                                                                                                                                                                                                           |
| `twopdotq`          | $$2 \times P_{\text{Incident}} \cdot q_{\text{Virtual}}$$                                                                          | Dot product of ion & virtual photon 4‐momenta (with factor 2).                                                                        | Key for DIS variables involving the virtual photon.                                                                                                                                                    |
| `s_q`               | $$M_{\text{Ion}}^2 + \texttt{twopdotq}$$                                                                                           | “Photon–Ion” CM energy squared (replacing electron with the photon).                                                                  | Useful in analyzing the hadronic system post virtual-photon emission.                                                                                                                                   |
| **Core DIS**        |                                                                                                                                   |                                                                                                                                       |                                                                                                                                                                                                           |
| `q2`                | Interpolated: \(\;Q^2 = Q2_{\text{Max}}\,u + Q2_{\text{Min}}\,(1 - u)\)                                                            | Negative four‐momentum transfer squared of the virtual photon.                                                                        | Central to DIS. Larger \(Q^2\) = more “hard” scattering.                                                                                                                                               |
| `xbj`               | $$x_{\mathrm{Bj}} = x_{\mathrm{Min}}^{\,1-u_v}\times x_{\mathrm{Max}}^{\,u_v}$$                                                   | Bjorken \(x\), fraction of target momentum carried by struck quark.                                                                   | Dimensionless. Sometimes simply read from final code’s calculation.                                                                                                                                    |
| `x_d`               | $$x_d = x_{\mathrm{Bj}} \times \frac{M_{\text{Proton}}}{M_{\text{Ion}}}$$                                                         | Rescaled \(x\) for the ion’s mass.                                                                                                    | Typically relevant if the target is heavier than a proton (nuclear corrections).                                                                                                                      |
| `y_d`               | $$y_d = \frac{Q^2}{x_d \times \texttt{twopdotk}}$$                                                                                 | Inelasticity parameter in the ion frame.                                                                                              | Reflects fraction of electron’s energy lost to the hadronic system.                                                                                                                                     |
| `yplus`             | $$1 + (1 - y_d)^2$$                                                                                                                | Factor often used in DIS cross sections: \(1 + (1-y)^2\).                                                                              | A purely dimensionless factor weighting certain structure function terms.                                                                                                                               |
| **Proton (Ion) Momentum** |                                                                                                                                   |                                                                                                                                       |                                                                                                                                                                                                           |
| `pdrest`            | $$\sqrt{\,P_{\text{Incident\_Rest}}(0)^2 + P_{\text{Incident\_Rest}}(1)^2 + P_{\text{Incident\_Rest}}(2)^2}\,}$$                 | Magnitude of proton/ion 3‐momentum in its rest frame.                                                                                 | Often near zero if you truly are in the ion’s own rest frame (for a proton, might be a small placeholder).                                                                                              |
| **Missing Mass**    |                                                                                                                                   |                                                                                                                                       |                                                                                                                                                                                                           |
| `mx2`               | $$(\texttt{PX\_Vertex}).M^2$$                                                                                                      | Squared invariant mass of leftover hadronic system.                                                                                   | Useful for checking exclusivity or missing particles. Negative values can appear if off‐shell or partial definitions.                                                                                  |
| **Spectator Kinematics** |                                                                                                                                   |                                                                                                                                       |                                                                                                                                                                                                           |
| `alphas`            | $$A_{\text{Beam}} \times \frac{\,pS_{\text{rest}} \cdot \cos(\theta_{\text{Recoil}})\;+\;E_{\text{Spectator\_Rest}}\,}{\,M_{\text{Ion}}\,}$$ | Light‐cone fraction of the spectator/recoil baryon.                                                                                    | Common in Sullivan‐process or deuteron‐spectator analyses, indicates how much beam momentum the spectator carries.                                                                                     |
| `pPerpS`            | $$pS_{\text{rest}} \times \sqrt{\,1 - \cos^2(\theta_{\text{Recoil}})\,}$$                                                         | Transverse momentum of the spectator in the rest frame.                                                                                | Tells how far from the beam axis the spectator is kicked.                                                                                                                                              |
| **Miscellaneous**   |                                                                                                                                   |                                                                                                                                       |                                                                                                                                                                                                           |
| `nu`                | \(\;E_e - E_{e'}\) or \(\;\frac{p \cdot q}{m_p}\)                                                                                  | Energy transfer to hadronic system.                                                                                                   | Another standard DIS variable, depending on the code’s approach.                                                                                                                                        |
| `tspectator` / `tprime` | Typically \(t = (p_{\text{inc}} - p_{\text{rec}})^2\), \(t^\prime = t - t_{\min}\)                                             | Momentum transfer to the recoil system, or “excess” above minimal \(t\).                                                              | Negative in HEP sign convention. Often used in exclusive/semi‐exclusive processes.                                                                                                                     |
| `tempvar`           | —                                                                                                                                 | Temporary or debug variable.                                                                                                           | Typically zero or not used in final analysis.                                                                                                                                                            |


- **Name Consistency**: In the MC code, variables were `invts.*`; you renamed them with `dis_...`. The definitions match 1–1, except for additional variables not shown in the snippet (like `tspectator`, `tempvar`).
- **Units**: All energies and momenta are typically in GeV (or GeV^2 if squared). Dimensionless variables (\(x_{\mathrm{Bj}}, y_d, \alpha_s, yplus\)) have no units.

