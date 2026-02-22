# ROCK SLOPE STABILITY CALCULATOR

## 1. Technical Overview and Scientific Context

### Purpose and Functionality of the Computational Tool
The **Rock Slope Stability Calculator** is a specialized computational tool designed to facilitate the hydraulic design and assessment of rock armoured slopes commonly used in coastal protection, breakwaters, groynes, revetments, etc.

Unlike simplified empirical tools that rely on single formulas and static coefficients, this calculator integrates the dynamic evolution of coastal engineering research spanning the last four decades. It performs a comprehensive analysis of hydraulic stability across different water depth zones (**Deep, Shallow, Very Shallow, and Swash zones**), selecting the most appropriate empirical formula based on the hydraulic regime.

The primary function of this tool is to calculate the required size ($D_{n50}$) and weight ($M_{50}$) of rock armor units. It synthesizes the foundational work of **Van der Meer (1988)** with contemporary re-evaluations concerning spectral wave periods, shallow water non-linearities, and advanced rock placement techniques. By automating the complex decision-making process required to select the correct stability formula, it reduces the risk of design errors in non-standard hydraulic conditions.

### Methodology: Multi-Model Consensus & Manual Override
The calculator implements a "Multi-Model Consensus" approach. Rather than relying on a single algorithm, it computes stability using eight state-of-the-art empirical formulas derived from extensive hydraulic model testing. It then evaluates the hydraulic context (specifically the Relative Depth ratio) to recommend the scientifically most accurate result.

**Manual Selection Mode:**
The software includes an interactive logic for expert users. After the initial analysis, the script asks: *"Do you want to manually choose a stability formula? (y/n)"*. If the user inputs **True/y**, the system immediately prompts the user to select one of the specific formulas from the list below:

1.  Hudson (1959)
2.  Van der Meer (2021)
3.  Van Gent Modified (2003)
4.  Van Gent Simplified (2003)
5.  Eldrup & Andersen (2019)
6.  Etemad-Shahidi (2020)
7.  Scaravaglione Mod. VG (2025)
8.  Scaravaglione Mod. ES (2025)

The logic follows a rigorously defined **4-Step Process**:
1.  **Input Processing:** Parsing wave data ($H_{m0}$, $T_{m-1,0}$), geometry (slopes, depths), and material properties (density, permeability).
2.  **Hydraulic Analysis:** Computing local wavelength, celerity, wave steepness, surf similarity (Iribarren number), and determining the breaker type (Surging vs. Plunging).
3.  **Stability Calculation:** Running multiple empirical models simultaneously.
4.  **Intelligent Selection:** Automatically detecting the "Hydraulic Zone 1-4" (Deep vs. Shallow vs. Swash) to recommend the most valid result based on the latest literature.

---

## 2. Epistemological Framework

The stability of granular coastal structures is a problem of immense physical complexity, governing the interaction between a stochastic hydrodynamic load and a porous, flexible structural matrix.

**From Static Coefficients to Dynamic Interaction:**
Early engineering approaches, such as the **Hudson formula** (1959), treated this interaction through a monolithic stability coefficient ($K_D$), which aggregated all unknown variables—permeability, storm duration, wave period, and breaking characteristics—into a single scalar. While computationally convenient, this approach failed to capture the physical nuance of failure mechanisms. It implicitly assumed that the wave period had no effect on stability, an assumption disproven by modern research which shows that longer wave periods (associated with higher spectral density) penetrate deeper into the core, altering pore pressures and lift forces.

**The Role of the Breaker Index:**
The transition to the **Van der Meer formulae** (1988) marked a paradigm shift, introducing the surf similarity parameter (or breaker index, $\xi$) as a governing variable to distinguish between failure modes.
* **Plunging Waves:** Characterized by high velocity impacts and void extraction.
* **Surging Waves:** Characterized by high run-up/run-down and phreatic surface lag, leading to uplift.
However, original derivations relied on the mean wave period ($T_m$) and were constrained largely to deep-water conditions where the Rayleigh distribution of wave heights holds valid.

**Spectral Characterization and Shallow Water Complexity:**
Contemporary coastal engineering has moved toward spectral characterization, utilizing the spectral energy period ($T_{m-1,0}$) to better account for bimodal seas and shallow water spectral flattening. In shallow water, the wave height distribution deviates from Rayleigh, and the shape of the spectrum changes due to triad wave interactions. This necessitates formulas that include parameters like the wave height ratio ($H_{2\%}/H_{s}$) or physical permeability coefficients ($C_p$) rather than notional ones.

**The "Permissible Damage" Philosophy:**
Modern design does not seek "zero movement" but rather "permissible damage" over a service life. This probabilistic approach acknowledges that granular structures are flexible. They can reshape under extreme loads without catastrophic failure, unlike rigid concrete structures. The calculator's use of the damage parameter ($S_d$) allows engineers to design for specific performance levels, from "Start of Damage" to "Intermediate Damage" (reshaping).

**The Frontier of Knowledge:**
Despite these advancements, the physics of extremely shallow water remains a frontier of research where standard assumptions often break down. As explicitly stated by **Van der Meer et al. (2024)**: 
> *"In summary, for shallow water conditions with $h/H_{m0,\mathrm{deep}} < 1.5$, there is currently no reliable method to describe the stability of rock slopes under wave attack."*

This software utilizes the newest 2025 research, specifically the works by Van der Meer et al., Van Gent et al. and Scaravaglione et al., to bridge this specific gap and provide valid predictions where traditional methods fail.

---

## 3. Hydrodynamic Input Parameters and Spectral Characterization

The accurate characterization of the hydraulic load is the single most critical step in stability analysis. The calculator requires a distinction between offshore conditions and local conditions at the toe of the structure. The transformation of waves as they propagate from deep water onto a foreshore alters not only their height but their period, shape, and statistical distribution.

### The Shift to Spectral Parameters
Historically, stability formulae utilized time-domain parameters derived from zero-crossing analysis, such as the significant wave height $H_{1/3}$ and the mean period $T_m$. While intuitive, these parameters become unstable in shallow water where wave breaking disrupts the coherence of the zero-crossing definition, often resulting in "split" waves or merged crests that skew statistical counts.

To address this, modern guidance, including the **EurOtop** manual and the revisited Van der Meer formulae, mandates the use of spectral parameters. This calculator is built fundamentally upon the spectral significant wave height ($H_{m0}$) and the spectral energy period ($T_{m-1,0}$).

#### Spectral Significant Wave Height ($H_{m0}$)
The parameter $H_{m0}$ is derived from the zeroth moment of the wave energy density spectrum, $m_0$:

$$H_{m0} = 4 \sqrt{m_0}$$

In deep water, for a standard sea state, $H_{m0} \approx H_{1/3}$. However, in shallow water, nonlinear shoaling creates peaked crests and flat troughs (Stokes V or Cnoidal waves). This nonlinearity causes $H_{1/3}$ (measured in the time domain) to increase relative to $H_{m0}$.

Recent research by **Eldrup and Andersen (2019)** and **Van der Meer et al. (2024)** indicates that using $H_{1/3}$ in stability formulae for nonlinear shallow water waves can lead to overly conservative or erratic results. $H_{m0}$ remains a more robust measure of the total energy contained in the wave train, regardless of the individual wave shape distortions. Consequently, this calculator prioritizes $H_{m0}$ as the primary input for stability calculations, particularly when the relative depth $h/H_{m0,deep} < 3.0$.

#### Spectral Energy Period ($T_{m-1,0}$)
The spectral period $T_{m-1,0}$ is defined by the ratio of the first negative moment ($m_{-1}$) and the zeroth moment ($m_0$) of the spectrum:

$$T_{m-1,0} = \frac{m_{-1}}{m_0} = \frac{\int_{0}^{\infty} f^{-1} S(f) df}{\int_{0}^{\infty} S(f) df}$$

This definition places greater weight on the lower frequency components of the spectrum. This is physically significant for rock stability because larger, longer-period waves contain more energy and generate stronger surficial flows than short-period chop. 

**Note on Jonswap Spectra:** For a standard Jonswap spectrum, the relationship between the spectral energy period $T_{m-1,0}$ and the peak period $T_p$ can be approximated by:
$$T_p \approx 1.1 \cdot T_{m-1,0}$$

### Wave Height Distribution and $H_{2\%}$
The parameter $H_{2\%}$ represents the wave height exceeded by 2% of the waves in a storm record. In deep water, wave heights generally follow a Rayleigh distribution, where the relationship is constant:

$$
\frac{H_{2\%}}{H_{s}} \approx 1.4
$$

In shallow water, this ratio changes dynamically. As waves enter the shoaling zone, the distribution broadens, and the ratio may drop to nearly 1.2 near the breaker line before recovering. The "Van Gent Modified" formula (2003) explicitly utilizes $H_{2\%}$. However, predicting $H_{2\%}$ accurately without physical modelling is fraught with uncertainty. Because of this epistemic uncertainty, this calculator generally favors the rewritten formulae based on $H_{m0}$ and $T_{m-1,0}$ for design, using $H_{2\%}$ based methods only for comparative verification.

### The Role of Foreshore and Water Depth
The geometry of the foreshore is a critical input. The calculator requires the local water depth at the toe of the structure ($h_{toe}$) and the foreshore slope ($m$). These define the **Relative Depth Ratio**:

$$R_h = \frac{h_{toe}}{H_{m0,deep}}$$

This dimensionless parameter acts as the primary switch in the calculator's algorithm, determining whether the structure is in a Deep, Shallow, or Extremely Shallow regime.

---

## 4. Intelligent Design Engine: Hydraulic Zones & Selection Logic

The calculator divides the coastal profile into four distinct zones based on the latest research by *Van der Meer et al. (2024)* and *Scaravaglione et al. (2025)*:

* **Zone 1: Deep / Intermediate Water ($h/H_{m0} > 3.0$)**
    * **Physics:** Rayleigh wave distribution ($H_{2\%}/H_{s} \approx 1.4$). No depth-induced breaking.
    * **Recommended Formula:** **[2] Van der Meer (2021)**.
    * **Rationale:** The rewritten formula using spectral period $T_{m-1,0}$ is the industry standard for non-depth-limited conditions.

* **Zone 2: Shallow Water ($1.5 < h/H_{m0} \le 3.0$)**
    * **Physics:** Wave shoaling and spectral truncation.
    * **Recommended Formula:** **[2] Van der Meer (2021)**.
    * **Rationale:** *Van der Meer et al. (2024)* confirmed the formula's validity extends down to $h/H_{m0} = 1.5$, preferring $H_{m0}$ over $H_{1/3}$ to handle nonlinearity.

* **Zone 3: Very Shallow Water ($0.5 < h/H_{m0} \le 1.5$)**
    * **Physics:** Surf zone. Severe breaking, bore-like flow, and horizontal stability trends.
    * **Recommended Formula:** **[8] Scaravaglione Mod. ES (2025)**.
    * **Rationale:** Standard formulas diverge here. *Scaravaglione et al. (2025)* demonstrated that modifying the Etemad-Shahidi approach prevents the dangerous over-prediction of stability common in this zone.

* **Zone 4: Extremely Shallow / Swash Zone ($h/H_{m0} \le 0.5$)**
    * **Physics:** High turbulence, air entrainment, and lack of buoyancy.
    * **Recommended Formula:** **[7] Scaravaglione Mod. VG (2025)**.
    * **Rationale:** Standard coefficients are unsafe here. This formula uses a recalibrated coefficient ($C_{VG}=3.3$) to account for swash zone instability.

---

## 5. Structural Parameters and Material Characterization

While the hydrodynamic inputs define the load, the structural inputs define the resistance. The stability of a rubble mound is not merely a function of rock mass but of the geometric configuration that dictates energy dissipation.

### Slope Angle ($\alpha$)
The slope angle is entered as $\cot\alpha$ (e.g., 1.5, 2.0, 3.0, 4.0). The slope angle is the primary determinant of the breaker type.
* **Steep Slopes ($\cot\alpha < 2.0$):** Promote surging or collapsing breakers. The wave energy is reflected or surges up the slope, creating significant lift forces during down-rush.
* **Gentle Slopes ($\cot\alpha > 3.0$):** Promote plunging breakers. The wave crest curls over and impacts the slope, creating high localized pressure shocks and scouring forces.

The calculator limits inputs to $1.5 \le \cot\alpha \le 6.0$, aligning with the valid range of the underlying experimental datasets.

### Notional Permeability ($P$)

The notional permeability factor $P$ is a dimensionless parameter that quantifies the hydraulic conductivity of the structure's core and underlayers. It ranges theoretically from 0.0 (completely impervious) to 1.0, but practical limits for rock structures are typically between 0.1 and 0.6.

* **P = 0.1 (Impermeable Core):** This value represents a structure where the armor layer rests directly on an impermeable base, such as a clay core, geotextile, or sand body. In this configuration, wave energy cannot penetrate the structure. Consequently, all energy must be either dissipated by friction on the outer slope or reflected. This results in high run-up and run-down velocities, which increases hydraulic instability.
* **P = 0.4 (Standard Breakwater):** This value is typical of many harbour breakwaters. It describes a standard cross-section consisting of an armor layer, a filter layer, and a semi-permeable core. This configuration achieves a balance between energy dissipation and penetration.
* **P = 0.5 (Permeable Core):** This value indicates a structure with an armor layer placed on a filter over a remarkably open core. Significant wave energy dissipates into the large void spaces, which helps relieve pore pressure during wave drawdown and stabilizes the armor layer.
* **P = 0.6 (Homogeneous Structure):** This value represents a breakwater composed entirely of armor-sized stone. While rarely built due to high costs, this homogeneous structure offers maximum energy dissipation and the highest hydraulic stability.

**Insight:** The influence of $P$ is non-linear and strictly regime-dependent. In the **plunging regime**, a permeable core ($P=0.5$) increases stability by dampening the impact shock of the breaking wave. In the **surging regime**, permeability becomes even more critical; as demonstrated in Equation 5, $P$ appears as an exponent to the breaker parameter. For impermeable structures ($P=0.1$), stability is nearly independent of the wave period in the surging regime, whereas for permeable structures, stability increases drastically with longer wave periods.

### Damage Definition ($S_d$)
Damage to a rock slope is defined not by catastrophic failure but by the displacement of a certain volume of material. The non-dimensional damage level $S_d$ is defined as:

$$S_d = \frac{A_e}{D_{n50}^2}$$

Where $A_e$ is the cross-sectional erosion area.
* **$S_d = 2$:** "Start of Damage." This corresponds to the displacement of roughly 0-5% of stones. It is the typical design limit for static stability.
* **$S_d = 8$ (Slope 1:2) / $S_d = 17$ (Slope 1:4):** "Failure." This corresponds to the exposure of the filter layer. Once the filter is exposed, rapid unraveling of the structure can occur.

The calculator utilizes the relationship $S_d \propto N^{0.5}$ to adjust stability calculations for storm duration. Users must input the number of waves $N$. The validity of the power law $S_d = f(H^5)$ requires that $S_d$ remains within the gradual failure zone ($S_d < 15$). Extreme damage predictions ($S_d > 20$) are flagged as invalid as the structure's geometry would have fundamentally changed.

---

## 6. Rock Slope Stability Formulas

This calculator incorporates eight specific stability formulas, representing the historical progression and current state-of-the-art in coastal engineering.

### [1] Hudson (1959)
The **Hudson formula (1959)** represents the **legacy standard** for rubble mound stability analysis, historically enshrined in the Shore Protection Manual (SPM) and the Coastal Engineering Manual (CEM). While largely superseded by modern spectral methods, it is included in this calculator to provide a **critical baseline for comparison** against historical designs. The fundamental equation calculates rock slopes stability by relating the design wave height to the weight of the armour units. Here, the design wave height is explicitly defined as **(the average of the highest 10% of waves)**, which corresponds mathematically to **$H_{1/10}$**. This formula relies heavily on the **stability coefficient, $K_D$**, a monolithic parameter derived from regular wave physical model tests at the U.S. Army Corps of Engineers Waterways Experiment Station. This single coefficient aggregates multiple complex physical variables—such as the interlocking degree of the armour units, the wave shape, and the permeability of the core—into a scalar value, simplifying the design process but obscuring the specific influence of each variable.

A critical evolution in the formula's application occurred with the publication of the **1984 Shore Protection Manual (SPM)**, which introduced significant conservatism compared to earlier 1977 guidance. In this pivotal update, the design wave height was shifted from the significant wave height to effectively **increasing the hydraulic load by 27%** to address failures observed in previous designs. Concurrently, the SPM (1984) drastically revised the recommended values for breaking wave conditions, **lowering the coefficient from 3.5 to 2.0** to account for the intense instability caused by turbulent breaking forces, while retaining a $K_D$ of **4.0 for non-breaking waves**. This sharp distinction created a discontinuity in design that modern research often challenges; notably, Van der Meer has argued that for structures with a permeable core, a $K_D$ of **4.0** is generally appropriate regardless of whether the waves are breaking or non-breaking, suggesting that the Hudson formula’s rigidity often leads to inefficient designs in permeable structures.

Despite its historical ubiquity, the Hudson formula possesses profound limitations when viewed through the lens of modern hydrodynamics. Its most significant deficiency is the **complete omission of the wave period and storm duration** from the stability calculation. By ignoring the wave period, the formula cannot physically distinguish between **"plunging" breakers**, which destroy slopes via high-velocity extraction during down-rush, and **"surging" breakers**, which lift armour units via pore-pressure buildup during up-rush. This lack of period dependence can lead to **dangerously unsafe designs in conditions of long-period swell** (where surging dominates) or overly conservative designs in short-period seas. Furthermore, because it lacks a damage parameter, it assumes a binary "fail/no-fail" condition rather than predicting the progressive damage evolution of a flexible rock slope. However, its continued use is justified by its vast library of $K_D$ values for artificial concrete units and its utility in preliminary sizing or forensic analysis of older coastal infrastructure.

To adapt this static legacy method to the calculator's advanced zonated logic, the software does not use a fixed coefficient but instead applies a **Dynamic Selection based on the active hydraulic zone**. The algorithm assesses the relative depth to adjust the coefficient, mimicking the physical destabilization observed in shallower waters. In **Deep Water**, the calculator applies the standard non-breaking $K_D$ of **4.0**. As the relative depth decreases into the **Shallow Water** zone, the coefficient is reduced to **3.5**. In the **Very Shallow Water** zone, it is further lowered to **3.0**, and finally, for the highly turbulent **Swash Zone**, a conservative $K_D$ of **2.0** is automatically applied. This dynamic adjustment ensures that even the legacy baseline reflects the increased hydraulic demand of the surf zone, aligning the 1959 formula with the breaking wave constraints introduced in the 1984 SPM.

**Equation:**
$$\frac{H_{1/10}}{\Delta D_{n50}} = \frac{1.27 H_{s}}{\Delta D_{n50}} = (K_D \cot \alpha)^{1/3}$$

**Explanation:**
* **Purpose:** Developed by the USACE to provide a simple, robust design tool for rubble mounds using regular waves. It aggregates complex physics into the **$K_D$ coefficient**.
* **Capabilities:** Excellent for preliminary sizing and comparison. It is mathematically simple and widely understood. It forms the basis of many regulatory guidelines.
* **Limitations:** It **ignores the wave period ($T$) and storm duration ($N$)**. It treats plunging and surging waves identically, which can lead to unsafe designs for long-period swells (surging) or overly conservative designs for short steep waves.
* **Historical Note on $K_D$:** In the Shore Protection Manual of 1984 (SPM, 1984), not only was the coefficient 1.27 introduced (shifting design wave from $H_{s}$ to $H_{1/10}$), but the value of $K_D$ for breaking waves was revised and decreased from **3.5 to 2.0**, while for non-breaking waves it remained **4.0**.
* **Recommendation:** In contrast to the SPM (1984), research by Van der Meer suggests it is recommended to use **$K_D=4$** for the design of structures with a permeable core, irrespective of whether the conditions are with breaking waves on the foreshore or not.

**Dynamic $K_D$ Selection:**
The calculator updates the stability coefficient $K_D$ based on the hydraulic zone to reflect increased instability in shallow water:
* **Deep Water ($h/H_{m0} > 3.0$):** $K_D = 4.0$
* **Shallow Water ($1.5 < h/H_{m0} \le 3.0$):** $K_D = 3.5$
* **Very Shallow Water ($0.5 < h/H_{m0} \le 1.5$):** $K_D = 3.0$
* **Swash Zone ($h/H_{m0} \le 0.5$):** $K_D = 2.0$

### [2] Van der Meer (2021)
The **Van der Meer (2021)** formulation represents the **contemporary industry standard** for the stability design of rock-armoured slopes, serving as a rigorously updated successor to the foundational 1988 equations that first introduced the critical distinction between plunging and surging failure modes. While the original 1988 work was derived using time-domain parameters such as the significant wave height and mean wave period, modern coastal engineering has universally shifted toward spectral analysis. To align with this paradigm, Van der Meer explicitly rewrote the stability formulae to utilize the **spectral significant wave height** and, most crucially, the **spectral energy period**. This specific wave period is defined by the ratio of the first negative moment and the zeroth moment of the energy spectrum, placing greater statistical weight on the **longer period waves** that contain the bulk of the energy and are responsible for the deepest hydraulic penetration into the structure's core. By adopting the formulation effectively eliminates the dependency on spectral shape, making the calculated stability valid for both standard single-peaked storms and **complex bimodal sea states** without requiring arbitrary conversion factors.

A defining capability of this formula is its comprehensive physical modeling of the wave-structure interaction, which it achieves by dynamically calculating a **critical surf similarity parameter**, to delineate the transition between two distinct hydrodynamic regimes. In the **"plunging" regime**, where waves break violently upon the slope, the formula models the dominant failure mechanism as the **extraction of armour units** during the rapid down-rush of water, a process heavily influenced by the external friction and the interlocking of the stones. In the **"surging" regime**, where waves do not break but rather rush up the face of the slope, the failure mechanism shifts to the **uplift and sliding of the armour layer** caused by high internal pore pressures and reduced effective weight. The formula accounts for this via the **Notional Permeability Factor**, which mathematically describes the capacity of the core and underlayers to dissipate hydraulic energy; a permeable core significantly relieves these destabilizing pressures compared to an impermeable core, thereby increasing the calculated stability, particularly in the surging regime.

Beyond the hydrodynamics, the Van der Meer (2021) formulation introduces a **performance-based design philosophy** that contrasts sharply with the binary "fail/no-fail" approach of the legacy Hudson method. It explicitly incorporates the **number of waves** and a **quantifiable damage level** into the stability equation, generally following a power law where damage accumulates relative to the square root of the storm duration. This allows engineers to design structures for specific performance targets, ranging from "Start of Damage" to "Intermediate Damage" (reshaping), providing a nuanced tool for lifecycle cost analysis. However, the formula retains a critical limitation: it relies on the assumption that the wave heights follow a **Rayleigh distribution**, where the ratio of the exceedance wave height to the significant wave height is approximately **1.4**. In shallow foreshores where depth-induced breaking truncates the wave distribution, this ratio decreases, potentially rendering the standard formula inaccurate if not corrected or substituted by shallow-water specific alternatives like those proposed by Van Gent or Scaravaglione.

**Surf Similarity Parameter ($\xi_{m-1,0}$):**
The breaker type is determined using the Iribarren number calculated with the spectral energy period:

$$
\xi_{m-1,0} = \frac{\tan \alpha}{\sqrt{s_{m-1,0}}} \quad \text{where} \quad s_{m-1,0} = \frac{2\pi H_{m0}}{g T_{m-1,0}^2}
$$

**Critical Surf Similarity Parameter ($\xi_{cr}$):**
The transition from plunging to surging waves is defined by the critical surf similarity parameter, calculated as:

$$
\xi_{cr} = \left[ \frac{c_{pl}}{c_{su}} \cdot P^{0.31} \cdot \sqrt{\tan \alpha} \right]^{\frac{1}{P+0.5}}
$$

Where $c_{pl} = 6.49$ and $c_{su} = 0.97$.

**Plunging Waves ($\xi_{m-1,0} < \xi_{cr}$):**
$$N_s = c_{pl} \cdot P^{0.18} \cdot \left(\frac{S}{\sqrt{N}}\right)^{0.2} \cdot \xi_{m-1,0}^{-0.5}$$

**Surging Waves ($\xi_{m-1,0} \ge \xi_{cr}$):**
$$N_s = c_{su} \cdot P^{-0.13} \cdot \left(\frac{S}{\sqrt{N}}\right)^{0.2} \cdot \sqrt{\cot \alpha} \cdot \xi_{m-1,0}^{P}$$

**Explanation:**
* **Purpose:** To provide a scientifically rigorous stability prediction that accounts for the influence of wave period, permeability, and storm duration in deep to intermediate water.
* **Capabilities:** It is the **current state-of-the-art for standard breakwater design**. By using $T_{m-1,0}$, it handles bimodal spectra effectively. It accurately predicts the transition from erosive (plunging) failure to sliding (surging) failure.
* **Limitations:** It assumes a **Rayleigh wave height distribution ($H_{2\%}/H_{s} \approx 1.4$)**. In shallow water where waves are depth-limited and the distribution is truncated, this assumption may lead to inaccuracies if not corrected.

### [3] Van Gent Modified (2003)
The **Van Gent Modified (2003)** formulation represents a targeted recalibration of the original Van der Meer equations, specifically engineered to address the hydrodynamic complexities of **shallow foreshores** where standard deep-water assumptions fail. Originating from the experimental insights of Smith et al. (2002) and rigorously extended by Van Gent et al. (2003), this method utilized a specialized dataset of 207 physical model tests conducted with 1:30 and 1:100 foreshore slopes. The primary scientific motivation was to correct the inaccuracies observed when applying Rayleigh-based stability models to depth-limited wave conditions. In deep water, the ratio of the exceedance wave height to the significant wave height is constant at approximately 1.4; however, as waves propagate onto a shallow foreshore, **depth-induced breaking truncates the wave height distribution**, significantly reducing this ratio. The Van Gent Modified formula explicitly integrates this ratio as a governing parameter with an inverse relationship (exponent -1), thereby dynamically adjusting the calculated stability based on the actual statistical distribution of the wave field at the structure's toe.

A further critical refinement in this formulation is the strict adherence to the **spectral energy period**, rather than the peak period or mean period used in earlier works. This choice is physically robust for shallow water environments where the energy spectrum often becomes double-peaked or flattened due to triad wave interactions and surf beat. By basing the surf similarity parameter on this spectral period, the formula ensures that the hydraulic response—whether plunging or surging—is characterized by the wave components that actually carry the bulk of the destructive energy. The formulation retains the dual-regime logic of Van der Meer, providing distinct coefficients for plunging waves (dominated by void extraction) and surging waves (dominated by uplift), but **recalibrates these coefficients to 8.4 and 1.3** respectively to better align with the shallow water experimental data. This makes the formula particularly valuable for design scenarios in the **"transitional" depth zone** where waves are shoaling heavily but not yet fully saturated bores, a region where other formulas may underestimate the required armor mass.

**Modified Plunging ($\xi_{m-1,0} < \xi_{cr}$):**
$$
\frac{H_{s}}{\Delta D_{n50}} = 8.4 P^{0.18} \left( \frac{S_d}{\sqrt{N}} \right)^{0.2} \left( \frac{H_{s}}{H_{2\%}} \right) \xi_{m-1,0}^{-0.5}
$$

**Modified Surging ($\xi_{m-1,0} \ge \xi_{cr}$):**
$$
\frac{H_{s}}{\Delta D_{n50}} = 1.3 P^{-0.13} \left( \frac{S_d}{\sqrt{N}} \right)^{0.2} \left( \frac{H_{s}}{H_{2\%}} \right) \sqrt{\cot \alpha} \xi_{m-1,0}^{P}
$$

*(Notes: The Rock Manual suggests using optimized coefficients closer to $c_{pl(mod)} = 8.4$ and $c_{s(mod)} = 1.3$ when the explicit wave height ratio $H_{2\%}/H_{s}$ is included. Exact formulation may vary slightly between subsequent publications or recalibrations).*

**Explanation:**
* **Purpose:** Developed to improve the accuracy of rock slope design in water depth conditions more limited than those in Van der Meer's original dataset. It incorporates the influence of foreshore breaking and water depth.
* **Capabilities:**
    * **Spectral Period:** Uses **$T_{m-1,0}$ instead of $T_m$ or $T_p$**, making it suitable for bimodal spectra or flat energy spectra occurring in shallow water with intense breaking.
    * **Wave Transformation:** explicit inclusion of the ratio **$H_{2\%}/H_{s}$** as a parameter to reflect wave transformation and the deviation from the Rayleigh distribution in shallow water.
* **Limitations:** Requires accurate knowledge of **$H_{2\%}$**, which is difficult to estimate without numerical or physical modeling.

### [4] Van Gent Simplified (2003)
The **Van Gent Simplified (2003)** formula was developed as a pragmatic design tool specifically tailored for coastal structures situated on **shallow foreshores**, where complex wave transformations render standard deep-water assumptions invalid. Unlike the comprehensive Van der Meer equations which rely on deep-water spectral characteristics, this formulation addresses the specific hydraulic reality of the **surf zone** where waves have already undergone significant shoaling and breaking before reaching the structure. It was derived from an extensive series of physical model tests focused on verifying the validity of stability predictions in conditions where the ratio of water depth to wave height is low, providing a specialized alternative for designs where the foreshore bathymetry plays a dominant role in filtering the incoming wave energy.

A central tenet of this formulation is the **deliberate omission of the wave period parameter** from the stability calculation. Van Gent and colleagues observed that in shallow water conditions, particularly where heavy breaking occurs on the foreshore, the wave energy spectrum becomes significantly flattened and broadband due to non-linear triad interactions. Statistical analysis of the experimental data indicated that the influence of the spectral period on rock stability became secondary to the dominant influence of the breaking wave height and the structure's geometry. This simplification eliminates the significant epistemic uncertainty often associated with defining a "correct" characteristic wave period in the surf zone, where spectra are often multi-peaked or severely distorted.

To further reduce subjectivity in the design process, the formula replaces the abstract "Notional Permeability" factor used in the Van der Meer equations with a measurable, physical geometric parameter: the **ratio of the nominal diameter of the core material to the nominal diameter of the armour layer**. By explicitly linking stability to the physical void structure defined by the filter and core sizing, the formula provides a more objective assessment of energy dissipation within the mound. It mathematically recognizes that a coarser core (a higher diameter ratio) allows for greater hydraulic penetration and pressure relief during wave drawdown, thereby increasing the stability of the outer armor layer against extraction forces.

While this simplification makes the formula highly robust for shallow water design where spectral shapes are ambiguous, it entails specific limitations that engineers must recognize. By excluding the wave period, the formula cannot explicitly account for the increased destructiveness of **long-period infragravity waves (surf beat)**, which often contain significant energy in the surf zone and can destabilize structures during the run-down phase. Consequently, while it is an excellent fallback for shallow foreshores, it tends to yield conservative results when applied to deep-water structures where the nuanced interaction of wave period and permeability—better captured by the full Van der Meer equation—plays a governing role in stability.

**Equation:**
$$\frac{H_{s}}{\Delta D_{n50}} = 1.75 \sqrt{\cot \alpha} \left(1 + \frac{D_{n50,core}}{D_{n50,arm}}\right) \left(\frac{S}{\sqrt{N}}\right)^{0.2}$$

**Explanation:**
* **Purpose:** To provide a robust, simplified formula for shallow foreshores where the wave spectrum is so flattened by breaking that defining a characteristic period becomes ambiguous or statistically insignificant regarding stability scatter.
* **Capabilities:** It eliminates the uncertainty of selecting the "correct" wave period in the surf zone. It introduces a physical permeability term based on the **ratio of core stone size to armor stone size ($D_{n50,core}/D_{n50,arm}$)**, which is more objective than the notional $P$.
* **Limitations:** By removing the wave period, it fails to account for the destructive power of **long-period infragravity waves**, which can be significant in shallow water. It tends to be conservative in deep water.

### [5] Eldrup & Andersen (2019)
The **Eldrup and Andersen (2019)** formulation emerged from the necessity to rectify stability predictions in **shallow water regimes where wave nonlinearity** becomes a governing factor. Traditional formulas, largely calibrated on deep-water Rayleigh distributions, often fail to capture the specific hydrodynamics of shoaling waves that have developed peaked crests and flat troughs (cnoidal or Stokes V profiles) but have not yet fully broken. This formulation was specifically developed to extend the validity of rock armour stability assessment into these highly nonlinear zones, particularly focusing on conditions characterized by low wave steepness—such as long-period swells entering shallow coastal waters—where standard methods were observed to yield inconsistent safety margins.

A distinct capability of this model is its direct utilization of the **spectral significant wave height** rather than time-domain parameters, ensuring a more stable descriptor of energy in distorted wave fields. The formula introduces a refined handling of the permeability influence and the surf similarity parameter, adjusting the exponents to better reflect the physical reality of wave run-up and run-down on permeable slopes under nonlinear loading. By differentiating the stability response between plunging and surging waves with coefficients specifically tuned for these shallow, non-breaking to pre-breaking conditions (**4.5 for plunging and 3.1 for surging**), it offers a more nuanced prediction of damage accumulation than the monolithic approaches of the past. It effectively bridges the gap between deep-water linear wave theory and the chaotic surf zone, providing a reliable design tool for the **"shoaling zone"** where wave asymmetry is most pronounced.

Despite its advancements in modeling nonlinearity, the formula is not without limitations. The experimental dataset upon which it was derived did not extensively cover the most extreme surf zone conditions where the relative depth drops below 1.5. Consequently, its application to the swash zone or conditions of intense, saturated breaking involves **extrapolation** and should be approached with caution. Furthermore, while it improves upon deep-water formulas for swell conditions, it requires accurate estimation of the local spectral wave height, which itself can be challenging to predict without advanced numerical wave transformation models in complex bathymetries. Therefore, it serves best as a specialized tool for shallow, non-breaking to lightly breaking regimes rather than a universal solution for all nearshore depths.

**Plunging:**
$$\frac{H_{m0}}{\Delta D_{n50}} = 4.5 \left(\frac{S}{\sqrt{N}}\right)^{0.2} 1.6^P \cdot \xi_{m-1,0}^{(0.4P - 0.67)}$$

**Surging:**
$$\frac{H_{m0}}{\Delta D_{n50}} = 3.1 \left(\frac{S}{\sqrt{N}}\right)^{0.2} P^{0.17} \cdot \min(\cot \alpha, 2)^{0.23}$$

**Explanation:**
* **Purpose:** To correct the over-prediction or under-prediction of stability when applying standard deep-water formulas to **nonlinear, shallow-water waves (cnoidal/Stokes V waves)**.
* **Capabilities:** It utilizes $H_{m0}$ directly. It provides a nuanced transition between breaking types that better fits low-steepness swell waves often encountered in coastal zones.
* **Limitations:** The experimental dataset did not cover the most extreme surf zone conditions ($h/H_{m0} < 1.5$), meaning it is an **extrapolation when applied to the swash zone**.

### [6] Etemad-Shahidi (2020)
The **Etemad-Shahidi (2020)** formulation represents a significant **data-driven** advancement in coastal engineering, aiming to bridge the predictive gap between deep-water stability and shallow-water transformation with a single, **unified model**. Derived using advanced M5' model tree data mining techniques on a massive dataset of 1,228 historical records, this formula was specifically calibrated to overcome the inconsistencies found when switching between disparate methods for different depth regimes. A defining feature of this approach is the replacement of the subjective "Notional Permeability" factor with a more rigorous **physical parameter**. This coefficient is directly calculated from the **ratio of the core's nominal diameter to the armor's nominal diameter**, providing an objective measure of hydraulic conductivity that eliminates the "engineering judgement" often required to estimate in the Van der Meer equations.

Furthermore, the model introduces a critical **foreshore slope correction factor**, representing the **tangent of the foreshore slope angle**. This term explicitly accounts for the wave transformation processes—such as shoaling and pre-breaking energy dissipation—that occur before the wave reaches the structure's toe. By integrating the foreshore geometry directly into the stability equation, the formula acknowledges that a wave traveling over a gentle slope (e.g., 1:100) will have a different spectral shape and destructive potential than a wave traveling over a steep slope (e.g., 1:10) or deep water, even if their height at the toe is identical. This makes it particularly robust for complex bathymetries where the foreshore influence is pronounced. While it offers high statistical accuracy across a broad range of conditions, recent evaluations by Scaravaglione et al. (2025) suggest that its standard coefficients may still **overestimate stability** in the **highly turbulent bore-flow** of the surf zone, necessitating the modifications seen in Formula [8] for very shallow waters.

**Plunging:**
$$N_s = 4.5 \cdot C_p \cdot N^{-0.1} \cdot S^{1/6} \cdot \xi_{m-1,0}^{-7/12} \cdot (1-3m)$$

**Surging:**
$$N_s = 3.9 \cdot C_p \cdot N^{-0.1} \cdot S^{1/6} \cdot \xi_{m-1,0}^{-1/3} \cdot (1-3m)$$

**Explanation:**
* **Purpose:** To provide a single, unified formula derived from data mining (M5 model tree) a vast dataset of historical tests, bridging deep and shallow water contexts.
* **Capabilities:** It replaces the notional permeability $P$ with a measurable physical coefficient $C_p$. Crucially, it includes the parameter $m$ (tangent of the *foreshore* slope, i.e., $m = 1.0 / \mathrm{foreshore\_cot}$), allowing it to account for wave transformation and shoaling effects occurring *before* the wave reaches the structure.
* **Limitations:** While statistically powerful, it lacks the specific calibration for the highly turbulent bore-flow of the surf zone, leading to **potential overestimation of stability in very shallow water** as shown by Scaravaglione (2025).

### [7] Scaravaglione Mod. VG (2025)
The **Scaravaglione Modified VG (2025)** formulation represents a specialized calibration of the original Van Gent equation, developed to address the catastrophic loss of stability often observed in the **extremely shallow water and swash zones**. In these hydrodynamic regimes, the incident waves have fully transitioned into **turbulent bores**, characterized by high-velocity fronts and significant **air entrainment**. Standard stability models, including the original Van Gent formula, typically rely on the assumption of hydrostatic buoyancy to counteract the gravity forces holding the armour units in place. However, in the swash zone, the heavy aeration of the water column significantly reduces the fluid density, thereby diminishing the buoyant force and increasing the effective weight of the stones that the hydrodynamic drag must overcome. Consequently, formulas calibrated for deeper water tend to dangerously underpredict the required armour mass in this zone, leading to unexpected failures during storm surges or high tide events.

To rectify this unsafe bias, Scaravaglione et al. (2025) proposed a radical adjustment to the stability coefficient, **increasing it from the original 1.75 to 3.3**. This substantial increase reflects the heightened instability caused by the combination of reduced buoyancy and the impulsive nature of the bore impact. Furthermore, while the original Van Gent formula deliberately removed the wave period to simplify design in the surf zone, this modification re-introduces a **dependence on wave steepness** raised to the power of 0.1. This term accounts for the physical reality that the momentum flux of the bore—and thus its capacity to dislodge armour units—remains weakly coupled to the incident wave period even after breaking. By reintegrating the spectral period influence, the formula provides a more precise description of the energy flux in the swash zone compared to purely depth-limited approximations.

The application of this formula is strictly bounded by the hydraulic conditions for which it was derived. It is exclusively valid for the **"Very Shallow" to "Swash" zones where the relative water depth is less than 0.5**. In these conditions, the structure is subjected to broken waves that behave more like channel flow than oscillatory waves. Applying this formula in deep or intermediate water would result in excessively conservative designs, yielding armour stones significantly larger than necessary due to the high coefficient. Therefore, it serves as a critical **"niche" tool** in the calculator's intelligent selection engine, activated only when the software detects that the structure is located in the most landward, high-risk zone of the coastal profile where standard physics no longer apply.

**Equation:**
$$\frac{H_{m0}}{\Delta D_{n50}} = 3.3 \sqrt{\cot \alpha} \left(1 + \frac{D_{n50,core}}{D_{n50}}\right) s_{m-1,0}^{0.1} \left(\frac{S}{\sqrt{N}}\right)^{0.2}$$

**Explanation:**
* **Purpose:** To address the specific physics of the **swash zone ($h/H_{m0} < 0.5$)** where air entrainment reduces buoyancy and flow becomes a high-velocity turbulent bore.
* **Capabilities:** It significantly increases the **stability coefficient (from 1.75 in Van Gent to 3.3)**, reflecting the reduced stability in this zone. It re-introduces a dependence on wave steepness ($s_{m-1,0}$) to better model the bore dynamics.
* **Limitations:** It is strictly a **"niche" formula for extremely shallow water**. Applying it in deep water would yield unrealistic armor units weights.

### [8] Scaravaglione Mod. ES (2025)
The **Scaravaglione Modified ES (2025)** formulation is a critical refinement of the Etemad-Shahidi unified model, specifically engineered to correct the systematic overestimation of stability observed in the **surf zone where waves are fully broken**. While the original Etemad-Shahidi formula successfully unified deep and shallow water predictions using a large dataset, recent experimental validation revealed that its standard coefficients tend to underpredict the damage caused by the highly turbulent, bore-like flow characteristic of the surf zone. In this regime, the incident wave energy has already been saturated by breaking, creating a **horizontal momentum flux** that is distinct from the orbital motion of shoaling waves. To address this, Scaravaglione et al. (2025) recalibrated the formula using a specialized dataset of broken waves, adjusting the **stability coefficient to 3.55** to better reflect the increased destructive potential of the bore.

A pivotal modification in this formula is the **decoupling of the wave steepness term from the Iribarren number**. In standard stability formulas, the surf similarity parameter combines slope angle and wave steepness into a single governing variable. However, in the surf zone where waves have already broken, the physical meaning of the Iribarren number—which predicts breaker type—diminishes. The Scaravaglione modification separates these influences, assigning a specific, reduced exponent of 0.05 to the wave steepness term. This mathematical adjustment acknowledges that once a wave has broken into a bore, its steepness has a much weaker influence on stability compared to the sheer mass and velocity of the water front. By isolating the steepness effect, the formula prevents the artificial inflation of stability that can occur when applying deep-water steepness relationships to broken waves.

This formulation is strictly intended for the specific hydraulic window of the surf zone, defined by a **relative water depth between 0.5 and 1.5**. It serves as a necessary bridge between the shoaling zone (covered by Van der Meer and Eldrup & Andersen) and the swash zone (covered by Scaravaglione Mod. VG). Applying this formula outside of this "broken wave" regime could lead to inaccuracies; for instance, in deep water, the reduced dependency on wave steepness might fail to capture the sensitivity of plunging breakers. Therefore, within the calculator's logic, this formula is automatically selected only when the intelligent design engine detects that the structure is located within the surf zone, ensuring that the stability assessment is driven by the physics of broken waves rather than generalized deep-water assumptions.

**Equation (Surging/Bore):**
$$\frac{H_{m0}}{\Delta D_{n50}} = 3.55 \cdot C_p \cdot N^{-1/10} \cdot (\cot \alpha)^{1/3} \cdot S^{1/6} \cdot s_{m-1,0}^{1/20}$$

**Explanation:**
* **Purpose:** To correct the Etemad-Shahidi formula for conditions where **waves are already broken ($0.5 < h/H_{m0} < 1.5$)**.
* **Capabilities:** It decouples the wave steepness term from the Iribarren number. In broken waves, the surf similarity parameter loses some physical meaning; this formula relies on specific coefficients (3.55) and steepness powers (0.05) calibrated to **preventing over-prediction of stability in the surf zone**.
* **Limitations:** Calibrated specifically for broken waves in the surf zone; not intended for deep water applications.

---

Here is the updated Markdown documentation with the expanded tables and formula definitions.

---

## 7. Armourstone Layers Design

Once the required nominal diameter  is calculated using the stability formulas above, the calculator aids in designing the physical armor layer geometry and specification. This process is governed by **EN 13383-1:2013**.

### Standard Rock Gradings (EN 13383)

The software includes a comprehensive library of standard rock gradings. The calculator automatically scans these lists to find a grading that strictly contains the required mass.

**Selection Criteria:**
The logic selects a standard grading class based on the calculated target mass ($M_{50,target}$) using the following rules:

1. **Strict Containment:** The target mass must fall strictly between the **Nominal Lower Limit (NLL)** and **Nominal Upper Limit (NUL)** of the class ($NLL < M_{50,target} < NUL$).
2. **Optimization:** If multiple classes satisfy the containment rule, the calculator selects the grading with the **narrowest range** (smallest difference between $NUL$ and $NLL$) to ensure the most specific fit.

**Limit Definitions & Formulas:**
For the selected grading, the software calculates the full range of mass limits used for quality control based on the standard EN 13383 definitions:

* **Nominal Lower Limit (NLL):** The standard lower mass boundary defined by the category name.
* **Nominal Upper Limit (NUL):** The standard upper mass boundary defined by the category name.
* **Representative $M_{50}$:** Calculated as $0.5 \times (NUL + NLL)$.
* **Extreme Lower Limit (ELL):** Calculated as $0.7 \times NLL$.
* **Extreme Upper Limit (EUL):** Calculated as $1.5 \times NUL$.

#### 1. Heavy Mass Armourstone (HMA)

These are used for large breakwaters and revetments where individual stone mass is the governing stability parameter.

| Grading Class | NLL (kg) | NUL (kg) | ELL (kg) | EUL (kg) |
| --- | --- | --- | --- | --- |
| **HMA 300-1000** | 300 | 1000 | 210 | 1500 |
| **HMA 1000-3000** | 1000 | 3000 | 700 | 4500 |
| **HMA 3000-6000** | 3000 | 6000 | 2100 | 9000 |
| **HMA 6000-10000** | 6000 | 10000 | 4200 | 15000 |
| **HMA 10000-15000** | 10000 | 15000 | 7000 | 22500 |

#### 2. Light Mass Armourstone (LMA)

These are standard gradings for smaller hydraulic structures or filter layers.

| Grading Class | NLL (kg) | NUL (kg) | ELL (kg) | EUL (kg) |
| --- | --- | --- | --- | --- |
| **LMA 5-40** | 5 | 40 | 3.5 | 60 |
| **LMA 10-60** | 10 | 60 | 7 | 90 |
| **LMA 15-120** | 15 | 120 | 10.5 | 180 |
| **LMA 40-200** | 40 | 200 | 28 | 300 |
| **LMA 60-300** | 60 | 300 | 42 | 450 |
| **LMA 15-300** | 15 | 300 | 10.5 | 450 |

#### 3. Coarse Gradings (CP - Coarse Particle)

These gradings are generally specified by sieve size (mm) but have defined mass characteristics used here as nominal limits for calculation.

| Grading Class | NLL (kg) | NUL (kg) | ELL (kg) | EUL (kg) |
| --- | --- | --- | --- | --- |
| **CP 45/125** | 0.4 | 1.2 | 0.28 | 1.8 |
| **CP 63/180** | 1.2 | 3.8 | 0.84 | 5.7 |
| **CP 90/250** | 3.1 | 9.3 | 2.17 | 13.95 |
| **CP 45/180** | 0.4 | 1.2 | 0.28 | 1.8 |
| **CP 90/180** | 2.1 | 2.8 | 1.47 | 4.2 |

### Custom Power Law Grading

For projects requiring non-standard rock sizes, or when the target mass falls outside all standard EN 13383 nominal limits (e.g., very large armor >15t or specific quarry yields), the calculator switches to a **Custom Grading** mode.

In this mode, the software generates theoretical mass limits based on the calculated target M50. It uses a power law distribution to estimate a theoretical  and  that would provide stability equivalent to the target mass, ensuring a constructible layer thickness is still calculated.

### Layer Geometry & Construction Formulas

The physical geometry of the armor layer is calculated based on the selected rock grading (Standard or Custom).

**1. Nominal Diameter ($D_{n50}$):**
Derived from the mass and rock density ($\rho_r$):

$$D_{n50} = \left(\frac{M_{50}}{\rho_r}\right)^{1/3}$$

**2. Layer Thickness ($t$):**
The perpendicular thickness of the armor layer (usually double layer, $n=2$):

$$t = n \cdot k_t \cdot D_{n50}$$

* $n$: Number of layers (Standard = 2).
* $k_t$: Layer thickness coefficient.
*  $k_t \approx 1.0$ for angular/blocky rock.
*  $k_t$ varies slightly with placement method (0.9 to 1.15 in EN 13383 standards depending on shape).

**3. Packing Density ($\phi$):**
This parameter is critical for construction quality control, defining how many stones are required per square meter of slope.

$$N_{units}/m^2 = n \cdot k_t \cdot (1 - n_v) \cdot D_{n50}^{-2}$$

* $n_v$: Volumetric porosity of the armor layer.
* Standard placement porosity is typically **30%** ($n_v = 0.30$).
* This implies that roughly 70% of the layer volume is solid rock, and 30% is void space for water dissipation.

This detailed breakdown ensures that the theoretical stability mass calculated by the formulas is translated into a practical, constructible engineering specification compliant with European standards.

---

## 8. Advanced Factors: Rock Shape and Placement

The stability of the armour layer is also a function of the inter-particle friction and interlocking, which are determined by the shape of the rocks and their placement method. The calculator code can be changed to integrate modification factors ($c_{pl}, c_{su}$) based on re-analysis of HR Wallingford datasets.

### Rock Shape Coefficients (Latham et al., 1988)
The standard formula assumes "fresh, rough, angular" rock. There are different lithologies or wear states:
* **Standard:** $c_{pl} = 1.0, c_{su} = 1.0$.
* **Rounded/Worn:** (e.g., river stone or eroded quarry run). Interlocking is reduced.
    * $c_{pl} = 0.90$ (10% reduction in plunging stability).
    * $c_{su} = 0.75$ (25% reduction in surging stability).
    * **Insight:** The reduction is more severe in surging waves because surging forces rely heavily on friction to prevent sliding, which rounded rocks lack.
* **Tabular:** (Flat rocks). If placed randomly, they can align with flow and lift. If placed flat, they are stable. The data suggests a slight increase in stability is permissible:
    * $c_{pl} = 1.10, c_{su} = 1.10$.

### Placement and Brittle Failure (Stewart et al., 2003)
The calculator code could also be updated to allow selection between **Random (Bulk) Placement** and **Dense/Interlocking Placement**.

* **Random Placement:** Follows the gradual damage law $S_d \propto H^5$. The calculator outputs a continuous damage curve.
* **Dense Placement:** Stones are oriented to maximize contact points.
    * **Effect:** Stability increases dramatically. The layer can withstand waves 1.5x to 2.0x higher than random placement without damage.
    * **Risk:** The failure mode becomes **Brittle**. There is no gradual accumulation of damage ($S_d=2 \to 4 \to 6$). Instead, the slope remains intact ($S_d \approx 0$) up to a high threshold, and then collapses completely ($S_d > 15$).

**Calculator Logic:** If Dense Placement is selected, the calculator does **not** simply output a smaller rock size. Instead, it calculates the failure threshold and applies a **Safety Factor of 1.5** to determining the design wave height. This ensures the structure operates in the "zero damage" zone, avoiding the catastrophic brittle failure region.

---

## 9. Variables glossary

* **$A_e$**: Cross-sectional erosion area of the profile.
* **$\alpha$**: Structure slope angle (angle with the horizontal).
* **$B_{Lc}$**: Blockiness coefficient of the armour stone units.
* **$C_p$**: Physical permeability coefficient (Etemad-Shahidi), defined as $C_p = [1 + (D_{n50,core}/D_{n50})^{0.3}]^{0.6}$.
* **$c_{pl}$**: Constructive coefficient for plunging waves (Van der Meer).
* **$c_{su}$**: Constructive coefficient for surging waves (Van der Meer).
* **$c_{VG}$**: Stability coefficient in the Van Gent formula.
* **$\cot\alpha$**: Cotangent of the structure seaward slope angle.
* **$\Delta$**: Relative buoyant density of the rock, defined as $\Delta = (\rho_r/\rho_w) - 1$.
* **$D_{15}$**: Nominal diameter exceeded by 15% of the sample by mass.
* **$D_{85}$**: Nominal diameter exceeded by 85% of the sample by mass.
* **$D_{n50}$**: Nominal diameter of armour rock, defined as $D_{n50} = (M_{50}/\rho_r)^{1/3}$.
* **$D_{n50,core}$**: Nominal diameter of the core material.
* **$g$**: Gravitational acceleration ($9.80665 m/s^2$).
* **$h$** or **$h_{toe}$**: Water depth at the toe of the structure.
* **$h_r$**: Relative water depth, defined as $H_{m0,deep}/h$ (Note: Scaravaglione uses this definition; calculator uses $R_h = h/H_{m0}$).
* **$H_{1/3}$**: Significant wave height derived from time-domain analysis (average of the highest 1/3 of waves).
* **$H_{1/10}$**: Average wave height of the highest 10% of waves.
* **$H_{1/50}$**: Average wave height of the highest 50 waves in the record (Vidal).
* **$H_{2\%}$**: Wave height exceeded by 2% of the waves in the record.
* **$H_{50}$**: Average wave height of the 50 highest waves (Etemad-Shahidi/Vidal notation).
* **$H_{m0}$**: Spectral significant wave height, defined as $4\sqrt{m_0}$.
* **$H_{m0,deep}$**: Spectral significant wave height in deep water.
* **$H_{s}$**: Significant wave height (general term, often used interchangeably with $H_{1/3}$ or $H_{m0}$ depending on context).
* **$K_D$**: Stability coefficient used in the Hudson formula.
* **$k_t$**: Layer thickness coefficient.
* **$LT$**: Length-to-thickness ratio of the armour stone units.
* **$m$**: Slope of the foreshore (tangent). Often used in the form $(1-3m)$.
* **$m_n$**: Spectral moments of the n-th order ($n = -1, 0, \dots$).
* **$M_{50}$**: Median mass of the armour rock grading (50% value by sample mass).
* **$N$** or **$N_w$**: Number of waves in the storm duration.
* **$N_{od}$**: Number of displaced units (damage parameter used by Thompson & Shuttler).
* **$N_s$**: Stability number, defined as $H_{s} / (\Delta D_{n50})$.
* **$n_v$**: Volumetric porosity of the armor layer.
* **$P$**: Notional permeability factor (Van der Meer).
* **$\Pi_0$**: Nonlinearity parameter (Goda).
* **$\phi$**: Packing density (units per $m^2$).
* **$R_c$**: Crest freeboard (height of crest above still water level).
* **$R_h$**: Relative depth ratio used in the calculator logic ($h_{toe}/H_{m0,deep}$).
* **$R_{eA}$**: Reynolds number at the armour layer.
* **$\rho_r$**: Mass density of the rock.
* **$\rho_w$**: Mass density of the water.
* **$S_d$** or **$S$**: Non-dimensional damage level ($A_e/D_{n50}^2$).
* **$s_{m-1,0}$**: Wave steepness based on $T_{m-1,0}$ at the toe of the structure ($2\pi H_{m0} / (g T_{m-1,0}^2)$).
* **$s_{m-1,0,o}$**: Deep water (offshore) spectral wave steepness.
* **$s_{om}$**: Deep water wave steepness based on mean period $T_m$.
* **$t$**: Perpendicular layer thickness.
* **$T_m$**: Mean wave period from time-domain analysis.
* **$T_{m-1,0}$**: Spectral energy wave period ($m_{-1}/m_0$).
* **$T_p$**: Peak wave period.
* **$\xi_m$**: Surf similarity parameter (Iribarren number) based on mean period $T_m$.
* **$\xi_{m-1,0}$**: Surf similarity parameter based on spectral period $T_{m-1,0}$.
* **$\xi_{cr}$**: Critical surf similarity parameter defining the transition between plunging and surging waves.

---

## 10. Bibliography & Scientific References

**1. U.S. Army Corps of Engineers. (1984).** *Shore Protection Manual. Vol. I & II.* Coastal Engineering Research Center, Vicksburg, MS.  
[https://usace.contentdm.oclc.org/digital/collection/p16021coll11/id/1934/](https://usace.contentdm.oclc.org/digital/collection/p16021coll11/id/1934/)

**2. Van der Meer, J. W. (1987).** "Stability of breakwater armour layers—Design formulae." *Coastal Engineering*, 11(3), 219-239.  
[https://doi.org/10.1016/0378-3839(87)90013-5](https://doi.org/10.1016/0378-3839(87)90013-5)

**3. Van der Meer, J. W. (1988).** *Rock Slopes and Gravel Beaches Under Wave Attack.* Doctoral Thesis, Delft University of Technology.  
[https://repository.tudelft.nl/islandora/object/uuid:404b5nec-hm](https://repository.tudelft.nl/islandora/object/uuid:404b5nec-hm)

**4. Van Gent, M. R. A. (1995).** *Wave interaction with permeable coastal structures.* Doctoral Thesis, Delft University of Technology.  
[https://repository.tudelft.nl/islandora/object/uuid:7bbff8e4-215d-4bfc-a3af-51cdecb754bd](https://repository.tudelft.nl/islandora/object/uuid:7bbff8e4-215d-4bfc-a3af-51cdecb754bd)

**5. U.S. Army Corps of Engineers (USACE). (2002).** "Coastal Engineering Manual." Engineer Manual 1110-2-1100, Washington, D.C.  
[https://www.publications.usace.army.mil/USACE-Publications/Engineer-Manuals/](https://www.publications.usace.army.mil/USACE-Publications/Engineer-Manuals/)

**6. Van Gent, M. R. A., Smale, A. J., & Kuiper, C. (2003).** "Stability of rock slopes with shallow foreshores." *Proceedings of Coastal Structures 2003*, Portland, OR, 100-112.  
[https://doi.org/10.1061/40733(147)9](https://doi.org/10.1061/40733(147)9)

**7. Van Gent, M. R. A. (2004).** "On the stability of rock slopes." *Environmentally Friendly Coastal Protection: Proceedings of the NATO Advanced Research Workshop*, Varna, Bulgaria.  
[https://doi.org/10.1007/1-4020-3301-X_12](https://doi.org/10.1007/1-4020-3301-X_12)

**8. CIRIA, CUR, CETMEF. (2007).** *The Rock Manual. The Use of Rock in Hydraulic Engineering.* (2nd edition). C683, CIRIA, London.  
[https://www.ciria.org/ItemDetail?iProductCode=C683](https://www.ciria.org/ItemDetail?iProductCode=C683)

**9. CEN (2013).** "EN 13383-1:2013 Armourstone - Part 1: Specification." European Committee for Standardization.  
[https://standards.iteh.ai/catalog/standards/cen/5f6b770f-38ba-4320-ad1a-83d0e29f7db2/en-13383-1-2013](https://standards.iteh.ai/catalog/standards/cen/5f6b770f-38ba-4320-ad1a-83d0e29f7db2/en-13383-1-2013)

**10. Eldrup, M. R., & Lykke Andersen, T. (2019).** "Extension of shallow water rock armour stability formulae to nonlinear waves." *Coastal Engineering*, 153, 103536.  
[https://doi.org/10.1016/j.coastaleng.2019.103536](https://doi.org/10.1016/j.coastaleng.2019.103536)

**11. Etemad-Shahidi, A., Bali, M., & Van Gent, M. R. A. (2020).** "On the stability of rock armored rubble mound structures." *Coastal Engineering*, 158, 103655.  
[https://doi.org/10.1016/j.coastaleng.2020.103655](https://doi.org/10.1016/j.coastaleng.2020.103655)

**12. Van der Meer, J. W. (2021).** "Rock armour slope stability under wave attack; the Van der Meer Formula revisited." *Journal of Coastal and Hydraulic Structures*, 1.  
[https://doi.org/10.48438/jchs.2021.0008](https://doi.org/10.48438/jchs.2021.0008)

**13. Van der Meer, J. W., Lykke Andersen, T., & Roge Eldrup, M. (2024).** "Rock Armour Slope Stability under Wave Attack in Shallow Water." *Journal of Coastal and Hydraulic Structures*, 4.  
[https://doi.org/10.59490/jchs.2024.0035](https://doi.org/10.59490/jchs.2024.0035)

**14. Scaravaglione, G., Marino, S., Francone, A., Leone, E., Damiani, L., Tomasicchio, G. R., Van Gent, M. R. A., & Saponieri, A. (2025).** "The influence of shallow water on rock armour stability." *Coastal Engineering*, 197, 104657.  
[https://doi.org/10.1016/j.coastaleng.2024.104657](https://doi.org/10.1016/j.coastaleng.2024.104657)