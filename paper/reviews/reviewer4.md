# Reviewer 4

> The present manuscript describes Oceananigans and its interface, giving examples on how to
> use its user interface. The paper is well written, within the scope of this journal, and
> useful to all those scientists that are willing to approach Oceananigans without knowing too
> much of the interface and the functioning of the code.

**Recommendation:** Publication after major revision.

## General Remarks

### 1. "Accelerate progress" claim not justified

> In the abstract and elsewhere, it is claimed that the enhanced interface and unmatched
> performance will ultimately accelerate progress in ocean and climate science. However,
> no proof is provided; the author compares the model to existing ones and demonstrates
> that similar results are produced more quickly. This will hardly accelerate progress.

Provide a concrete example of something achievable with Oceananigans that cannot be simulated
with other models (e.g. unprecedented resolution or physical processes). Otherwise, soften the
claim and present Oceananigans as an incredible pedagogical tool.

### 2. Reduce performance claims

> In many places (73, 90-113, 639-40, Conclusions) it is stated that Oceananigans enables
> efficient computation. This claim is well grounded (Silvestri et al. 2023, 2025), but it
> is not the focus of this paper.

Reduce or rephrase these claims, especially in conclusions -- address solely the results shown.

### 3. I/O not mentioned

> The author never mentions input/output. In Silvestri et al. (2025), I/O is disabled in
> performance measurements. I/O is one of the main bottlenecks for ESMs, especially on GPUs.
> This subject is not mentioned anywhere, not even in future work, which starts being suspicious.

If there is ongoing work on I/O, mention it. Otherwise, clearly state that the lack of a
competitive I/O component is a major flaw compared to traditional models.

### 4. Package relationships unclear

> The author mentions at least three other packages: Oceanostics ("companion package"),
> ClimaOcean ("coupled modeling package"), and ClimaSeaIce (never defined). These are not
> definitions.

Clearly state the relationship between the main library and these packages.

### 5. Separate future work from results

> Pretty much everywhere the author mixes results and future work, and then an additional
> "future work" section is in the conclusions.

Consolidate all future work into a single subsection.

## Detailed Review

### Typos and small fixes

- **L74:** "languages" -> "language"
- **L110:** "argue" -> "argues"

### Performance section (L98-113)

Reframe this paragraph. Move to a separate subsection where it is clear from the first sentence
that performance is discussed in Silvestri et al. (2025), then proceed with updated numbers.

- **L102:** Gives the impression computational results will appear in section 3.2.3, but they don't.
- **L106-110:** "Nodes" is a bad unit of measure (implies the reader knows cores per node).
  Always use cores/tasks. Also, the cited paper says OM4p25 reaches 12 SYPD on 4,671 cores,
  which doesn't simply extrapolate to the comparison made. Provide a more verifiable comparison
  or update the source.
- **L108-109:** Important comparison -- is it fair to compare 1 GPU node to 1 CPU node?
  What about energy consumption? (Same issue in Silvestri et al. 2025.)

### Specific line comments

- **L129-130:** "New parameterizations are typically prototyped [in different codes]" is purely
  speculative. Reviewer's experience: one designs a targeted test case, not necessarily in a
  different language.
- **L146:** What does "pure Julia" mean? Is Oceananigans not "pure" Julia?
- **L210-211:** Not clear why custom forcings/BCs can't be achieved with JAX. The cited paper
  doesn't clarify this.
- **L219:** Documentation reference missing.
- **L242:** "A second package" -- what is the first package?
- **L247-248:** What name change? From ClimaOcean to what? Is this relevant?
- **L261:** "Fundamentally" confuses the reader. Rephrase in more absolute terms.
- **L288-289:** "Expression tree of discrete calculus" -- simplify/clarify. Are you referring to
  the computation of derivatives and arithmetic on line 18, listing 1?
- **L406-411:** Anecdotal remark. Remove; mention in future work if desired.
- **L412:** Missing connection between listing 4 and the rest of the paragraph (similar to L443).
- **L451:** On L286 (same NonHydrostaticModel), time is non-dimensional; here it is in seconds.
  Is this due to the Simulation function (not shown in listing 6)?
- **L454-460:** Important to highlight model weakness, but consider moving to future work.
- **L473:** How long does the simulation take? Add wall-clock duration or SYPD.
- **L499:** Figure 5 is barely discussed in the text. Justify its presence.
- **L501:** How are solid boundaries treated? What is behind the "grid-fitted immersed boundary
  method"? Does the model use partial (vs full) bottom-cell representation of topography, with
  piecewise-constant topography? Clarify staircase vs shaved cells (L679-680).
- **L559 (Section 3.2.1):** Since Oceananigans can run LES, it would be convincing to use the
  non-hydrostatic kernel as a reference for the computations shown (e.g. Figure 7).
- **L614:** ClimaOcean appears suddenly. Clarify separation/connection with Oceananigans.
- **L628:** ClimaSeaIce is never introduced. Provide a reference and description.
- **L639-640:** "Reduces the cost of high-resolution simulations" -- numbers are shown but not
  justified (no scaling curve, performance analysis, I/O discussion). Remove from conclusions.
- **L647-648:** On how many cores?
- **L777-786:** Clarify whether the non-linear free surface (e.g. z\* coordinate) is functional.

### Figures

- **Fig 1:** Two pictures without subfigure labels. Use (a)/(b) or reference "left"/"right" in text.
- **Fig 1, 2:** Missing legend and axis labels.
- **Fig 7:** Legend is confusing. Use the acronyms from the text/caption; don't cite papers in
  the legend.
