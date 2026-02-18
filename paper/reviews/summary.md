# Review Summary

**Due date:** March 17th

**Decision:** Major revision

**Editor:** Florian Lemarie

Three reviews received. All supportive, one gave full thumbs up. Editor emphasizes the need to
clearly distinguish what is currently implemented vs. prospective, and to provide more guidance
on which numerical/physical choices to prioritize depending on intended application scale.

## Actionable comments (motivate new code or simulations in this repo)

1. **I/O performance** (R4 #3): I/O is never mentioned. It is a major bottleneck for ESMs on GPUs.
   If there is ongoing work, mention it; otherwise acknowledge the gap.

2. **Immersed boundary / bottom cells** (R4 #501): Clarify whether the model uses partial vs full
   bottom-cell representation. Document the "grid-fitted immersed boundary method" clearly.

3. **Non-hydrostatic reference simulation** (R4 #559): Use the non-hydrostatic kernel as a
   reference for the LES-capable simulations in section 3.2.1 (e.g. Figure 7).

4. **Non-linear free surface / z\* coordinate status** (R4 #777-786, Editor): Clarify whether
   z\* is functional. The paper cites Silvestri et al. (2025) for the free surface, but that only
   covers a linear free surface; Appendix A appears to describe a nonlinear free surface.

5. **Vertical advection treatment** (Editor): Not discussed anywhere; editor flags this as
   essential for a z-coordinate model.

6. **Lateral boundary conditions at coast** (Editor): Not discussed; editor flags as essential.

7. **Tripolar grid north-fold boundary condition** (R5): How is this handled? Not covered in
   Silvestri (2025).

8. **Simulation timing / SYPD** (R4 #473): Add wall-clock duration or SYPD for the examples shown.

## Opinion / framing comments (about claims and positioning)

1. **"Accelerate progress" claim not justified** (R4 #1): No proof that Oceananigans enables
   something that *cannot* be done with other models. Either provide a concrete example of
   unprecedented capability, or soften the claim.

2. **"Fastest ocean modelling engine" unqualified** (R5): Needs rigorous like-for-like comparison
   at comparable physical fidelity; probably requires community benchmarks.

3. **Performance claims overrepresented** (R4 #2): Performance is not the focus of this paper;
   reduce or rephrase repeated claims, especially in conclusions.

4. **GPU vs CPU node comparison fairness** (R4 #108-109): Is it fair to compare 1 GPU node to
   1 CPU node? No justification given in Silvestri et al. (2025) either.

5. **"Prototyped in different languages" is speculative** (R4 #129-130): Claim that new
   parameterizations are "typically prototyped" in different codes/languages is not universally true.

6. **Mathematical background purpose unclear** (R5): If novel, needs validation via comparison;
   if not novel, why include it?

7. **Paper scope / target audience** (Editor, R5): Unclear whether the paper is about the
   programming interface or validation. Editor suggests targeting current/future users with
   guidance on when to use which options.

## Paper writing comments (not relevant to code)

### Structural

- Separate future work from results into a single dedicated subsection (R4 #5)
- Clarify relationship between Oceananigans, ClimaOcean, ClimaSeaIce, Oceanostics (R4 #4)
- Introduce ClimaOcean before first use at line 614 (R4 #614)
- Introduce and define ClimaSeaIce with a reference (R4 #628)
- Provide more detail on forcing data mechanism (JRA55 download/caching) (R5)
- Describe how the interface can be extended when users need new capabilities (R5)
- Describe arithmetic / discrete calculus abstractions in more detail (R5)

### Line-level edits

- L74: "languages" -> "language" (R4)
- L98-113: Reframe performance paragraph; move to separate subsection pointing to Silvestri et al. (R4)
- L102: Misleading reference to section 3.2.3 (R4)
- L106-110: Use cores/tasks not "nodes"; verify the comparison numbers and cite properly (R4)
- L110: "argue" -> "argues" (R4)
- L146: Clarify what "pure Julia" means (R4)
- L210-211: Clarify why this can't be achieved with JAX (R4)
- L219: Documentation reference missing (R4)
- L242: "A second package" -- what is the first? (R4)
- L247-248: Name change from ClimaOcean to what? Clarify or remove (R4)
- L261: "Fundamentally" is confusing; rephrase (R4)
- L288-289: Simplify "expression tree of discrete calculus" (R4)
- L406-411: Anecdotal remark; remove or move to future work (R4)
- L412: Missing connection between listing 4 and the rest of the paragraph (R4)
- L451: Inconsistency between non-dimensional time (L286) and seconds here (R4)
- L454-460: Model weakness discussion; consider moving to future work (R4)
- L499: Figure 5 nearly uncommented in text; justify its presence (R4)
- L639-640: Remove "reduces the cost" claim from conclusion (not supported in this paper) (R4)
- L647-648: Specify how many cores (R4)

### Figures

- Fig 1: Use subfigure labels (a), (b) or refer to "left" / "right" in text (R4)
- Fig 1, 2: Missing legend and axis labels (R4)
- Fig 7: Legend should use acronyms from text/caption, not paper citations (R4)
