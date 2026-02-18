# Editor Comments

**Editor:** Florian Lemarie

## Decision

Major revision. Three supportive reviews received (one full thumbs up).

> Please give this process the care needed to fully address the comments. Also, please give
> the reviewers and editors your thanks in the acknowledgment section.

## Key Points

1. The paper lies at the intersection of technical documentation, a survey of methodologies,
   and a perspective on future developments. The extremely broad scope makes it unclear what
   the reader should expect, and interesting points are addressed only superficially.

2. **Target audience appears to be current/future users.** Go beyond a catalog-like description
   and provide perspective on which options to prioritize depending on the intended application.
   A key feature is the ambition to "model motions from millimeters to millennia" -- provide a
   synthesis of numerical and physical choices depending on the scales being modeled, and for
   each regime, outline the corresponding perspectives. For example: when to use vector-invariant
   form vs flux form of the momentum equation.

3. **Distinguish implemented vs. prospective.** It is important to unambiguously distinguish what
   is currently available and properly tested (and under which circumstances) from what is
   prospective. For example, Silvestri et al. (2025) is cited for the free surface but only
   discusses a *linear* free surface (z-coordinate), while Appendix A appears to deal with a
   *nonlinear* free surface.

4. **Missing topics.** Covering so many topics raises the question of why certain elements are
   not discussed:
   - Treatment of vertical advection
   - Nature of lateral boundary conditions at the coast (essential in a z-coordinate model)
