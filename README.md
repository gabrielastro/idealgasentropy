_JMJ-V!_

# idealgasentropy
Entropie (und etwas mehr) eines idealen (nicht perfekten!) H+He-Gases wie in D'Angelo & Bodenheimer (2013) ausrechnen
|| Compute the entropy (and a bit more) of an ideal (but not perfect!) H+He gas as in D'Angelo & Bodenheimer (2013)

(c) Gabriel-Dominique Marleau, Uni Tübingen,
    with parts taken from [https://github.com/andrewcumming/gasgiant](gasgiant) from Andrew Cumming / David Berardo

v.1.0, 31.10.2019: Initial commit
v.1.1, 03.03.2023: Making a python package for easy import

## Main functions
Total specific entropy (per unit mass) _s_ of the mixture, in units of k_B/baryon:
- stot_proMasse_rhoT(Y, rho, T)        : the main entropy function, but...
- stot_proMasse_rhoT_bystro(Y, rho, T) : ... this is the same but should be fast (быстро) because repetitions are avoided

The same functions exist with ..._PT() to give (P,T) as arguments (simple wrappers)

More functions:
- deladDAB13_rhoT(Y, rho, T)           : adiabatic gradient (dlnT/dlnP)\_{const s}
- muDAB13(Y, rho, T)                   : mean molecular weight, dimensionless
- Dichte_PT(Y, P, T)                   : get rho from P and T; taken from David Berardo from https://github.com/andrewcumming/gasgiant and (very) adapted
 
### Conventions
- rho: density,              always in g/cm^3
- T:   temperature,          always in K
- Y:   helium mass fraction, always dimensionless

## Notes
- Ortho:Para ratio is not treated explicitly: Only the limit T >> T_rot ~ 85 K is currently implemented

- No metals. If there should be metals, adding them to helium offers an approximate treatment

- Written and tested in Python 2 only

## Comments
- For most of the relevant rho--T plane, delad from here and from D'Angelo & Bodenheimer (2013; hereafter DAB13), using their functions as implemented in Pluto (Vaidya et al. 2015), agree to better than 5%. However, there are some ~20% differences.

- At the lowest pressure of Saumon, Chabrier & van Horn (1995), the gas should be ideal yet there are small piecewise-constant offsets in the entropy. This could be due to a mistake on my part (despite checking and despite the relative simplicity of the equations) or to non-ideal (interaction, not degeneracy) effects being important in SCvH. Their approach is very different from the simple non-perfect ideal gas.

- Some documentation needs to be done... The entropy formula is an "original derivation"; I have not seen it elsewhere but have also not searched, and is relatively easy to derive.

## Comparisons
- There are some plots in Abb/ showing the goodness of the match of this implementation to DAB13. Why it is not perfect, is not clear.

- There are also some comparison to SCvH (see comments above).

All comments, questions, suggestions for improvements, etc. are welcome! Please write to me at uni-tuebingen.de with gabriel.marleau in front.
