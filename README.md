_JMJ-V!_

# idealgasentropy
Entropie (und etwas mehr) eines idealen (nicht perfekten!) H+He-Gases wie in D'Angelo & Bodenheimer (2013) ausrechnen
|| Compute the entropy (and a bit more) of an ideal (but not perfect!) H+He gas as in D'Angelo & Bodenheimer (2013)

v.1: 31.10.2019 (c) Gabriel-Dominique Marleau, Uni Tübingen

## Main functions
Total specific entropy (per unit mass) _s_ of the mixture, in units of k_B/baryon:
- stot_proMasse_rhoT(Y, rho, T)        : the main entropy function, but...
- stot_proMasse_rhoT_bystro(Y, rho, T) : ... this is the same but should be faster because repetitions are avoided

More functions:
- deladDAB13_rhoT(Y, rho, T)           : adiabatic gradient (dlnT/dlnP)\_{const s}
- muDAB13(Y, rho, T)                   : mean molecular weight, dimensionless
 
### Conventions
- rho: density,         always in g/cm^3
- T:   temperature,     always in K
- Y:   helium fraction, always dimensionless

## Notes
- Ortho:Para ratio is not treated explicitly: Only the limit T >> T_rot ~ 85 K is implemented currently
- No metals. If there should be metals, adding them to helium offers an approximate treatment.

## Comments
- For most of the relevant rho--T plane, delad from here and from D'Angelo & Bodenheimer (2013), using their functions as implemented in Pluto (Vaidya et al. 2015), agree to better than 5%. However, there are some ~20% differences

- At the lowest pressure of Saumon, Chabrier & van Horn (1995), the gas should be ideal yet there are small piecewise-constant offsets in the entropy. This could be due to a mistake on my part (despite checking and despite the relative simplicity of the equations) or to non-ideal (interaction, not degeneracy) effects being important in SCvH. Their approach is very different from the simple non-perfect ideal gas.

- Some documentation needs to be done... The entropy formula is an "original derivation"; I have not seen it elsewhere but have also not searched, and is relatively easy to derive.


All comments, questions, suggestions for improvements, etc. are welcome! Please write to me at uni-tuebingen.de with gabriel.marleau in front.
