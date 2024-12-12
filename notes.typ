#set page(height: auto)

#let Bi = "Bi"
#let Fo = "Fo"
#let erf = "erf"
#let erfc = "erfc"
#let Re = "Re"
#let Nu = "Nu"
#let Pr = "Pr"

= Definitions
/ Adiabatic Process: A process in which no heat transfer occurs
/ Isothermal: All points on a surface are the same temperature
/ Diffuse: No prefered direction for outgoing rays
/ Specular: Prefered direction for outgoing rays
/ Black Body: Perfect Emitter and Absorber ($epsilon = 1, alpha = 1$)
/ Gray Body: Properties are independent of wavelength
/ Irradiance (G): Flux of energy that irradiates a surface
/ Radiosity (J): total flux of radiative energy away from a surface
/ Reflectivity ($rho$): How much of the irradiance is reflected
/ Absorptibity ($alpha$): How much of the irradiance is absorbed
/ Transmissivity ($tau$): How much of the irradiance passes through the material

= When Assumptions Can Be Made
- Steady State: 
- 1D Conduction: When conduction is significant in only a single direction
- $tau = 0$: material is opaque or thick (>1 wavelength of light)

= Governing Equation

$ E_"in" - E_"out" + E_"gen" = E_"st" $

/*

CONDUCTION

*/
= Conduction
=== Heat Diffusion Equation
Cartesian Coordinates:
$ frac(partial, partial x) (k frac(partial T, partial x)) + frac(partial, partial y) (k frac(partial T, partial y)) + frac(partial, partial z) (k frac(partial T, partial z)) + dot(q)^''' = rho c_p frac(partial T, partial t) $

Cylindrical Coordinates:
$ 1/r frac(partial, partial r) (k r frac(partial T, partial r)) +  1/r^2 frac(partial, partial phi.alt) (k frac(partial T, partial phi.alt)) + frac(partial, partial z) (k frac(partial T, partial z)) + dot(q)^''' = rho c_p frac(partial T, partial t) $

Spherical Coordinates:
$ 1/r^2 frac(partial, partial r) (k r^2 frac(partial T, partial r)) + frac(1, r^2 sin^2(theta)) frac(partial, partial phi.alt) (k frac(partial T, partial phi.alt)) + frac(1, r^2 sin(theta)) frac(partial, partial phi.alt) (k sin theta frac(partial T, partial theta)) + dot(q)^''' = rho c_p frac(partial T, partial t) $

=== Boundary Conditions
1. Constant Surface Temperature: $T(0, t) = T_s$
2. Constant Surface Heat Flux \ 
    a. Finite Heat Flux: $-k frac(partial T, partial x) |_(x=0) = q_s^''$ \
    b. Adiabatic or Insulated Surface: $-k frac(partial T, partial x) |_(x=0) = 0$
3. Convection Surface Condtion: $-k frac(partial T, partial x) |_(x=0) = h[T_infinity - T(0, t)]$

=== Thermal Diffusivity
The ratio of how a material conducts thermal energy to how well it stores thermal energy.
$ alpha = frac(k, rho c_p) $

=== Circuit Analogy
Cartesian Coordinates:
$ R_"cond" =  frac(L, k A_s) $

Cylindrical Coordinates:
$ R_"cond" = frac(ln(r_2 / r_1), 2 pi k L) $

Spherical Coordinates:
$ R_"cond" = frac(1, 4 pi k) (1/r_1 - 1/r_2) $

=== Overall Heat Transfer Coefficient
Used for composite materials such as walls with different materials.
$ U = frac(1, R_"tot" A_s) $

=== Contact Resistance
To account for gaps due to surface roughness between mating surfaces.
$ R_"tc"^'' = frac(T_A - T_B, q_x^'') $

=== Porous Materials
These materials have pockets of liquid which greatly affects the thermal conductivity, so the average thermal conductivity is used to ease the calculates called $k_"eff"$.

Porosity (void fraction): $epsilon$
Fluid Thermal Conductivity: $k_f$
Solid Thermal Conductivity: $k_s$

$ k_("eff","min") = frac(1, (1-epsilon)/k_s + epsilon/k_f) $

$ k_("eff","max") = epsilon k_f + (1 - epsilon) k_s $

$ k_"eff" = [frac(k_f + 2 k_s - 2 epsilon (k_s - k_f), k_f + 2 k_s + epsilon (k_s - k_f) )] k_s $

=== Fins (Extended Surfaces)
(Page 118 has table of fin efficiencies for various shapes Basic Heat and Mass Transfer)

$ m = (frac(h P, k A_c))^(1/2) = (frac(2 h, k t))^(1/2) = (frac(4 h, k D))^(1/2) $

Fin Effectiveness: 
$ epsilon_f = (frac(k P, h A_c))^(1/2) $

Fin Efficiency:
$ eta_f = frac(tanh(m L), m L) $

Resistance:
$ R_(t,f) = frac(1, h A_f eta_f) $

== Transient Conduction
=== Biot Number
Ratio of thermal resistances between conduction and convection. If Bi << 1 then the conduction resistance is much less than the convective resistance and it is safe to assume uniform temperature distribution.

$ Bi equiv frac(accent(h, macron) L, k) $

=== Fourier Number
Dimensionless Time
$ Fo equiv frac(alpha t, L_c^2) $

=== Lumped Thermal Capacitance
Condition: $Bi < 0.1$ (error associated with method is negligible)

$ theta / theta_i = frac(T - T_infinity, T_i - T_infinity) = exp(-Bi dot Fo) $

=== Problem-Solving Strategy
(Page 209 Basic Heat and Mass Transfer)

+ Calculate Bi and if Bi < 0.1 then use Lumped Thermal Capacitance
+ Calculate Fourier Number. If Fo < 0.05 then Eqn 3.61 (pg. 194)
+ Calculate Fourier Number. If 0.05 < Fo < 0.2 then Eqn 3.72 and 3.73 (pg. 205)
+ Calculate Fourier Number. If Fo > 0.2 then Eqn 3.75 - 3.77 (pg. 206)

/*

CONVECTION

*/

= Convection
=== Reynold's Number
$ Re_x = frac(rho u_infinity x, mu) = frac(u_infinity x, nu) $

=== Nusselt Number
$ Nu equiv frac(h L, k_f) $

Flat plate in Laminar Flow:
$ Nu = 0.664 Re_L^(1/2) Pr^(1/3) $

=== Prandtl Number
$ Pr equiv nu / k $

=== Convective Heat Transfer Coefficient
$ h = -frac(k, T_w - T_infinity) frac(partial T, partial y) |_(y=0) $


/*

RADIATION

*/

= Radiation
=== Stefan-Boltzmann Law
$ d q = d U + P d VV $

Energy Density ($epsilon.alt$) \
Internal Energy ($U$) = $epsilon.alt VV$  \
Pressure ($P$) = $frac(1,3) epsilon.alt$ \

=== Spectral Energy Density
$ E(lambda) = frac(k_B T, lambda^4) $
$ E(f) = frac(k_B T, c^3) f^3 $

Boltzmann Const ($k_B$) = $1.38 times 10^(-23) frac(J, K)$

=== Total Radiative Energy
$ E = integral_0^infinity E(lambda) d lambda = integral_0^infinity E(f) d f $

=== Emissivity (Emittance)
The fraction of emissive power a real body ($E(T)$) emits compared to a black body ($E_b (T)$).
$ epsilon equiv frac(E (T), E_b (T)) $

=== Emissive Power
$ E(T) = epsilon sigma_"SB" T^4 $

=== Kirchoff's Law
A body in thermodynamic equilibrium must emit as much energy as it absorbs in each direction at each wavelength. This is to avoid violating the 2nd Law of Thermodynamics.
$ epsilon_lambda (T, theta, phi.alt) = alpha_lambda (T, theta, phi.alt) $

Diffuse Form:
$ epsilon_lambda (T) = alpha_lambda (T) $

Diffuse and Gray Form:
$ epsilon(T) = alpha(T) $

=== Wiens Displacement Relation
$ lambda_"peak" T = "const" = 2898 mu m K $

=== View Factors
The view factor is a correction to the Stefan-Boltzmann constant for black bodies. The view factor is a function of the surface area of two materials, the angle between the normals of each surface and the ray of radiation ($beta_1, beta_2$), and the distance between the two surfaces. (See page 561 in AHTT for table of view factors)

$ F_(1-2) = frac(1, A_1) integral_A_1 integral_A_2 frac(cos(beta_1) cos(beta_2), pi s^2) d A_2 d A_1 $

$ Q_"net"(1-2) = A_1 F_(1-2) sigma_"SB" (T_1^4 - T_2^4) $

$ A_1 F_(1-2) = A_2 F_(2-1) $

=== Irradiance
The flux of energy that irradiates a surface.
$ rho + alpha + tau = 1 $

=== Radiosity
The total flux of radiative energy away from a surface.
$ J = E + rho G  = epsilon E_b + rho G $

=== Circuit Analogy
Surface Resistance for a diffuse, gray surface:
$ R_"surf" = frac(1 - epsilon, epsilon A) $

Geometrical Resistance for a diffuse, gray surface:
$ R_"geo" = frac(1, A_1 F_(1-2)) $

Total Heat Flux between two diffuse, gray surfaces:
$ Q_"net"_(1-2) = frac(sigma_"SB" (T_1^4 - T_2^4), R_"surf"_1 + R_"geo"_(1-2) + R_"surf"_2) $

= Heat Pipes

$ L_c = V / A_s $

= General Equations of Usefulness
=== Hyperbolic Functions
$ sinh x = frac(e^x - e^(-x), 2) $
$ cosh x = frac(e^x + e^(-x), 2) $
$ tanh x = frac(e^x - e^(-x), e^x + e^(-x)) $

=== Surface Area of a Sphere
$ A_s = 4 pi r^2 $

=== Volume of a Sphere
$ V = 4/3 pi r^3 $

=== Surface Area of a Cylinder
$ A_s = 2 pi r h $

=== Volume of a Cylinder
$ V = pi r^2 h $

=== Error Function
$ erf(eta) = (2/pi^(1/2)) integral_0^eta e^(-u^2) d u $

Complimentary Error Function:
$ erfc(eta) = 1 - erf(eta) $