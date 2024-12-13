#set page(height: auto)

#let Bi = "Bi"
#let Fo = "Fo"
#let erf = "erf"
#let erfc = "erfc"
#let Re = "Re"
#let Nu = "Nu"
#let Pr = "Pr"
#let dt = "dt"
#let dx = "dx"
#let dy = "dy"

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

= Constants
=== Boltzmann Constant
$ k_B = 1.38 times 10^(-23) J/K $

=== Stefan-Boltzmann Constant
$ sigma_"SB" = 5.67 times 10^(-8) frac(W, m^2 K^4) $

= Governing Equation

$ E_"in" - E_"out" + E_"gen" = E_"st" $

#pagebreak()

/*

CONDUCTION

*/
= Conduction
=== Fourier Law of Heat Conduction
$ q^'' = -k frac(partial T, partial x) $

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

Total Surface Efficiency:
$ eta_t = 1 - A_f/A (1 - eta_f) $

Resistance of a Fin:
The resistance for a single fin
$ R_(t,f) = frac(1, h A_f eta_f) $

Thermal Resistance of Finned Surface:
The resistance for an entire finned surface
$ R = frac(1, h A eta_t) $

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

#pagebreak()

/*

CONVECTION

*/

= Convection
$ q_"local" = h A(T_s - T_infinity) $
$ q_"global" = U A (T_s - T_infinity) $

=== Thermal Resistance
$ R = frac(1, h A) $

=== Reynold's Number
$ Re_x = frac(rho u_infinity x, mu) = frac(u_infinity x, nu) $

=== Nusselt Number
$ Nu equiv frac(h L, k_f) $

Flat plate in Laminar Flow:
$ Nu = 0.664 Re_L^(1/2) Pr^(1/3) $

=== Prandtl Number
$ Pr equiv nu / k $

=== Convective Heat Transfer Coefficient
$ h = -frac(k_f, T_s - T_infinity) frac(partial T, partial y) |_(y=0) $

#pagebreak()

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
Saturated liquids evaporate by absorbing heat from a higher temperature and saturated vapors condense by releasing heat to a lower temperature.

$ Delta P_"capillary" - Delta P_"gravitational" = Delta P_L + Delta P_V  $

=== Effective Length of a Heat Pipe
$ L_"eff" = L_A + 1/2 (L_E + L_C) $

=== Capillary Pressure
Pore Radius ($r_p$)
Surface Tension ($gamma$)

$ Delta P_"capillary" = 2 gamma [ cos(theta_E)/ r_p - cos(theta_C) / r_p ] $

=== Gravitational Pressures
Angle of Heat Pipe ($phi.alt$)

$ Delta P_"gravitational" = rho g (L_"eff" sin(phi.alt)) $

=== Darcy Relation
Flow Rate ($Q_f$)
Permeability ($Kappa_p$)
Effective Wick Area ($A_w$)

$ Q_f = - Kappa_p / mu A_w [frac(Delta P_L, L_"eff")] $

=== Liquid Phase Change Pressure Drop
$ Delta P_L = frac(mu L_"eff", Kappa_p A_w) (dot(m), rho_L) $

=== Vapor Related Pressure Drop
$ Delta P_V = (1/2 rho_v v^2) (64 / Re) [frac(L_"eff", 4 A_w / rho)] $

=== Heat Pipe Equation
$ q dot L_"eff" = (frac(rho_L gamma h_"fg", mu_L)) [Kappa_p A_w] [2 / r_p] $

=== Weber Number
$ "We" = frac(rho v^2, gamma) l $

=== Issues with Heat Pipes
+ Sonic limit (liquid metal heat pipes)
+ Entrainment -> dry out
+ Boiling limitation (bubbles in wick) 
+ Chokig of vapor flux

#pagebreak()

/*

HEAT EXCHANGERS

*/
= Heat Exchangers
$ CC_h = dot(m)_h C_h $
$ CC_c = dot(m)_c C_c $

=== Power Lost
$ d q = - dot(m)_h C_h d T_h $

=== Power Gained
$ d q = dot(m)_c C_c d T_c $

=== Log-Mean Temperature Difference
$ Delta T_"LMTD" = frac( Delta T_2 - Delta T_1, ln(frac(Delta T_2, Delta T_1)) ) $

=== Heat Flux
$ ln(frac(Delta T_2, Delta T_1)) = - U A [ frac(1, dot(m)_c C_c) + frac(1, dot(m)_h C_h) ] $

$ q = U A Delta T_"LMTD" $

=== Max Heat Flux
$ q_"max" = CC_"min" (T_"hi" - T_"ci") $

=== Heat Exchanger Effectiveness
$ epsilon = q / q_"max" $

=== Number of Transfer Units
$ "NTU" = frac(U A, CC_min) $

== Parallel Flow
$ Delta T_1 = T_"hi" - T_"ci" $
$ Delta T_2 = T_"ho" - T_"co" $


== Counter Flow
$ Delta T_1 = T_"hi" - T_"co" $
$ Delta T_2 = T_"ho" - T_"ci" $

#pagebreak()

/*

THERMAL ELECTRICS

*/
= Thermal Electrics
$ m = R/RR $
$ Z = frac(S^2, K RR) $

=== Seebeck Effect
$ S = frac(Delta V, Delta T) $

== Heat Engine

=== Carnot Efficiency
$ eta_"carnot" = frac(T_H - T_C, T_H) $

=== Carzon-Ahlborn Efficiency
$ eta_(c - a) = 1 - sqrt(T_C / T_H) $

=== Thermal Conductance
$ K = frac(k_1 A_1, l_1) + frac(k_2 A_2, l_2) $

=== Resistance of Materials
$ RR = frac(rho_1 l_1, A_1) + frac(rho_2 l_2, A_2) $

=== Current
$ I = frac(S(T_H - T_C), R + RR) $

=== Thermal Electric Efficiency
$ eta_"TE" = eta_"carnot" [frac( (m / (m+1) ), 1 + (K RR)/S^2 ((m+1) / T_H) - 1/2 eta_"carnot" (1 / (m+1)) )] $

=== Geometric Constraint
$ sqrt(frac(k_2 rho_1, k_1 rho_2)) = A_1 / A_2 $

=== Optimal m
$ m_"optimal" = sqrt(1 + 1/2 Z (T_H + T_C)) $

== Peltier Cooler
$ q_Pi = q_o + q_T = Pi I $
$ q_T = 1/2 I^2 RR + K (T_H - T_C) $

=== Coefficient of Perfomance
$ "COP" = frac(T_C , T_H - T_C) $

=== Peltier Coefficient
$ Pi = S T_C $

=== Critical Current
$ I_C = frac(S T_C, RR) $

#pagebreak()

/*

GENERAL EQUATIONS

*/
= General Equations of Usefulness

$ L_c = V / A_s $

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

=== Complimentary Error Function:
$ erfc(eta) = 1 - erf(eta) $

=== Mass Flow Rate
$ dot(m) = Q_f dot rho_f $

#pagebreak()

/*

Differential Equations

*/
= Differential Equations
=== First Order ODE
$ frac(partial f, partial t) + p(t) f(t) = g(t) $

1. Find $mu(t)$
$ mu(t) = exp(integral p(t) d t) $

2. Multiply by $mu(t)$
$ frac(d, d t) (mu(t) f(t)) = mu(t) g(t) $

3. Integrate Both Sides
$ integral frac(d, dt) (mu(t) f(t)) dt = integral mu(t) g(t) dt $
$ mu(t) f(t) = integral mu(t) g(t) dt $
$ f(t) = 1 / mu(t) integral mu(t) g(t) dt $

=== Second Order ODE
$ frac(partial^2 y(t), partial t^2) + p(t) frac(partial y(t), partial t) + q(t) y(t) = g(t) $

1. Find the two roots (assume $y = exp(r t)$)
$ frac(partial y(t), partial t) = r exp(r t) $
$ frac(partial^2 y(t), partial t^2) = r^2 exp(r t) $

$ r_1, r_2 = frac( -p(t) plus.minus sqrt( p(t)^2 - 4 dot q(t) ), 2 ) $

2. Subsitute roots into general solution
$ y(t) = c_1 exp(r_1 t) + c_2 exp(r_2 t) $

3. Use Initial Conditions to find $c_1$ and $c_2$