import math
import random
import textwrap

import numpy as np
import streamlit as st

# ------------------ Page config ------------------ #
st.set_page_config(
    page_title="Fluids Companion for Dushto",
    layout="wide",
    initial_sidebar_state="expanded",
)

WRAP_WIDTH = 95


def wrap(text: str) -> str:
    return textwrap.fill(text.strip(), WRAP_WIDTH)


# ------------------ Personalized intro ------------------ #

def welcome_page():
    st.header("🌊 Fluid Mechanics – Personal space for Juairiya")
    st.write(
        wrap(
            """
This app was made by meeee just for you, so you can walk through fluid
mechanics calmly and ace your exam.

Use it like a tutor:
- Start with **Concept Tutor** when a topic feels difficult.
- Use **Calculators** to check your numbers on problems.
- Practice with **Guided Problems & Quizzes** until things become easier.
- skim **Exam Mode** to see if anything still feels off.

Take your time. You Got This Dushtoooo!!!
"""
        )
    )

    st.subheader("How to study with this tool (suggested plan)")
    st.markdown(
        """
1. **Pick today’s lecture topic** in the sidebar.  
2. In **Concept Tutor**, read the *Intuition* first, then the math.  
3. Do at least **one Guided Practice** problem for that topic.  
4. Use **Calculators** to verify any head losses / Reynolds / hydrostatics.  
5. End the session with **3–5 quiz questions** to lock it in.

Small consistent steps will make it amazing. You are absolutely going to
cracking this course.
"""
    )


# ------------------ Deep Concept Map ------------------ #

TOPICS = {
    "1. Fluid Basics & Pressure Field": {
        "1.1 Fluid vs Solid & Basic Properties": {
            "Intuition": """
A solid resists shear: you push sideways, it deforms a bit then stops.
A fluid never "settles" under shear: even a tiny shear stress causes
continuous deformation (flow). That's the essence of a fluid.

Key properties:
- Density ρ [kg/m³]: mass per unit volume.
- Specific weight γ = ρ g [N/m³]: weight per unit volume.
- Dynamic viscosity μ [Pa·s]: resistance to sliding of fluid layers.
- Kinematic viscosity ν = μ / ρ [m²/s]: "viscosity per density".
- Vapor pressure p_vap: pressure where liquid ↔ vapor in equilibrium:
  if local pressure drops below p_vap, cavitation (local boiling) can
  occur, damaging pumps/propellers.
""",
            "Math": r"""
Units & typical values at ~20°C:
- Water: ρ ≈ 1000 kg/m³, μ ≈ 1.0×10⁻³ Pa·s, ν ≈ 1.0×10⁻⁶ m²/s.
- Air:   ρ ≈ 1.2 kg/m³,    μ ≈ 1.8×10⁻⁵ Pa·s, ν ≈ 1.5×10⁻⁵ m²/s.

Dynamic viscosity (Newtonian fluid, simple shear):

    τ = μ (du/dy)

where
- τ is shear stress [Pa],
- u is velocity parallel to the surface [m/s],
- y is coordinate normal to surface [m].

Kinematic viscosity:

    ν = μ / ρ   [m²/s]

Specific weight:

    γ = ρ g   [N/m³],   g ≈ 9.81 m/s²
""",
            "Rigor": r"""
The Newtonian fluid model assumes a linear relationship between shear
stress and velocity gradient:

    τ = μ (du/dy)

This is a constitutive law: it connects the kinematics (velocity field)
to the dynamics (stress). Non-Newtonian fluids are those where the
relationship is not linear or depends on shear rate history. For most
undergrad mechanical engineering fluids problems, we assume the fluid
is Newtonian (e.g., water, air, simple oils) and incompressible unless
stated otherwise.
""",
        },
        "1.2 Newtonian vs Non-Newtonian": {
            "Intuition": """
Newtonian fluids: viscosity μ is constant for a given temperature.
If you double the shear rate, shear stress doubles.

Non-Newtonian fluids: viscosity effectively changes with shear rate
(shear-thinning ketchup, shear-thickening cornstarch-water mixtures,
Bingham plastics like toothpaste needing a yield stress).
""",
            "Math": r"""
Newtonian:

    τ = μ (du/dy).

Non-Newtonian (examples, just conceptually):

- Power-law fluid:   τ = K (du/dy)ⁿ
  where n < 1 → shear-thinning, n > 1 → shear-thickening.
- Bingham plastic:   τ = τ_y + μ_p (du/dy) for |τ| > τ_y,
                     no flow if |τ| < τ_y.
""",
            "Rigor": """
In this course, almost all derivations and standard formulas (pipe flow,
boundary layers, etc.) assume a Newtonian fluid, so τ is linearly
proportional to the velocity gradient. For exam purposes, you mainly
need to recognize "Newtonian" means τ ∝ du/dy and non-Newtonian means
that simple proportionality is not valid.
""",
        },
        "1.3 Pressure Variation in a Static Fluid": {
            "Intuition": """
In a static fluid there is no shear, only normal stress (pressure).
Pressure increases with depth because the fluid below must support the
weight of all the fluid above it.

At the same depth in a connected static fluid, pressure is the same
(Pascal's law), regardless of container shape.
""",
            "Math": r"""
Take z as upward. Force balance on a small fluid element:

    dp/dz = -ρ g

For constant density ρ:

    p(z) = p_ref + ρ g (z_ref - z)

Gauge pressure (relative to atmosphere):

    p_g = p - p_atm

Hydrostatic pressure difference between two depths z₁ and z₂:

    p₂ - p₁ = ρ g (z₁ - z₂)
""",
            "Rigor": """
Consider a differential fluid element (a small "box" of fluid). In
vertical equilibrium, the sum of vertical forces is zero:

    p(z+dz) A - p(z) A - ρ g A dz = 0

Divide by A dz and take the limit dz → 0:

    dp/dz = -ρ g

Integrating gives the hydrostatic distribution. This derivation assumes
no motion (static fluid) and constant density (incompressible). More
complex cases (e.g. compressible gas) integrate with ρ = ρ(p).
""",
        },
        "1.4 Manometers": {
            "Intuition": """
Manometers use columns of fluid to measure pressure differences. When
fluids are static and connected, pressure at points on the same
horizontal level in a continuous fluid is equal.

To solve manometer problems, you 'walk' through the fluid, adding ρ g Δz
when going down and subtracting when going up.
""",
            "Math": r"""
Basic idea for a U-tube manometer with a single fluid:

    p_A + ρ g h_A = p_B + ρ g h_B

If two fluids (ρ₁, ρ₂) are present:

    p_A + ρ₁ g h₁ + ρ₂ g h₂ = p_B

where h₁, h₂ are vertical height differences along different segments.

General method:
1) Choose a reference point.
2) Write pressure at that point from each side by moving through the
   fluid columns, adding/subtracting ρ g Δz.
3) Set expressions equal and solve for unknown pressure.
""",
            "Rigor": """
The manometer relation is just repeated application of hydrostatics:
in each static fluid region, dp/dz = -ρ g. Over a vertical distance
Δz, p changes by ρ g Δz. At interfaces, pressure is continuous if no
surface tension or membranes are involved. Matching pressures at a
common level allows you to relate pressures at A and B algebraically.
""",
        },
    },

    "2. Hydrostatics: Forces & Buoyancy": {
        "2.1 Forces on Plane Surfaces": {
            "Intuition": """
Because pressure increases linearly with depth, a submerged surface
feels a force distribution that increases with depth. The resultant
force is the integral of pressure over area, and it acts at a point
deeper than the centroid, called the center of pressure.
""",
            "Math": r"""
For a plane surface of area A in a liquid of density ρ with its
centroid at depth h_cg below the free surface:

    F_R = ρ g h_cg A

is the resultant hydrostatic force.

Depth of center of pressure h_cp:

    h_cp = I_g / (h_cg A) + h_cg

where I_g is the second moment of area of the surface about a horizontal
axis through the free surface, parallel to the surface.

For a vertical rectangular surface of width b and height H whose top
edge is at depth h_top:

    h_cg = h_top + H/2

    I_g (about free surface) = (b/3) [ (h_top + H)³ - h_top³ ]
""",
            "Rigor": r"""
Take a differential area dA at depth h. Pressure there is p = ρ g h.
Differential force:

    dF = p dA = ρ g h dA

Resultant force:

    F_R = ∬_A ρ g h dA = ρ g h_cg A

To find line of action, equate moments about the free surface:

    M = ∬ h (ρ g h dA) = ρ g ∬ h² dA

and define:

    F_R h_cp = M  ⇒  h_cp = ( ∬ h² dA ) / ( ∬ h dA )

Using the parallel-axis theorem, this gives the standard formula

    h_cp = I_g / (h_cg A) + h_cg
""",
        },
        "2.2 Buoyancy & Archimedes": {
            "Intuition": """
Buoyancy is the net vertical force on a submerged or floating body due
to the hydrostatic pressure distribution. Archimedes' principle says:

The buoyant force equals the weight of the displaced fluid.

For a floating object, weight of the object = buoyant force.
""",
            "Math": r"""
Buoyant force:

    F_B = ρ_fluid g V_disp

where V_disp is the volume of displaced fluid.

For a fully submerged body:

    V_disp = V_body

For a floating body:

    W_body = F_B  ⇒  m_body g = ρ_fluid g V_disp
                 ⇒  V_disp = m_body / ρ_fluid
""",
            "Rigor": """
Consider the pressure forces on the surface of the body. The vertical
components add up to a net upward force equal to the weight of the
displaced fluid. This is seen by comparing with the weight of a
"phantom" fluid occupying the same volume as the body. This rigorous
integral result simplifies beautifully to Archimedes' principle.
""",
        },
        "2.3 Stability of Floating Bodies (Metacentric)": {
            "Intuition": """
A floating body can be stable or unstable when tilted. Stability depends
on the relative positions of the center of gravity G, center of buoyancy
B, and the metacenter M.

- If M is above G (GM > 0), small tilts produce a righting moment
  → stable.
- If M is below G (GM < 0), small tilts lead to overturning
  → unstable.
""",
            "Math": r"""
For small angles:

    GM = BM - BG

where:
- G is center of gravity,
- B is center of buoyancy (centroid of displaced volume),
- M is metacenter, found from waterline geometry.

For a rectangular floating prism (ship-style approximation):

    BM = I_waterline / V_disp

where I_waterline is the second moment of area of the waterplane
about the tilt axis, V_disp is displaced volume.

Stability criterion:

    GM > 0 ⇒ stable
    GM < 0 ⇒ unstable
""",
            "Rigor": """
When the body is tilted slightly, the shape of the displaced volume
changes. The new center of buoyancy B' shifts relative to B. The
metacenter M is defined as the intersection of the vertical line
through B' with the original vertical line through B for small angles.

The restoring (or overturning) moment is:

    M_restoring = W (GM) sin(θ)

For small θ, sign of GM determines whether the moment restores the
body to equilibrium (stable) or not (unstable).
""",
        },
    },

    "3. Flow Kinematics & Differential View": {
        "3.1 Eulerian Description & Streamlines": {
            "Intuition": """
Instead of following individual fluid particles (Lagrangian), we use
the Eulerian description: we define fields u(x,y,z,t), v(x,y,z,t),
w(x,y,z,t) and p(x,y,z,t). At each point in space and time, the flow
has a velocity and pressure.

Streamlines: lines everywhere tangent to the instantaneous velocity
field. In steady flow, streamlines are also pathlines (particle paths).
""",
            "Math": r"""
Velocity field in 3D:

    \mathbf{u}(x,y,z,t) = u \hat{i} + v \hat{j} + w \hat{k}.

Streamlines at an instant t satisfy:

    \frac{dx}{u} = \frac{dy}{v} = \frac{dz}{w}.

Steady flow:

    ∂/∂t = 0   for all flow properties at a fixed point.
""",
            "Rigor": """
Streamlines are defined such that the tangent vector at each point is
aligned with the velocity vector. In mathematical terms, for a curve
x(s), y(s), z(s) parameterized by s, the tangent satisfies:

    (dx/ds, dy/ds, dz/ds) ∥ (u, v, w).

Hence the differential relation:

    dx/u = dy/v = dz/w.

This is purely geometric and does not require any dynamical equations.
""",
        },
        "3.2 Material Derivative & Acceleration": {
            "Intuition": """
Acceleration of a fluid particle is not just the local time derivative
of velocity at a point. A particle moves through space, so even if the
flow is steady at each fixed point, the particle may speed up or slow
down as it enters regions with different velocity.

This is captured by the material derivative.
""",
            "Math": r"""
Velocity field u(x,y,z,t). Material derivative of any scalar φ:

    Dφ/Dt = ∂φ/∂t + u ∂φ/∂x + v ∂φ/∂y + w ∂φ/∂z
           = ∂φ/∂t + (\mathbf{u} ⋅ ∇)φ.

Acceleration:

    \mathbf{a} = D\mathbf{u}/Dt
               = ∂\mathbf{u}/∂t + (\mathbf{u}⋅∇)\mathbf{u}.

Steady flow (∂/∂t = 0) still can have nonzero acceleration via the
convective term (\mathbf{u}⋅∇)\mathbf{u}.
""",
            "Rigor": """
In a Lagrangian view, a fluid particle's position X(t) evolves, and its
velocity is u(X(t), t). Acceleration is d/dt[u(X(t), t)].

By chain rule:

    d u / dt = ∂u/∂t + (dx/dt) ∂u/∂x + (dy/dt) ∂u/∂y + (dz/dt) ∂u/∂z
             = ∂u/∂t + u ∂u/∂x + v ∂u/∂y + w ∂u/∂z.

Extending to vector form gives the full material derivative. This is a
key link between Eulerian and Lagrangian descriptions.
""",
        },
    },

    "4. Control Volume: Continuity, Energy, Momentum": {
        "4.1 Continuity Equation (Integral)": {
            "Intuition": """
Continuity is just conservation of mass. For a control volume in steady,
incompressible flow, the mass flow rate that goes in must equal the mass
flow rate that comes out.

For a single inlet and outlet with constant density:

    V₁ A₁ = V₂ A₂.
""",
            "Math": r"""
General integral continuity (any fluid):

    d/dt ∫_CV ρ dV + ∑_out ∫_A ρ V_n dA - ∑_in ∫_A ρ V_n dA = 0.

Steady, incompressible, uniform velocity at each section:

    ∑_in (ρ V_n A) = ∑_out (ρ V_n A)

or

    ∑_in (V_n A) = ∑_out (V_n A).
""",
            "Rigor": """
Starting from conservation of mass for a system and applying the Reynolds
Transport Theorem gives the continuity equation for a control volume.
The terms represent the rate of mass accumulation in the CV plus net
outflow across its boundaries. Setting this sum to zero enforces that
mass is neither created nor destroyed.
""",
        },
        "4.2 Energy Equation & Bernoulli": {
            "Intuition": """
The energy equation for a control volume balances energy per unit weight
between inlet(s) and outlet(s), plus what pumps, turbines, and losses
do. Bernoulli is just the simplified, ideal special case along a
streamline with no pumps, no turbines, and no losses.
""",
            "Math": r"""
Steady incompressible flow between 1 and 2 (single streamline):

    p₁/γ + V₁²/(2g) + z₁ + h_p
  = p₂/γ + V₂²/(2g) + z₂ + h_t + h_L

where:
- p/γ is pressure head,
- V²/(2g) is velocity head,
- z is elevation head,
- h_p is pump head added,
- h_t is turbine head removed,
- h_L is total head loss between 1 and 2.

Ideal Bernoulli (no pump, no turbine, no losses, along a streamline):

    p₁/γ + V₁²/(2g) + z₁ = p₂/γ + V₂²/(2g) + z₂.
""",
            "Rigor": """
Derive from the steady-flow energy equation for a control volume by
choosing CV boundaries through inlets and outlets and including shaft
work terms (pumps/turbines). Assume incompressible fluid, negligible
heat transfer, and negligible internal energy changes to reduce to the
head form. Along a streamline in an inviscid flow (no losses and no
shaft work), this further reduces to Bernoulli's equation.
""",
        },
        "4.3 Momentum Equation (Control Volume)": {
            "Intuition": """
The linear momentum equation is Newton's 2nd law for a control volume:
net external force on the fluid equals net rate of change of momentum
out minus in. It is the tool for forces on pipe bends, nozzles, and
jets hitting plates.
""",
            "Math": r"""
General integral momentum in direction s:

    ∑ F_{s,ext} =
        ∑_out (ρ Q V_s)_out - ∑_in (ρ Q V_s)_in

where
- F_{s,ext} are external forces (pressure, weight, reaction),
- Q is volumetric flow rate,
- V_s is the component of velocity in direction s.

Example (horizontal bend, incompressible, single inlet/outlet):

    F_x =
        ρ Q V_{2,x} - ρ Q V_{1,x}

Reaction on bend = -F_x (equal and opposite).
""",
            "Rigor": """
Start from conservation of momentum for a system and apply the Reynolds
Transport Theorem with extensive property B = linear momentum. The
result relates external forces to the net outflow of momentum across
the control surface. Decomposing into component directions yields the
scalar momentum equations often used in pipe bend and jet problems.
""",
        },
    },

    "5. Dimensional Analysis & Similarity": {
        "5.1 Buckingham π & Dimensionless Groups": {
            "Intuition": """
Dimensional analysis reduces the number of variables by grouping them
into dimensionless combinations (π groups). These reveal which
dimensionless numbers control the physics (Re, Fr, Ma, etc.) and help
design experiments and scale models.
""",
            "Math": r"""
Buckingham π theorem:
If a physical relation involves n variables and k fundamental dimensions
(M, L, T, ...), you can form n - k independent dimensionless groups
(π₁, π₂, ..., π_{n-k}).

Example: drag on a sphere:
Variables: F_D, ρ, μ, V, D
Dimensions:
- F_D: [M L T⁻²]
- ρ:   [M L⁻³]
- μ:   [M L⁻¹ T⁻¹]
- V:   [L T⁻¹]
- D:   [L]

n = 5, k = 3 ⇒ 2 π groups.

Common choice:
- π₁ = C_D = F_D / (½ ρ V² A)
- π₂ = Re = ρ V D / μ

So we can write:

    C_D = f(Re).
""",
            "Rigor": """
Form a dimension matrix with each row corresponding to a fundamental
dimension and each column to a variable. The null space of this matrix
gives exponent vectors for dimensionless products of variables. The
Buckingham π theorem guarantees that all physically meaningful
relationships can be expressed in terms of a complete set of π groups.
""",
        },
        "5.2 Similarity (Geometric, Kinematic, Dynamic)": {
            "Intuition": """
Similarity between a model and prototype requires:
- Geometric similarity: same shapes, scaled lengths.
- Kinematic similarity: similar velocity patterns (correct dimensionless
  numbers like Re, Fr, Ma).
- Dynamic similarity: force ratios match (again via dimensionless
  numbers).

You can't match all dimensionless numbers simultaneously for complex
flows, so you choose the ones that matter most (e.g., Re for pipe
flows, Fr for free surface flows).
""",
            "Math": r"""
- Geometric similarity:
  length ratios same in all directions, shapes are similar.

- Kinematic similarity:
  e.g., matching Reynolds number:
      Re_model = Re_prototype
      (ρ V L / μ)_model = (ρ V L / μ)_prototype

- Dynamic similarity:
  key dimensionless groups (Re, Fr, Ma, etc.) are the same for model
  and prototype.

For free-surface flows (like ship models), Froude number often dominates:

    Fr = V / √(g L)
""",
            "Rigor": """
The equations of motion for fluids (Navier–Stokes, etc.) can be written
in dimensionless form using characteristic scales. The resulting
dimensionless groups multiply terms in the equations. Dynamic similarity
between model and prototype requires matching these groups, ensuring
the dimensionless equations are identical. In practice, compromises are
made by prioritizing the most influential groups for the problem at hand.
""",
        },
    },

    "6. Internal Flow: Laminar, Turbulent, Losses": {
        "6.1 Reynolds Number & Regimes": {
            "Intuition": """
Reynolds number compares inertial forces to viscous forces. High Re:
inertia dominates and flow tends to be turbulent. Low Re: viscosity
dominates and flow stays laminar.

In a circular pipe, the thresholds are roughly:
- Re < 2000: laminar
- 2000 < Re < 4000: transitional
- Re > 4000: turbulent
""",
            "Math": r"""
Reynolds number in a pipe:

    Re = ρ V D / μ,

where
- ρ is density [kg/m³],
- V is average velocity [m/s],
- D is pipe diameter [m],
- μ is dynamic viscosity [Pa·s].
""",
            "Rigor": """
The Reynolds number arises naturally when the Navier–Stokes equations
are non-dimensionalized with characteristic velocity V and length D.
It multiplies the inertial term relative to the viscous term. The
critical Reynolds number for transition is determined experimentally,
not from first principles, and depends on disturbances and geometry.
""",
        },
        "6.2 Fully Developed Laminar Pipe Flow": {
            "Intuition": """
In fully developed laminar flow, the velocity profile in a circular
pipe is parabolic: zero at the wall (no-slip) and maximum at the
center. Pressure drops linearly along the pipe, and viscous effects
control the motion.
""",
            "Math": r"""
Navier–Stokes solution for steady, fully developed, laminar, incompressible
flow in a circular pipe of radius R and length L:

Velocity profile:

    u(r) = (Δp / (4 μ L)) (R² - r²),

where r is radial coordinate from center.

Max velocity at r=0:

    u_max = (Δp / (4 μ L)) R².

Average velocity:

    V = (1/A) ∫_A u dA = u_max / 2.

Volumetric flow rate:

    Q = (π R⁴ / (8 μ L)) Δp.

Darcy friction factor for laminar flow:

    f = 64 / Re.
""",
            "Rigor": """
Simplify Navier–Stokes: axisymmetric, steady, fully developed (no axial
variation in velocity), no body forces in axial direction. The only
non-zero velocity component is u_z(r). This reduces N-S to an ODE in r,
which integrates to give the parabolic profile. Wall shear stress
relates to pressure gradient, leading to Hagen–Poiseuille law and
f = 64/Re in the Darcy–Weisbach formulation.
""",
        },
        "6.3 Turbulent Flow & Moody Chart": {
            "Intuition": """
Turbulent flow is chaotic, with strong mixing and a much flatter
velocity profile near the center of the pipe. Friction factors are
higher than laminar, and depend on both Reynolds number and relative
roughness ε/D.

The Moody chart summarizes decades of experiments: it's the main tool
for turbulent friction factors.
""",
            "Math": r"""
Darcy–Weisbach equation for head loss:

    h_L = f (L/D) (V² / (2g)).

Laminar:   f = 64 / Re.
Turbulent: f = f(Re, ε/D) from Moody chart or correlations.

Swamee–Jain explicit correlation for turbulent flow:

    f = 0.25
        / [ log₁₀( ε/(3.7D) + 5.74/Re⁰·⁹ ) ]²

where ε is absolute roughness, D is diameter.
""",
            "Rigor": """
In turbulent flow, the time-averaged Navier–Stokes equations involve
Reynolds stresses, which require turbulence modeling for closure.
Analytical solutions are not available in general, so friction factors
are obtained by experiments and empirical correlations. The Moody chart
is a graphical representation of these empirical relationships.
""",
        },
        "6.4 Major & Minor Losses, Pipe Networks": {
            "Intuition": """
Energy (head) is lost due to friction along straight pipe lengths
(major losses) and at fittings, valves, entrances, exits, and changes
in area (minor losses). Pipe networks combine these to determine how
flow splits and what total head is required from pumps.
""",
            "Math": r"""
Major loss (Darcy–Weisbach):

    h_{L,major} = f (L/D) (V²/(2g)).

Minor losses:

    h_{L,minor} = K (V²/(2g)),

where K is a loss coefficient depending on the type of fitting.

Total head loss:

    h_L,total = ∑ h_{L,major} + ∑ h_{L,minor}.

Series pipes (same Q): head losses add.
Parallel pipes (same ΔH): flow splits such that head loss in each
branch is equal; use continuity at junctions plus equal-loss conditions.
""",
            "Rigor": """
Minor loss coefficients are obtained from experiments or detailed
turbulence modeling. In network analysis, each loop or path is treated
with energy equations and continuity at junctions. The Hardy–Cross
method iteratively adjusts flow rates in loops to satisfy conservation
of mass and energy simultaneously.
""",
        },
    },

    "7. External Flow & Boundary Layers (Intro)": {
        "7.1 Boundary Layer Concept": {
            "Intuition": """
In high-Reynolds-number external flows, viscous effects are confined to
a thin region near the wall: the boundary layer. Outside it, the flow
is nearly inviscid. The boundary layer starts at the leading edge, grows
downstream, and may transition from laminar to turbulent. Separation
of the boundary layer can cause huge drag.
""",
            "Math": r"""
Order-of-magnitude estimates for laminar boundary layer on a flat plate:

    Re_x = ρ U∞ x / μ

BL thickness approximately:

    δ(x) ~ 5 x / √(Re_x).

Transition to turbulence near Re_x ≈ 5×10⁵ for a smooth flat plate
(rough value; depends on disturbance level).
""",
            "Rigor": """
By non-dimensionalizing the Navier–Stokes equations for flow over a
flat plate at high Re and focusing on the thin region near the wall,
terms of smaller magnitude are neglected, leading to the boundary-layer
equations (Prandtl). These are simpler than full N-S and can be solved
approximately to obtain velocity profiles and wall shear stress.
""",
        },
        "7.2 Drag on Bluff Bodies & Drag Coefficient": {
            "Intuition": """
Objects in a fluid experience drag due to pressure differences and
viscous shear. We package all the messy details into a dimensionless
drag coefficient C_D, which is obtained experimentally.

Drag force formula:

    F_D = ½ ρ V² C_D A,

where A is reference area (e.g., projected frontal area).
""",
            "Math": r"""
Drag coefficient definition:

    C_D = F_D / (½ ρ V² A).

For a sphere, cylinder, airfoil, etc., C_D depends on Reynolds number
and, for some bodies, angle of attack and surface roughness. Data is
usually presented as curves C_D(Re) from experiments.
""",
            "Rigor": """
Drag arises from both pressure (form drag) and viscous shear (skin
friction). Resolving and integrating the pressure and shear stress
distributions over the body surface gives F_D. In practice this is
rarely done from first principles in undergrad fluids; instead, we
use experimentally measured C_D and the drag formula.
""",
        },
    },

    "8. Open-Channel & Compressible (Intro Level)": {
        "8.1 Froude Number & Critical Flow": {
            "Intuition": """
In shallow-water flows (open-channel), the Froude number compares flow
speed to the speed of small surface waves. It plays a role similar to
Mach number in compressible flow.

- Fr < 1: subcritical (wave speed > flow speed) – disturbances can move upstream.
- Fr = 1: critical flow – min specific energy.
- Fr > 1: supercritical (flow outruns waves).
""",
            "Math": r"""
For a rectangular channel with depth y:

    Fr = V / √(g y).

Classification:
- Fr < 1: subcritical
- Fr = 1: critical
- Fr > 1: supercritical
""",
            "Rigor": """
Derive from specific energy:

    E = y + V²/(2g),

and the wave speed in shallow water a = √(g y). The condition for
minimum specific energy for a given discharge leads to Fr = 1 as the
critical condition. Momentum and energy analyses of hydraulic jumps
use Fr to relate conjugate depths.
""",
        },
        "8.2 Mach Number & Isentropic Nozzles (Very Intro)": {
            "Intuition": """
When gas speeds approach the speed of sound, compressibility matters.
Mach number M = V/a compares flow speed to sound speed. In an
isentropic nozzle with sufficiently low downstream pressure, flow can
accelerate to M=1 at the throat (choked flow), then to supersonic
downstream of the throat in a diverging section.
""",
            "Math": r"""
Mach number:

    M = V / a,   a = √(γ R T).

Isentropic relations (ideal gas):

    T0/T = 1 + (γ - 1)/2 M²
    p0/p = (T0/T)^{γ/(γ - 1)}

Area–Mach relation (no need to memorize exact form unless required):

    A/A* = (1/M) [ (2/(γ+1)) (1 + (γ-1)/2 M²) ]^{(γ+1)/(2(γ-1))}
""",
            "Rigor": """
Assuming steady, 1D, adiabatic, inviscid flow of an ideal gas with
no shaft work, combine continuity, momentum, and energy equations and
the ideal gas law. Imposing isentropic conditions (no entropy
generation) leads to the isentropic flow relations between temperature,
pressure, Mach number, and area ratio.
""",
        },
    },
}


# ------------------ Concept Tutor UI ------------------ #

def concept_tutor_page():
    st.header("🧠 Concept Tutor")

    col1, col2 = st.columns([1, 2])
    with col1:
        main_topics = list(TOPICS.keys())
        main_topic = st.selectbox("Main topic", main_topics)
        subtopics = list(TOPICS[main_topic].keys())
        subtopic = st.selectbox("Subtopic", subtopics)

    topic_data = TOPICS[main_topic][subtopic]

    with col2:
        st.subheader(subtopic)
        with st.expander("Intuition (read this first)", expanded=True):
            st.write(wrap(topic_data["Intuition"]))
        with st.expander("Math & key formulas", expanded=False):
            st.code(topic_data["Math"])
        with st.expander("Derivation / deeper idea", expanded=False):
            st.write(wrap(topic_data["Rigor"]))


# ------------------ Calculators ------------------ #

def reynolds_calculator():
    st.subheader("Reynolds Number & Flow Regime")

    col1, col2 = st.columns(2)
    with col1:
        rho = st.number_input("Density ρ [kg/m³]", value=1000.0)
        mu = st.number_input("Dynamic viscosity μ [Pa·s]", value=1.0e-3, format="%.3e")
        D = st.number_input("Pipe diameter D [m]", value=0.05)
        V = st.number_input("Average velocity V [m/s]", value=1.5)

    if st.button("Compute Re and regime", key="re_btn"):
        Re = rho * V * D / mu
        st.write(f"**Re = {Re:,.1f}**")

        if Re < 2000:
            regime = "Laminar (Re < 2000)"
            color = "🟢"
        elif Re < 4000:
            regime = "Transitional (2000 < Re < 4000)"
            color = "🟡"
        else:
            regime = "Turbulent (Re > 4000)"
            color = "🔴"
        st.write(f"{color} Flow regime: **{regime}**")


def pipe_head_loss_calculator():
    st.subheader("Pipe Head Loss (Darcy–Weisbach)")

    col1, col2 = st.columns(2)
    with col1:
        rho = st.number_input("Density ρ [kg/m³]", value=1000.0, key="rho_pipe")
        mu = st.number_input(
            "Dynamic viscosity μ [Pa·s]", value=1.0e-3, format="%.3e", key="mu_pipe"
        )
        D = st.number_input("Diameter D [m]", value=0.05, key="D_pipe")
        L = st.number_input("Length L [m]", value=10.0, key="L_pipe")
        V = st.number_input("Mean velocity V [m/s]", value=1.5, key="V_pipe")
        eps = st.number_input(
            "Roughness ε [m] (e.g. 1.5e-5 for commercial steel)",
            value=1.5e-5,
            format="%.1e",
            key="eps_pipe",
        )

    if st.button("Compute head loss", key="pipe_btn"):
        g = 9.81
        Re = rho * V * D / mu

        if Re < 2000:
            f = 64.0 / Re
            regime = "Laminar"
        else:
            f = 0.25 / (
                math.log10((eps / (3.7 * D)) + (5.74 / (Re**0.9))) ** 2
            )
            regime = "Turbulent (Swamee–Jain)"

        hL = f * (L / D) * (V**2 / (2 * g))

        st.write(f"Re = {Re:,.1f} → **{regime}**")
        st.write(f"Friction factor **f = {f:.5f}**")
        st.write(f"Head loss **hₗ = {hL:.4f} m**")


def hydrostatic_force_calculator():
    st.subheader("Hydrostatic Force on a Vertical Rectangular Gate")

    col1, col2 = st.columns(2)
    with col1:
        rho = st.number_input("Fluid density ρ [kg/m³]", value=1000.0, key="rho_hs")
        b = st.number_input("Width b [m]", value=2.0, key="b_hs")
        H = st.number_input("Height H [m]", value=3.0, key="H_hs")
        h_top = st.number_input(
            "Depth of top edge below free surface h_top [m]",
            value=1.0,
            key="htop_hs",
        )

    if st.button("Compute hydrostatic force", key="hs_btn"):
        g = 9.81
        A = b * H
        h_cg = h_top + H / 2
        F_R = rho * g * h_cg * A
        I_g = b / 3 * ((h_top + H) ** 3 - h_top**3)
        h_cp = I_g / (h_cg * A)

        st.write(f"Area A = {A:.2f} m²")
        st.write(f"Centroid depth h_cg = {h_cg:.2f} m")
        st.write(f"Resultant force **F_R = {F_R:,.1f} N**")
        st.write(f"Center of pressure depth **h_cp = {h_cp:.3f} m**")


def hydraulic_jump_calculator():
    st.subheader("Hydraulic Jump (Rectangular Channel)")

    col1, col2 = st.columns(2)
    with col1:
        y1 = st.number_input("Upstream depth y₁ [m]", value=0.3, key="y1_hj")
        V1 = st.number_input("Upstream velocity V₁ [m/s]", value=6.0, key="V1_hj")

    if st.button("Compute hydraulic jump", key="hj_btn"):
        g = 9.81
        Fr1 = V1 / math.sqrt(g * y1)
        y2 = 0.5 * y1 * (-1 + math.sqrt(1 + 8 * Fr1**2))

        st.write(f"Froude number upstream: **Fr₁ = {Fr1:.2f}**")
        if Fr1 <= 1:
            st.warning("Fr₁ ≤ 1 → no classical hydraulic jump (flow not supercritical).")
        else:
            st.write(f"Conjugate depth downstream: **y₂ ≈ {y2:.3f} m**")


def pump_sizing_calculator():
    st.subheader("Pump Sizing (Head & Power)")

    col1, col2 = st.columns(2)
    with col1:
        rho = st.number_input("Fluid density ρ [kg/m³]", value=1000.0, key="rho_pump")
        Q = st.number_input("Flow rate Q [m³/s]", value=0.02, key="Q_pump")
        delta_z = st.number_input("Elevation increase Δz [m]", value=10.0, key="dz_pump")
        delta_p = st.number_input(
            "Outlet pressure - inlet pressure Δp [Pa]", value=200000.0, key="dp_pump"
        )
        hL = st.number_input("Estimated total head loss hₗ [m]", value=5.0, key="hL_pump")
        eta = st.slider("Pump efficiency η", min_value=0.4, max_value=0.95, value=0.7)

    if st.button("Compute pump requirements", key="pump_btn"):
        g = 9.81
        gamma = rho * g
        H_static = delta_z + delta_p / gamma
        H_total = H_static + hL
        P_hydraulic = rho * g * Q * H_total
        P_shaft = P_hydraulic / eta

        st.write(f"Static head: **H_static = {H_static:.2f} m**")
        st.write(f"Total required head: **H_total = {H_total:.2f} m**")
        st.write(f"Hydraulic power: **{P_hydraulic/1000:.2f} kW**")
        st.write(f"Required shaft power: **{P_shaft/1000:.2f} kW**")


def calculators_page():
    st.header("🧮 Calculators")
    calc = st.selectbox(
        "Pick a calculator",
        [
            "Reynolds number & regime",
            "Pipe head loss",
            "Hydrostatic force on gate",
            "Hydraulic jump (rectangular channel)",
            "Pump sizing",
        ],
    )

    if calc == "Reynolds number & regime":
        reynolds_calculator()
    elif calc == "Pipe head loss":
        pipe_head_loss_calculator()
    elif calc == "Hydrostatic force on gate":
        hydrostatic_force_calculator()
    elif calc == "Hydraulic jump (rectangular channel)":
        hydraulic_jump_calculator()
    elif calc == "Pump sizing":
        pump_sizing_calculator()


# ------------------ Visuals Page ------------------ #

def visuals_page():
    st.header("📈 Visual Intuition")

    st.subheader("Laminar Pipe Flow – Velocity Profile")
    st.markdown(
        "Adjust radius and pressure gradient and see how the parabolic profile changes."
    )

    col1, col2 = st.columns(2)
    with col1:
        R = st.slider("Pipe radius R [m]", 0.02, 0.15, 0.05, step=0.01)
        dpdx = st.slider(
            "Pressure gradient -Δp/L [Pa/m] (positive number)", 100.0, 5000.0, 1000.0
        )

    mu = 1.0e-3
    # axial direction z; using Hagen–Poiseuille u(r) = (Δp / (4 μ L))(R² - r²)
    # here dpdx represents Δp/L (positive), so u(r) = (dpdx / (4 μ))(R² - r²)
    r_vals = np.linspace(0, R, 100)
    u_vals = (dpdx / (4 * mu)) * (R**2 - r_vals**2)

    with col2:
        st.line_chart(
            {"u(r) [m/s]": u_vals},
            x=None,
        )
        st.caption("Velocity is maximum at the center and zero at the wall (no-slip).")


# ------------------ Guided Practice Problems ------------------ #

def hydrostatic_practice():
    st.subheader("Hydrostatics Practice – Vertical Gate")

    if "hs_prob" not in st.session_state:
        st.session_state.hs_prob = None

    if st.button("New random problem", key="hs_new"):
        b = random.choice([1.5, 2.0, 2.5])
        H = random.choice([2.0, 3.0, 4.0])
        h_top = random.choice([0.5, 1.0, 1.5])
        rho = 1000.0
        g = 9.81
        st.session_state.hs_prob = (b, H, h_top, rho, g)

    prob = st.session_state.hs_prob
    if prob is None:
        st.info("Click 'New random problem' to generate a scenario.")
        return

    b, H, h_top, rho, g = prob
    st.write(
        wrap(
            f"A vertical rectangular gate of width b = {b} m and height H = {H} m "
            f"is submerged in water (ρ = 1000 kg/m³). The top of the gate is "
            f"h_top = {h_top} m below the free surface. "
            "Find (1) the resultant force on the gate and (2) the depth of the "
            "center of pressure below the free surface."
        )
    )

    col1, col2 = st.columns(2)
    with col1:
        F_user = st.number_input("Your F_R [N]", value=0.0, key="hs_F_user")
    with col2:
        hcp_user = st.number_input("Your h_cp [m]", value=0.0, key="hs_hcp_user")

    if st.button("Check answers", key="hs_check"):
        A = b * H
        h_cg = h_top + H / 2
        F_true = rho * g * h_cg * A
        I_g = b / 3 * ((h_top + H) ** 3 - h_top**3)
        h_cp_true = I_g / (h_cg * A)

        st.write(f"Correct F_R = {F_true:,.1f} N")
        st.write(f"Correct h_cp = {h_cp_true:.3f} m")

        def feedback(user, true, name):
            if true == 0:
                return
            err = abs(user - true) / abs(true)
            if err < 0.05:
                st.success(f"{name}: within 5% ✔")
            elif err < 0.15:
                st.warning(f"{name}: within 15% — close, re-check arithmetic.")
            else:
                st.error(f"{name}: large error — re-check formula/units.")

        feedback(F_user, F_true, "F_R")
        feedback(hcp_user, h_cp_true, "h_cp")

    with st.expander("Show full solution"):
        A = b * H
        h_cg = h_top + H / 2
        F_true = rho * g * h_cg * A
        I_g = b / 3 * ((h_top + H) ** 3 - h_top**3)
        h_cp_true = I_g / (h_cg * A)
        st.markdown(
            f"""
**1. Geometry**

- b = {b} m, H = {H} m  
- Area A = bH = {A:.2f} m²  
- Top depth h_top = {h_top} m  
- Centroid depth h_cg = h_top + H/2 = {h_cg:.2f} m  

**2. Resultant force**

F_R = ρ g h_cg A  
= {rho:.0f} × {g:.2f} × {h_cg:.2f} × {A:.2f}  
= {F_true:,.1f} N

**3. Center of pressure**

I_g (about free surface) = (b/3)[(h_top+H)³ − h_top³] = {I_g:.3f} m⁴

h_cp = I_g / (h_cg A) = {h_cp_true:.3f} m below free surface.
"""
        )


def pipe_flow_practice():
    st.subheader("Pipe Flow Practice – Head Loss & Regime")

    if "pipe_prob" not in st.session_state:
        st.session_state.pipe_prob = None

    if st.button("New random problem", key="pipe_new"):
        D = random.choice([0.03, 0.05, 0.08])
        L = random.choice([10.0, 20.0, 30.0])
        V = random.choice([1.0, 1.5, 2.5])
        eps = random.choice([1.5e-5, 4.5e-5])
        rho = 1000.0
        mu = 1e-3
        st.session_state.pipe_prob = (D, L, V, eps, rho, mu)

    prob = st.session_state.pipe_prob
    if prob is None:
        st.info("Click 'New random problem' to generate a scenario.")
        return

    D, L, V, eps, rho, mu = prob
    st.write(
        wrap(
            f"Water (ρ = 1000 kg/m³, μ = 1.0×10⁻³ Pa·s) flows in a circular pipe "
            f"of diameter D = {D} m and length L = {L} m with average velocity "
            f"V = {V} m/s. Pipe roughness ε = {eps:.1e} m. "
            "Find Re, friction factor f, and head loss h_L."
        )
    )

    col1, col2, col3 = st.columns(3)
    with col1:
        Re_user = st.number_input("Your Re", value=0.0, key="pf_Re_user")
    with col2:
        f_user = st.number_input("Your f", value=0.0, key="pf_f_user")
    with col3:
        hL_user = st.number_input("Your h_L [m]", value=0.0, key="pf_hL_user")

    if st.button("Check answers", key="pf_check"):
        g = 9.81
        Re = rho * V * D / mu
        if Re < 2000:
            f = 64.0 / Re
        else:
            f = 0.25 / (
                math.log10((eps / (3.7 * D)) + (5.74 / (Re**0.9))) ** 2
            )
        hL = f * (L / D) * (V**2 / (2 * g))

        st.write(f"Re = {Re:,.1f}")
        st.write(f"f = {f:.5f}")
        st.write(f"h_L = {hL:.4f} m")

        def fb(user, true, name):
            if true == 0:
                return
            err = abs(user - true) / abs(true)
            if err < 0.05:
                st.success(f"{name}: within 5% ✔")
            elif err < 0.15:
                st.warning(f"{name}: within 15% — decent.")
            else:
                st.error(f"{name}: re-check formula and units.")

        fb(Re_user, Re, "Re")
        fb(f_user, f, "f")
        fb(hL_user, hL, "h_L")


def hydraulic_jump_practice():
    st.subheader("Hydraulic Jump Practice – Conjugate Depth")

    if "hj_prob" not in st.session_state:
        st.session_state.hj_prob = None

    if st.button("New random problem", key="hj_new"):
        y1 = random.choice([0.2, 0.3, 0.4])
        V1 = random.choice([5.0, 6.0, 7.0])
        st.session_state.hj_prob = (y1, V1)

    prob = st.session_state.hj_prob
    if prob is None:
        st.info("Click 'New random problem' to generate a scenario.")
        return

    y1, V1 = prob
    g = 9.81
    st.write(
        wrap(
            f"In a horizontal rectangular channel, upstream depth y₁ = {y1} m "
            f"and velocity V₁ = {V1} m/s. Assume a hydraulic jump forms. "
            "Find the conjugate depth y₂ using the standard momentum relation."
        )
    )

    y2_user = st.number_input("Your y₂ [m]", value=0.0, key="hj_y2_user")

    if st.button("Check answer", key="hj_check"):
        Fr1 = V1 / math.sqrt(g * y1)
        y2_true = 0.5 * y1 * (-1 + math.sqrt(1 + 8 * Fr1**2))

        st.write(f"Fr₁ = {Fr1:.2f}")
        st.write(f"y₂ ≈ {y2_true:.3f} m")

        if y2_true != 0:
            err = abs(y2_user - y2_true) / abs(y2_true)
            if err < 0.05:
                st.success("y₂: within 5% ✔")
            elif err < 0.15:
                st.warning("y₂: within 15% — okay, but refine algebra.")
            else:
                st.error("y₂: large error — re-check the jump formula.")


def bernoulli_venturi_practice():
    st.subheader("Bernoulli Practice – Venturi Flow Rate")

    if "bv_prob" not in st.session_state:
        st.session_state.bv_prob = None

    if st.button("New random problem", key="bv_new"):
        D1 = random.choice([0.10, 0.12])
        D2 = random.choice([0.05, 0.06])
        dp = random.choice([5000.0, 8000.0, 12000.0])
        rho = 1000.0
        st.session_state.bv_prob = (D1, D2, dp, rho)

    prob = st.session_state.bv_prob
    if prob is None:
        st.info("Click 'New random problem' to generate a scenario.")
        return

    D1, D2, dp, rho = prob
    A1 = math.pi * D1**2 / 4
    A2 = math.pi * D2**2 / 4
    numerator = 2 * dp / rho
    denom = (1.0 / A2**2) - (1.0 / A1**2)
    Q_true = math.sqrt(numerator / denom)

    st.write(
        wrap(
            f"Water (ρ = {rho} kg/m³) flows through a horizontal venturi. "
            f"D₁ = {D1} m, D₂ = {D2} m, and pressure drop Δp = {dp:.0f} Pa. "
            "Neglecting losses, find the flow rate Q."
        )
    )

    Q_user = st.number_input("Your Q [m³/s]", value=0.0, key="bv_Q_user")

    if st.button("Check answer", key="bv_check"):
        st.write(f"Correct Q ≈ {Q_true:.4f} m³/s")
        if Q_true != 0:
            err = abs(Q_user - Q_true) / abs(Q_true)
            if err < 0.05:
                st.success("Q: within 5% ✔")
            elif err < 0.15:
                st.warning("Q: within 15% — okay.")
            else:
                st.error("Q: large error — revisit Bernoulli + continuity.")


def practice_page():
    st.header("🧪 Guided Practice")
    mode = st.selectbox(
        "Choose a practice type",
        [
            "Hydrostatics – vertical gate",
            "Pipe flow – head loss & regime",
            "Hydraulic jump – conjugate depth",
            "Bernoulli – venturi flow rate",
        ],
    )

    if mode == "Hydrostatics – vertical gate":
        hydrostatic_practice()
    elif mode == "Pipe flow – head loss & regime":
        pipe_flow_practice()
    elif mode == "Hydraulic jump – conjugate depth":
        hydraulic_jump_practice()
    elif mode == "Bernoulli – venturi flow rate":
        bernoulli_venturi_practice()


# ------------------ Quiz Page ------------------ #

QUIZ_BANK = {
    "Hydrostatics": [
        {
            "q": "In a static fluid of constant density, which statement is true?",
            "options": [
                "Pressure is constant everywhere.",
                "Pressure increases linearly with depth.",
                "Pressure is higher at the free surface than at the bottom.",
                "Pressure depends only on container shape, not depth.",
            ],
            "answer": 1,
            "expl": "For constant density, dp/dz = -ρg, so p increases linearly with depth.",
        },
        {
            "q": "The center of pressure on a vertical plane surface is:",
            "options": [
                "Always at the centroid.",
                "Always above the centroid.",
                "Always below the centroid.",
                "Independent of the centroid.",
            ],
            "answer": 2,
            "expl": "Because pressure increases with depth, the resultant force acts below the centroid.",
        },
    ],
    "Bernoulli & CV": [
        {
            "q": "Bernoulli’s equation in its simplest form assumes:",
            "options": [
                "Inviscid, incompressible, steady flow along a streamline.",
                "Viscous and compressible flow with pumps.",
                "Unsteady turbulent flow with heat transfer.",
                "Any flow as long as mass is conserved.",
            ],
            "answer": 0,
            "expl": "Classical Bernoulli is derived for inviscid, incompressible, steady flow along a streamline with no shaft work and negligible losses.",
        },
        {
            "q": "In a control volume, continuity for a steady incompressible flow implies:",
            "options": [
                "Sum of pressures in = sum of pressures out.",
                "Sum of velocities in = sum of velocities out.",
                "Sum of mass flow rates in = sum of mass flow rates out.",
                "Sum of energies in = sum of energies out.",
            ],
            "answer": 2,
            "expl": "Continuity enforces mass conservation: total mass flow in equals total mass flow out in steady state.",
        },
    ],
    "Internal Flow": [
        {
            "q": "For laminar flow in a circular pipe, the Darcy friction factor f is:",
            "options": [
                "Constant and independent of Re.",
                "Proportional to Re.",
                "Given by 64/Re.",
                "Found only from the Moody chart.",
            ],
            "answer": 2,
            "expl": "For fully developed laminar flow in a pipe, f = 64/Re exactly.",
        },
        {
            "q": "In turbulent pipe flow, the friction factor f depends primarily on:",
            "options": [
                "Only the pipe diameter.",
                "Only the pipe length.",
                "Reynolds number and relative roughness ε/D.",
                "Fluid temperature only.",
            ],
            "answer": 2,
            "expl": "The Moody chart relates f to both Re and ε/D in turbulent flow.",
        },
    ],
}


def quiz_page():
    st.header("📝 Quick Quizzes")

    topic = st.selectbox("Pick a quiz topic", list(QUIZ_BANK.keys()))

    if "quiz_state" not in st.session_state or st.session_state.quiz_state.get("topic") != topic:
        st.session_state.quiz_state = {"topic": topic, "index": 0, "score": 0, "answered": False}

    qlist = QUIZ_BANK[topic]
    idx = st.session_state.quiz_state["index"]
    score = st.session_state.quiz_state["score"]
    qdata = qlist[idx]

    st.subheader(f"Question {idx + 1} of {len(qlist)}")
    st.write(qdata["q"])
    choice = st.radio("Select one:", qdata["options"], index=None)

    if st.button("Check answer", key="quiz_check"):
        if choice is None:
            st.warning("Choose an option first.")
        else:
            correct_idx = qdata["answer"]
            correct_opt = qdata["options"][correct_idx]
            if choice == correct_opt:
                st.success("Correct ✅")
                st.session_state.quiz_state["score"] += 1
            else:
                st.error(f"Not quite. Correct answer: '{correct_opt}'.")
            st.info(qdata["expl"])
            st.session_state.quiz_state["answered"] = True

    col1, col2 = st.columns(2)
    with col1:
        if st.button("Next question"):
            st.session_state.quiz_state["index"] = (idx + 1) % len(qlist)
            st.session_state.quiz_state["answered"] = False

    with col2:
        st.write(f"Score so far: **{score} / {len(qlist)}**")


# ------------------ Exam Mode ------------------ #

def exam_mode_page():
    st.header("📘 Exam Mode – Concept Checklist")

    st.write(
        wrap(
            "Use these as a mental checklist the night before an exam. "
            "If you can explain each one without looking, you are gooddddd."
        )
    )

    with st.expander("Hydrostatics & Buoyancy"):
        st.markdown(
            """
- Explain why pressure increases with depth mathematically and physically.
- Derive the expression for resultant force on a vertical rectangular gate.
- State Archimedes’ principle and apply it to a floating body.
- Explain the difference between stable and unstable equilibrium for a floating body.
"""
        )

    with st.expander("Control Volume & Bernoulli"):
        st.markdown(
            """
- Write the steady incompressible continuity equation for a control volume.
- Explain each term in the energy equation and when Bernoulli is valid.
- Use the momentum equation to find the force on a pipe bend.
- Explain physically why velocity and pressure trade off in Bernoulli.
"""
        )

    with st.expander("Internal Flow & Losses"):
        st.markdown(
            """
- Define Reynolds number and typical thresholds for laminar/turbulent flow.
- For laminar pipe flow, sketch the velocity profile and recall f = 64/Re.
- Use the Darcy–Weisbach equation to compute head loss.
- Combine major and minor losses in a simple pipe network.
"""
        )

    with st.expander("Dimensional Analysis & External Flow"):
        st.markdown(
            """
- Use Buckingham π to derive dimensionless groups for drag on a sphere.
- Explain geometric, kinematic, and dynamic similarity.
- Define drag coefficient and interpret C_D vs Re curves.
- Describe qualitatively what a boundary layer is and why separation matters.
"""
        )


# ------------------ Main ------------------ #

def main():
    st.sidebar.title("Fluids With Dushto")
    st.sidebar.markdown(
        "A quiet space to make ENME-style fluids easyyyyy."
    )

    page = st.sidebar.radio(
        "Go to",
        ["Welcome", "Concept Tutor", "Calculators", "Visuals", "Guided Practice", "Quizzes", "Exam Mode"],
    )

    st.sidebar.markdown("---")
    st.sidebar.markdown(
        "Study tip: pick *one* idea, read the intuition, then try a problem. "
        "Small steps are enough."
    )

    if page == "Welcome":
        welcome_page()
    elif page == "Concept Tutor":
        concept_tutor_page()
    elif page == "Calculators":
        calculators_page()
    elif page == "Visuals":
        visuals_page()
    elif page == "Guided Practice":
        practice_page()
    elif page == "Quizzes":
        quiz_page()
    elif page == "Exam Mode":
        exam_mode_page()


if __name__ == "__main__":
    main()
