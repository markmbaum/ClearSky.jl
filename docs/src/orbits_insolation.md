# Orbits and Insolation

A collection of functions are available for working with general elliptical orbits and insolation patterns for planets orbiting stars.

Function arguments are defined below

| Argument | Definition | Units |
| -------: | :--------- | :---- |
| `a` | [semi-major axis](https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes) | m |
| `e` | [eccentricity](https://en.wikipedia.org/wiki/Orbital_eccentricity) | dimensionless |
| `E` | [eccentric anomaly](https://en.wikipedia.org/wiki/Eccentric_anomaly) | radians |
| `f` | [true anomaly](https://en.wikipedia.org/wiki/True_anomaly) or stellar longitude | radians |
| `γ` | [obliquity](https://en.wikipedia.org/wiki/Axial_tilt) | radians |
| `m` | star mass | kg |
| `p` | [precession angle](https://en.wikipedia.org/wiki/Axial_precession) | radians |
| `ϕ` | latitude | radians |
| `ϕₛ` | substellar latitude | radians |
| `rₐ` | apoapsis distance | m |
| `rₚ` | periapsis distance | m |
| `t` | time | seconds |
| `T` | orbital period | seconds |

The mass `m`, most precisely, should be the sum of the star mass and planet mass, ``m_s + m_p``. For most cases the planet mass is negligible, however, and ``m_s + m_p \approx m_s``.

The precession angle `p` is defined so that when ``p=0``, the northern hemisphere is tilted directly toward the star at periapsis. This means that northern summer occurs when planet is closet to the star. Different values of ``p ∈ [0,2π]`` control when in the orbital path the equinoxes and solstices occur.

-----

```@autodocs
Modules = [ClearSky]
Pages   = ["orbital.jl", "insolation.jl"]
```
