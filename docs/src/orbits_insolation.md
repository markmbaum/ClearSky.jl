# Orbits and Insolation

A collection of functions are available for working with general elliptical orbits and insolation patterns for planets orbiting stars.

Function arguments are defined below

| Argument | Definition | Units |
| -------: | :--------- | :---- |
| `a` | [semi-major axis](https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes) | m |
| `e` | [eccentricity](https://en.wikipedia.org/wiki/Orbital_eccentricity) | - |
| `E` | [eccentric anomaly](https://en.wikipedia.org/wiki/Eccentric_anomaly) | rad |
| `f` | [true anomaly](https://en.wikipedia.org/wiki/True_anomaly) or stellar longitude | rad |
| `γ` | [obliquity](https://en.wikipedia.org/wiki/Axial_tilt) | rad |
| `m` | star mass | kg |
| `p` | [precession angle](https://en.wikipedia.org/wiki/Axial_precession) | rad |
| `θ` | latitude | rad |
| `θₛ` | substellar latitude | rad |
| `rₐ` | apoapsis distance | m |
| `rₚ` | periapsis distance | m |
| `t` | time | sec |
| `T` | orbital period | sec |

The mass `m`, most precisely, should be the sum of the star mass and planet mass, ``m_s + m_p``. For most cases the planet mass is negligible, however, and ``m_s + m_p \approx m_s``.

The precession angle `p` is defined so that when ``p=0``, the northern hemisphere is tilted directly toward the star at periapsis. This means that northern summer occurs when planet is closet to the star. Different values of ``p ∈ [0,2π]`` control when in the orbital path the equinoxes and solstices occur. For example, if ``p = π/2``, the vernal equinox occurs at periapsis and the northern hemisphere is moving into summer.

-----

```@autodocs
Modules = [ClearSky]
Pages   = ["orbits.jl", "insolation.jl"]
```
