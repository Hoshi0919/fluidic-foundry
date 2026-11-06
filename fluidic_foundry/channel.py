"""
Microchannel geometry and hydrodynamic resistance modeling.
Implements Hagen-Poiseuille rectangular duct flow resistance and
two-phase Bretherton lubrication droplet excess resistance.
"""

import math
from dataclasses import dataclass
from typing import Optional


@dataclass
class ChannelGeometry:
    """
    Rectangular microfluidic channel parameters.
    Dimensions in SI units (meters).
    """
    length: float  # Channel length L (m)
    width: float   # Channel width w (m)
    height: float  # Channel height h (m)

    @property
    def cross_section_area(self) -> float:
        """Cross-sectional area A = w * h (m^2)."""
        return self.width * self.height

    @property
    def hydraulic_diameter(self) -> float:
        """Hydraulic diameter D_h = 4A / P = 2*w*h / (w+h) (m)."""
        return 2.0 * self.width * self.height / (self.width + self.height)

    @property
    def aspect_ratio(self) -> float:
        """Aspect ratio alpha = min(w, h) / max(w, h) <= 1.0."""
        return min(self.width, self.height) / max(self.width, self.height)


def poiseuille_resistance(
    geom: ChannelGeometry,
    viscosity: float = 1.0e-3,  # Dynamic viscosity mu of continuous phase (Pa*s, water ~ 1e-3)
    series_terms: int = 10,
) -> float:
    """
    Calculate exact Hagen-Poiseuille hydrodynamic flow resistance for a rectangular channel.
    R_hyd = Delta_P / Q (Unit: Pa * s / m^3).

    Analytical solution from Navier-Stokes creeping flow:
    R = (12 * mu * L) / (h_min^3 * w_max * (1 - (192 / pi^5) * alpha * sum_{n=1,3,...} (1/n^5) * tanh(n*pi / (2*alpha))))
    where alpha = h_min / w_max <= 1.
    """
    if geom.length <= 0 or geom.width <= 0 or geom.height <= 0:
        raise ValueError("Channel dimensions must be strictly positive.")

    w = max(geom.width, geom.height)
    h = min(geom.width, geom.height)
    alpha = h / w

    # Fourier series term
    sum_val = 0.0
    for i in range(series_terms):
        n = 2 * i + 1
        term = (1.0 / (n**5)) * math.tanh((n * math.pi) / (2.0 * alpha))
        sum_val += term

    geometric_factor = 1.0 - (192.0 / (math.pi**5)) * alpha * sum_val
    if geometric_factor <= 0:
        # Fallback to classical 0.630 approximation if numerical edge
        geometric_factor = max(1e-5, 1.0 - 0.630 * alpha)

    resistance = (12.0 * viscosity * geom.length) / (w * (h**3) * geometric_factor)
    return resistance


def droplet_excess_resistance(
    base_resistance: float,
    droplet_length: float,
    channel_length: float,
    resistance_ratio: float = 1.5,
) -> float:
    """
    Calculate additional hydrodynamic resistance imparted by a confined droplet/bubble.
    Delta R_drop = resistance_ratio * (droplet_length / channel_length) * base_resistance.
    
    In microfluidic bubble logic (Prakash & Gershenfeld 2007), a single bubble
    entering a branch typically increases branch resistance by 100% to 250% (ratio 1.0~2.5),
    inverting the minimum resistance path for following droplets.
    """
    if channel_length <= 0:
        return 0.0
    confinement = min(1.0, max(0.0, droplet_length / channel_length))
    return base_resistance * resistance_ratio * confinement
