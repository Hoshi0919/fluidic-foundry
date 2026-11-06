"""
Fluidic Foundry: Microfluidic Bubble Logic & Capillary Computing Simulation Engine.
Based on low-Reynolds-number two-phase hydrodynamic resistance switching (Prakash & Gershenfeld, Science 2007).
"""

from .channel import ChannelGeometry, poiseuille_resistance, droplet_excess_resistance
from .network import HydrodynamicNetwork, FluidicNode, ChannelEdge
from .simulator import FluidicSimulator, Droplet
from .gates import create_and_or_gate, create_not_inverter_gate, create_half_adder_circuit

__version__ = "1.0.0"
__author__ = "Hoshi <hoshi0919@outlook.com>"

__all__ = [
    "ChannelGeometry",
    "poiseuille_resistance",
    "droplet_excess_resistance",
    "HydrodynamicNetwork",
    "FluidicNode",
    "ChannelEdge",
    "FluidicSimulator",
    "Droplet",
    "create_and_or_gate",
    "create_not_inverter_gate",
    "create_half_adder_circuit",
]
