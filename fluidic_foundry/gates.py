"""
Standard microfluidic logic gate circuits and canonical topologies.
Constructs AND/OR gates, NOT gates, and Half Adders.
"""

from typing import Tuple, Dict, Any
from .network import HydrodynamicNetwork


def create_and_or_gate(
    driving_pressure: float = 1200.0, # Inlet pressure in Pascals (~1.2 kPa)
    channel_height: float = 50.0e-6,   # 50 micrometers
    base_width: float = 50.0e-6,       # 50 micrometers
) -> Tuple[HydrodynamicNetwork, Dict[str, str]]:
    """
    Constructs an asymmetric bifurcation AND/OR logic gate (Prakash & Gershenfeld 2007).
    
    Topology:
      [In_A] ---> (J_in) --- Branch_OR (R1, lower resistance) ----> [Out_OR]
      [In_B] ---/        --- Branch_AND (R2, higher resistance) ---> [Out_AND]
    
    Optimal sizing:
      Branch OR width = 1.15 * base_width -> R1 ~ 6.95e12
      Branch AND width = 0.95 * base_width -> R2 ~ 1.01e13
      Empty: R2 / R1 = 1.45 (OR branch dominates)
      Loaded with droplet: R1_eff / R2 = 1.55 (AND branch dominates)
    """
    net = HydrodynamicNetwork(continuous_viscosity=1.0e-3)

    # Coordinates (meters)
    net.add_node("in_a", x=0.0e-3, y=1.0e-3, fixed_pressure=driving_pressure)
    net.add_node("in_b", x=0.0e-3, y=-1.0e-3, fixed_pressure=driving_pressure)
    net.add_node("j_in", x=1.5e-3, y=0.0e-3, fixed_pressure=None)
    net.add_node("out_or", x=3.5e-3, y=0.8e-3, fixed_pressure=0.0)
    net.add_node("out_and", x=3.5e-3, y=-0.8e-3, fixed_pressure=0.0)

    # Input channels
    net.add_channel(
        "ch_in_a", "in_a", "j_in",
        length=1.5e-3, width=base_width, height=channel_height
    )
    net.add_channel(
        "ch_in_b", "in_b", "j_in",
        length=1.5e-3, width=base_width, height=channel_height
    )

    # Branch OR: Wider (w = 1.15 * base) -> lower resistance
    net.add_channel(
        "ch_branch_or", "j_in", "out_or",
        length=2.0e-3, width=base_width * 1.15, height=channel_height
    )

    # Branch AND: Narrower (w = 0.95 * base) -> higher resistance
    net.add_channel(
        "ch_branch_and", "j_in", "out_and",
        length=2.0e-3, width=base_width * 0.95, height=channel_height
    )

    ports = {
        "in_a": "ch_in_a",
        "in_b": "ch_in_b",
        "out_or": "out_or",
        "out_and": "out_and",
        "branch_or": "ch_branch_or",
        "branch_and": "ch_branch_and",
    }
    return net, ports


def create_not_inverter_gate(
    driving_pressure: float = 1200.0,
    channel_height: float = 50.0e-6,
    base_width: float = 50.0e-6,
) -> Tuple[HydrodynamicNetwork, Dict[str, str]]:
    """
    Constructs an all-fluidic NOT inverter gate with a continuous clock stream.
    """
    net = HydrodynamicNetwork(continuous_viscosity=1.0e-3)

    # Nodes
    net.add_node("in_clk", x=0.0e-3, y=0.8e-3, fixed_pressure=driving_pressure)
    net.add_node("in_ctrl", x=0.0e-3, y=-0.8e-3, fixed_pressure=driving_pressure)
    net.add_node("j_clk", x=1.5e-3, y=0.8e-3, fixed_pressure=None)
    net.add_node("j_ctrl", x=1.5e-3, y=-0.8e-3, fixed_pressure=None)
    net.add_node("out_not", x=3.5e-3, y=1.2e-3, fixed_pressure=0.0)
    net.add_node("out_dump", x=3.5e-3, y=0.3e-3, fixed_pressure=0.0)
    net.add_node("out_waste", x=3.5e-3, y=-0.8e-3, fixed_pressure=0.0)

    # Clock inlet channel
    net.add_channel("ch_clk", "in_clk", "j_clk", length=1.5e-3, width=base_width, height=channel_height)
    # Control inlet channel
    net.add_channel("ch_ctrl", "in_ctrl", "j_ctrl", length=1.5e-3, width=base_width, height=channel_height)

    # Cross-bridge coupling channel
    net.add_channel("ch_bridge", "j_ctrl", "j_clk", length=1.6e-3, width=base_width * 0.7, height=channel_height)

    # Out NOT: favored under default clock
    net.add_channel("ch_out_not", "j_clk", "out_not", length=2.0e-3, width=base_width * 1.15, height=channel_height)
    # Out Dump: higher resistance
    net.add_channel("ch_out_dump", "j_clk", "out_dump", length=2.0e-3, width=base_width * 0.95, height=channel_height)
    # Control waste
    net.add_channel("ch_out_waste", "j_ctrl", "out_waste", length=2.0e-3, width=base_width, height=channel_height)

    ports = {
        "in_clk": "ch_clk",
        "in_ctrl": "ch_ctrl",
        "out_not": "out_not",
        "out_dump": "out_dump",
        "out_waste": "out_waste",
    }
    return net, ports


def create_half_adder_circuit(
    driving_pressure: float = 1200.0,
    channel_height: float = 50.0e-6,
    base_width: float = 50.0e-6,
) -> Tuple[HydrodynamicNetwork, Dict[str, str]]:
    """
    Constructs a 1-Bit Fluidic Half Adder.
    Outputs: Sum (XOR / single drop) and Carry (AND / double drop).
    """
    net = HydrodynamicNetwork(continuous_viscosity=1.0e-3)

    net.add_node("in_a", x=0.0e-3, y=1.2e-3, fixed_pressure=driving_pressure)
    net.add_node("in_b", x=0.0e-3, y=-1.2e-3, fixed_pressure=driving_pressure)
    net.add_node("j_stage1", x=1.5e-3, y=0.0e-3, fixed_pressure=None)
    net.add_node("j_sum", x=3.0e-3, y=0.8e-3, fixed_pressure=None)
    net.add_node("j_carry", x=3.0e-3, y=-0.8e-3, fixed_pressure=None)
    net.add_node("out_sum", x=4.5e-3, y=0.8e-3, fixed_pressure=0.0)
    net.add_node("out_carry", x=4.5e-3, y=-0.8e-3, fixed_pressure=0.0)

    # Feed channels
    net.add_channel("ch_in_a", "in_a", "j_stage1", length=1.5e-3, width=base_width, height=channel_height)
    net.add_channel("ch_in_b", "in_b", "j_stage1", length=1.5e-3, width=base_width, height=channel_height)

    # First branch to Sum (lower resistance)
    net.add_channel("ch_to_sum", "j_stage1", "j_sum", length=1.5e-3, width=base_width * 1.15, height=channel_height)
    # Second branch to Carry (higher resistance)
    net.add_channel("ch_to_carry", "j_stage1", "j_carry", length=1.5e-3, width=base_width * 0.95, height=channel_height)

    # Final runs
    net.add_channel("ch_out_sum", "j_sum", "out_sum", length=1.5e-3, width=base_width, height=channel_height)
    net.add_channel("ch_out_carry", "j_carry", "out_carry", length=1.5e-3, width=base_width, height=channel_height)

    ports = {
        "in_a": "ch_in_a",
        "in_b": "ch_in_b",
        "out_sum": "out_sum",
        "out_carry": "out_carry",
    }
    return net, ports
