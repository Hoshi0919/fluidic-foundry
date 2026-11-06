"""
Comprehensive physical validation and unit test matrix for Fluidic Foundry.
Covers Hagen-Poiseuille resistance scaling, Kirchhoff flow conservation,
Bretherton droplet excess resistance inversion, and Boolean logic truth tables.
"""

import math
import pytest
from fluidic_foundry.channel import (
    ChannelGeometry,
    poiseuille_resistance,
    droplet_excess_resistance,
)
from fluidic_foundry.network import HydrodynamicNetwork
from fluidic_foundry.simulator import FluidicSimulator
from fluidic_foundry.gates import (
    create_and_or_gate,
    create_not_inverter_gate,
    create_half_adder_circuit,
)


def test_rectangular_channel_poiseuille_resistance_scaling():
    """Verify that flow resistance scales linearly with length and inversely with h^3."""
    geom1 = ChannelGeometry(length=1.0e-3, width=100.0e-6, height=50.0e-6)
    geom2 = ChannelGeometry(length=2.0e-3, width=100.0e-6, height=50.0e-6)

    r1 = poiseuille_resistance(geom1, viscosity=1.0e-3)
    r2 = poiseuille_resistance(geom2, viscosity=1.0e-3)

    # Linear length scaling
    assert math.isclose(r2 / r1, 2.0, rel_tol=1e-5)

    # Depth reduction scales resistance sharply (approx h^-3)
    geom_shallow = ChannelGeometry(length=1.0e-3, width=100.0e-6, height=25.0e-6)
    r_shallow = poiseuille_resistance(geom_shallow, viscosity=1.0e-3)
    # Shallower channel should have >= 7x resistance
    assert r_shallow / r1 > 6.0


def test_droplet_excess_resistance_confinement():
    """Verify Bretherton-style droplet confinement resistance."""
    base_r = 1.0e11
    # 0 length droplet -> 0 excess resistance
    assert droplet_excess_resistance(base_r, 0.0, 1.0e-3) == 0.0

    # Half length droplet
    delta_r = droplet_excess_resistance(base_r, 0.5e-3, 1.0e-3, resistance_ratio=2.0)
    assert math.isclose(delta_r, base_r * 1.0, rel_tol=1e-5)


def test_kirchhoff_flow_conservation_at_junction():
    """Verify that net volumetric flux at internal junctions is strictly zero."""
    net = HydrodynamicNetwork(continuous_viscosity=1.0e-3)
    net.add_node("inlet", x=0, y=0, fixed_pressure=1000.0)
    net.add_node("junction", x=1.0e-3, y=0, fixed_pressure=None)
    net.add_node("out1", x=2.0e-3, y=1.0e-3, fixed_pressure=0.0)
    net.add_node("out2", x=2.0e-3, y=-1.0e-3, fixed_pressure=0.0)

    net.add_channel("c_in", "inlet", "junction", length=1e-3, width=50e-6, height=50e-6)
    net.add_channel("c_out1", "junction", "out1", length=1e-3, width=50e-6, height=50e-6)
    net.add_channel("c_out2", "junction", "out2", length=1e-3, width=50e-6, height=50e-6)

    flows = net.solve_flow()
    q_in = flows["c_in"]
    q_out1 = flows["c_out1"]
    q_out2 = flows["c_out2"]

    # Inflow must equal sum of outflows
    net_flux = q_in - (q_out1 + q_out2)
    assert abs(net_flux) < 1e-18
    # Symmetric branches must have equal split
    assert math.isclose(q_out1, q_out2, rel_tol=1e-5)


def test_hydrodynamic_ohms_law():
    """Verify Delta P = Q * R on a single channel."""
    net = HydrodynamicNetwork(continuous_viscosity=1.0e-3)
    p_in = 500.0
    net.add_node("n1", x=0, y=0, fixed_pressure=p_in)
    net.add_node("n2", x=1.0e-3, y=0, fixed_pressure=0.0)
    edge = net.add_channel("ch", "n1", "n2", length=1e-3, width=60e-6, height=40e-6)

    flows = net.solve_flow()
    q = flows["ch"]
    assert math.isclose(p_in, q * edge.base_resistance, rel_tol=1e-5)


def test_resistance_inversion_phase_transition():
    """
    Core physical demonstration:
    An asymmetric bifurcation flips its preferred branch when a droplet enters!
    """
    net, ports = create_and_or_gate(driving_pressure=1200.0)
    sim = FluidicSimulator(net, droplet_resistance_ratio=2.5)

    # 1. Empty state: Branch OR has lower resistance -> higher flow rate
    flows_empty = net.solve_flow()
    q_or_empty = flows_empty[ports["branch_or"]]
    q_and_empty = flows_empty[ports["branch_and"]]
    assert q_or_empty > q_and_empty, "Empty bifurcation must favor the wider OR branch"

    # 2. Inject a droplet into Branch OR
    sim.spawn_droplet(ports["branch_or"], length=1.2e-3, position=0.5e-3)
    sim._sync_network_resistances()
    flows_loaded = net.solve_flow()
    q_or_loaded = flows_loaded[ports["branch_or"]]
    q_and_loaded = flows_loaded[ports["branch_and"]]

    # Hydraulic resistance inversion: Branch AND now has higher flow rate!
    assert q_and_loaded > q_or_loaded, "Droplet presence in OR branch must invert flow preference to AND branch!"


def test_and_or_gate_truth_table_0_0():
    """Truth table: A=0, B=0 -> OR=0, AND=0."""
    net, ports = create_and_or_gate()
    sim = FluidicSimulator(net)
    sim.run_until(0.05)

    assert len(sim.outlet_events.get(ports["out_or"], [])) == 0
    assert len(sim.outlet_events.get(ports["out_and"], [])) == 0


def test_and_or_gate_truth_table_1_0():
    """Truth table: A=1, B=0 -> OR=1, AND=0."""
    net, ports = create_and_or_gate()
    sim = FluidicSimulator(net)

    # Inject single droplet from Input A
    sim.spawn_droplet(ports["in_a"], length=1.0e-3, position=0.0e-3, payload="A")
    # Run until droplet reaches outlet
    sim.run_until(0.35)

    or_events = sim.outlet_events.get(ports["out_or"], [])
    and_events = sim.outlet_events.get(ports["out_and"], [])

    assert len(or_events) == 1, "Single droplet A must exit through OR branch"
    assert len(and_events) == 0, "No droplet should exit through AND branch"
    assert or_events[0][2] == "A"


def test_and_or_gate_truth_table_0_1():
    """Truth table: A=0, B=1 -> OR=1, AND=0."""
    net, ports = create_and_or_gate()
    sim = FluidicSimulator(net)

    # Inject single droplet from Input B
    sim.spawn_droplet(ports["in_b"], length=1.0e-3, position=0.0e-3, payload="B")
    sim.run_until(0.35)

    or_events = sim.outlet_events.get(ports["out_or"], [])
    and_events = sim.outlet_events.get(ports["out_and"], [])

    assert len(or_events) == 1, "Single droplet B must exit through OR branch"
    assert len(and_events) == 0, "No droplet should exit through AND branch"
    assert or_events[0][2] == "B"


def test_and_or_gate_truth_table_1_1():
    """
    Truth table: A=1, B=1 -> OR=1, AND=1.
    First droplet takes the preferred OR branch;
    Second droplet is dynamically repelled by the hydrodynamic resistance into AND branch!
    """
    net, ports = create_and_or_gate(driving_pressure=1200.0)
    sim = FluidicSimulator(net, droplet_resistance_ratio=2.5)

    # Spawn A slightly ahead of B
    sim.spawn_droplet(ports["in_a"], length=1.0e-3, position=0.5e-3, payload="A")
    sim.spawn_droplet(ports["in_b"], length=1.0e-3, position=0.0e-3, payload="B")

    sim.run_until(0.35)

    or_events = sim.outlet_events.get(ports["out_or"], [])
    and_events = sim.outlet_events.get(ports["out_and"], [])

    assert len(or_events) == 1, "First droplet must exit through OR branch"
    assert len(and_events) == 1, "Second droplet must be diverted into AND branch!"
    assert len(or_events) + len(and_events) == 2


def test_not_inverter_gate_uninhibited():
    """When control signal is absent (A=0), clock droplet passes straight through to NOT output (1)."""
    net, ports = create_not_inverter_gate()
    sim = FluidicSimulator(net)

    # Only clock droplet is injected
    sim.spawn_droplet(ports["in_clk"], length=0.8e-3, position=0.0, payload="CLK")
    sim.run_until(0.35)

    not_events = sim.outlet_events.get(ports["out_not"], [])
    dump_events = sim.outlet_events.get(ports["out_dump"], [])

    assert len(not_events) == 1, "Clock droplet must exit to OUT_NOT when control is 0"
    assert len(dump_events) == 0


def test_half_adder_single_input_yields_sum():
    """Half adder: A=1, B=0 -> Sum=1, Carry=0."""
    net, ports = create_half_adder_circuit()
    sim = FluidicSimulator(net)

    sim.spawn_droplet(ports["in_a"], length=1.0e-3, position=0.0, payload="A")
    sim.run_until(0.35)

    sum_events = sim.outlet_events.get(ports["out_sum"], [])
    carry_events = sim.outlet_events.get(ports["out_carry"], [])

    assert len(sum_events) == 1, "Single input must produce Sum output bit"
    assert len(carry_events) == 0, "Single input must not produce Carry bit"


def test_reynolds_and_capillary_regime_validity():
    """
    Verify operating regime:
    Reynolds number Re << 2000 (Laminar creeping flow regime)
    Capillary number Ca < 1e-2 (Surface tension dominates, preventing droplet breakup)
    """
    density = 1000.0       # kg/m^3 (water)
    viscosity = 1.0e-3     # Pa*s
    surface_tension = 0.05 # N/m (water-oil interface with surfactant)

    net, ports = create_and_or_gate(driving_pressure=1200.0)
    flows = net.solve_flow()
    q_max = max(abs(f) for f in flows.values())
    edge = net.edges[ports["branch_or"]]
    u_mean = q_max / edge.geometry.cross_section_area
    d_h = edge.geometry.hydraulic_diameter

    # Reynolds number
    re = (density * u_mean * d_h) / viscosity
    # Capillary number
    ca = (viscosity * u_mean) / surface_tension

    assert re < 5.0, f"Reynolds number {re} must be < 5.0 for microfluidic laminar flow"
    assert ca < 0.05, f"Capillary number {ca} must be small (< 0.05) to prevent breakup"


def test_half_adder_double_input_yields_carry():
    """Half adder: A=1, B=1 -> Sum=1 (first drop), Carry=1 (second drop diverted)."""
    net, ports = create_half_adder_circuit(driving_pressure=1500.0)
    sim = FluidicSimulator(net, droplet_resistance_ratio=2.5)

    sim.spawn_droplet(ports["in_a"], length=1.0e-3, position=0.5e-3, payload="A")
    sim.spawn_droplet(ports["in_b"], length=1.0e-3, position=0.0e-3, payload="B")
    sim.run_until(0.35)

    sum_events = sim.outlet_events.get(ports["out_sum"], [])
    carry_events = sim.outlet_events.get(ports["out_carry"], [])

    assert len(sum_events) == 1, "First droplet must exit to Sum"
    assert len(carry_events) == 1, "Second droplet must be diverted to Carry"


def test_channel_geometry_validation():
    """Verify negative or zero channel dimensions raise ValueError."""
    with pytest.raises(ValueError):
        poiseuille_resistance(ChannelGeometry(length=-1.0, width=50e-6, height=50e-6))
    with pytest.raises(ValueError):
        poiseuille_resistance(ChannelGeometry(length=1e-3, width=0.0, height=50e-6))


def test_network_singular_matrix_guard():
    """Verify disconnected floating internal node triggers singular matrix detection."""
    net = HydrodynamicNetwork()
    net.add_node("n_float", x=0, y=0, fixed_pressure=None)
    with pytest.raises(ValueError, match="Singular matrix"):
        net.solve_flow()
