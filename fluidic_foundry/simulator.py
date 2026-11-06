"""
Lagrangian droplet transport engine and dynamic bifurcation logic simulator.
Couples discrete droplet advection with dynamic network resistance feedback.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from .network import HydrodynamicNetwork


@dataclass
class Droplet:
    droplet_id: int
    edge_id: str
    position: float       # Position along edge (meters, 0.0 to edge.geometry.length)
    length: float         # Slug length (meters)
    volume: float         # Volume (m^3)
    payload: str = "bit"  # Logical payload or chemical marker
    color: str = "#38bdf8"
    active: bool = True
    reached_outlet: Optional[str] = None
    exit_time: Optional[float] = None


class FluidicSimulator:
    """
    Coupled two-phase fluidic logic simulation.
    """
    def __init__(
        self,
        network: HydrodynamicNetwork,
        droplet_resistance_ratio: float = 2.5,
        slip_factor: float = 1.1,
    ):
        self.network = network
        self.droplet_resistance_ratio = droplet_resistance_ratio
        self.slip_factor = slip_factor

        self.droplets: List[Droplet] = []
        self.next_droplet_id: int = 1
        self.current_time: float = 0.0

        # Event tracking: outlet_node_id -> list of (droplet_id, exit_time, payload)
        self.outlet_events: Dict[str, List[Tuple[int, float, str]]] = {}

        # Initial baseline solve
        self._sync_network_resistances()
        self.network.solve_flow()

    def spawn_droplet(
        self,
        edge_id: str,
        length: float,
        position: float = 0.0,
        payload: str = "bit",
        color: str = "#38bdf8",
    ) -> Droplet:
        """Inject a droplet into a specified channel."""
        if edge_id not in self.network.edges:
            raise KeyError(f"Channel '{edge_id}' does not exist.")
        edge = self.network.edges[edge_id]
        volume = length * edge.geometry.cross_section_area
        droplet = Droplet(
            droplet_id=self.next_droplet_id,
            edge_id=edge_id,
            position=min(edge.geometry.length, max(0.0, position)),
            length=length,
            volume=volume,
            payload=payload,
            color=color,
            active=True,
        )
        self.next_droplet_id += 1
        self.droplets.append(droplet)
        self._sync_network_resistances()
        self.network.solve_flow()
        return droplet

    def step(self, dt: float):
        """
        Advance physics by time step dt (seconds).
        """
        if dt <= 0:
            return

        # 1. Solve current flow field
        self._sync_network_resistances()
        self.network.solve_flow()

        # 2. Advect each active droplet
        for drop in self.droplets:
            if not drop.active:
                continue

            edge = self.network.edges[drop.edge_id]
            u_mean = edge.flow_rate / edge.geometry.cross_section_area
            v_drop = self.slip_factor * u_mean

            drop.position += v_drop * dt

            # Check boundary crossings
            if drop.position >= edge.geometry.length:
                self._handle_node_arrival(drop, edge.to_node)
            elif drop.position < 0.0:
                self._handle_node_arrival(drop, edge.from_node)

        self.current_time += dt

    def run_until(self, target_time: float, max_dt: float = 0.0005):
        """Run simulation until reaching target_time."""
        while self.current_time < target_time:
            dt = min(max_dt, target_time - self.current_time)
            self.step(dt)

    def _sync_network_resistances(self):
        """Update each channel's effective resistance based on resident droplets."""
        edge_drop_counts: Dict[str, int] = {eid: 0 for eid in self.network.edges}
        for drop in self.droplets:
            if drop.active:
                edge_drop_counts[drop.edge_id] = edge_drop_counts.get(drop.edge_id, 0) + 1

        for eid, count in edge_drop_counts.items():
            edge = self.network.edges[eid]
            drop_len = 0.0
            if count > 0:
                lengths = [d.length for d in self.droplets if d.active and d.edge_id == eid]
                drop_len = sum(lengths) / len(lengths)
            self.network.update_channel_resistance(
                eid,
                droplet_count=count,
                droplet_length=drop_len,
                ratio=self.droplet_resistance_ratio,
            )

    def _handle_node_arrival(self, drop: Droplet, node_id: str):
        """
        Handle droplet transition at a junction/node.
        Applies hydrodynamic minimum-resistance / maximum-outflow bifurcation selection.
        """
        candidate_edges = []
        for eid in self.network.adjacency[node_id]:
            if eid == drop.edge_id:
                continue
            edge = self.network.edges[eid]
            if edge.from_node == node_id and edge.flow_rate > 1e-15:
                candidate_edges.append((eid, edge.flow_rate))
            elif edge.to_node == node_id and edge.flow_rate < -1e-15:
                candidate_edges.append((eid, -edge.flow_rate))

        if not candidate_edges:
            # Droplet reached terminal sink
            drop.active = False
            drop.reached_outlet = node_id
            drop.exit_time = self.current_time
            if node_id not in self.outlet_events:
                self.outlet_events[node_id] = []
            self.outlet_events[node_id].append((drop.droplet_id, self.current_time, drop.payload))
            # Sync resistance since droplet left
            self._sync_network_resistances()
            self.network.solve_flow()
            return

        # Microfluidic bifurcation rule:
        # Droplet enters branch with MAXIMUM forward volumetric flow rate
        candidate_edges.sort(key=lambda item: item[1], reverse=True)
        chosen_edge_id = candidate_edges[0][0]

        drop.edge_id = chosen_edge_id
        chosen_edge = self.network.edges[chosen_edge_id]
        if chosen_edge.from_node == node_id:
            drop.position = 0.0
        else:
            drop.position = chosen_edge.geometry.length

        # Immediately update flow field so trailing droplets perceive new hydrodynamic resistance
        self._sync_network_resistances()
        self.network.solve_flow()
