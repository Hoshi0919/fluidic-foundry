"""
Hydrodynamic flow network solver based on Kirchhoff's laws for fluidic circuits.
Solves node pressures and channel flow rates for arbitrary 2D microfluidic topologies.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
from .channel import ChannelGeometry, poiseuille_resistance, droplet_excess_resistance


@dataclass
class ChannelEdge:
    edge_id: str
    from_node: str
    to_node: str
    geometry: ChannelGeometry
    base_resistance: float
    current_resistance: float
    droplet_count: int = 0
    flow_rate: float = 0.0  # Volumetric flow rate Q (m^3/s). Positive from -> to, negative to -> from.


@dataclass
class FluidicNode:
    node_id: str
    x: float  # Physical coordinate (m)
    y: float
    fixed_pressure: Optional[float] = None  # None for internal nodes, float for boundary ports (Pa)
    external_flux: float = 0.0              # Net injected volumetric flow rate (m^3/s)
    pressure: float = 0.0                   # Calculated pressure (Pa)


class HydrodynamicNetwork:
    """
    Microfluidic circuit graph representing channels and junctions.
    """
    def __init__(self, continuous_viscosity: float = 1.0e-3):
        self.viscosity: float = continuous_viscosity
        self.nodes: Dict[str, FluidicNode] = {}
        self.edges: Dict[str, ChannelEdge] = {}
        self.adjacency: Dict[str, List[str]] = {}  # node_id -> list of edge_ids

    def add_node(
        self,
        node_id: str,
        x: float,
        y: float,
        fixed_pressure: Optional[float] = None,
        external_flux: float = 0.0,
    ) -> FluidicNode:
        if node_id in self.nodes:
            raise ValueError(f"Node '{node_id}' already exists.")
        node = FluidicNode(
            node_id=node_id,
            x=x,
            y=y,
            fixed_pressure=fixed_pressure,
            external_flux=external_flux,
            pressure=fixed_pressure if fixed_pressure is not None else 0.0,
        )
        self.nodes[node_id] = node
        self.adjacency[node_id] = []
        return node

    def add_channel(
        self,
        edge_id: str,
        from_node: str,
        to_node: str,
        length: float,
        width: float,
        height: float,
        resistance_override: Optional[float] = None,
    ) -> ChannelEdge:
        if edge_id in self.edges:
            raise ValueError(f"Channel '{edge_id}' already exists.")
        if from_node not in self.nodes or to_node not in self.nodes:
            raise KeyError(f"Both nodes '{from_node}' and '{to_node}' must exist in network.")

        geom = ChannelGeometry(length=length, width=width, height=height)
        if resistance_override is not None:
            r0 = resistance_override
        else:
            r0 = poiseuille_resistance(geom, self.viscosity)

        edge = ChannelEdge(
            edge_id=edge_id,
            from_node=from_node,
            to_node=to_node,
            geometry=geom,
            base_resistance=r0,
            current_resistance=r0,
        )
        self.edges[edge_id] = edge
        self.adjacency[from_node].append(edge_id)
        self.adjacency[to_node].append(edge_id)
        return edge

    def update_channel_resistance(self, edge_id: str, droplet_count: int, droplet_length: float, ratio: float = 1.8):
        """Update channel resistance due to presence of droplets."""
        edge = self.edges[edge_id]
        edge.droplet_count = droplet_count
        if droplet_count <= 0:
            edge.current_resistance = edge.base_resistance
        else:
            total_drop_len = droplet_count * droplet_length
            delta_r = droplet_excess_resistance(edge.base_resistance, total_drop_len, edge.geometry.length, ratio)
            edge.current_resistance = edge.base_resistance + delta_r

    def solve_flow(self) -> Dict[str, float]:
        """
        Solve Kirchhoff network equations for unknown node pressures and channel flow rates.
        Returns mapping of edge_id -> flow_rate (m^3/s).
        """
        # Identify internal (unknown pressure) nodes vs Dirichlet boundary nodes
        all_node_ids = list(self.nodes.keys())
        unknown_nodes = [nid for nid in all_node_ids if self.nodes[nid].fixed_pressure is None]
        n_unknown = len(unknown_nodes)

        if n_unknown == 0:
            # All pressures fixed
            self._update_edge_flows()
            return {eid: e.flow_rate for eid, e in self.edges.items()}

        node_to_idx = {nid: i for i, nid in enumerate(unknown_nodes)}

        # Build linear system A * P = b
        A = [[0.0 for _ in range(n_unknown)] for _ in range(n_unknown)]
        b = [0.0 for _ in range(n_unknown)]

        for i, u_id in enumerate(unknown_nodes):
            u_node = self.nodes[u_id]
            b[i] += u_node.external_flux  # External flux into node

            for eid in self.adjacency[u_id]:
                edge = self.edges[eid]
                cond = 1.0 / edge.current_resistance
                v_id = edge.to_node if edge.from_node == u_id else edge.from_node

                A[i][i] += cond
                v_node = self.nodes[v_id]
                if v_node.fixed_pressure is not None:
                    # Boundary node with fixed potential
                    b[i] += cond * v_node.fixed_pressure
                else:
                    # Internal unknown node
                    j = node_to_idx[v_id]
                    A[i][j] -= cond

        # Solve linear system using Gaussian elimination with partial pivoting
        p_sol = self._solve_linear_system(A, b)
        for i, u_id in enumerate(unknown_nodes):
            self.nodes[u_id].pressure = p_sol[i]

        self._update_edge_flows()
        return {eid: e.flow_rate for eid, e in self.edges.items()}

    def _update_edge_flows(self):
        """Compute flow rates on all channels based on resolved node pressures."""
        for edge in self.edges.values():
            pu = self.nodes[edge.from_node].pressure
            pv = self.nodes[edge.to_node].pressure
            edge.flow_rate = (pu - pv) / edge.current_resistance

    @staticmethod
    def _solve_linear_system(A: List[List[float]], b: List[float]) -> List[float]:
        """Gaussian elimination with partial pivoting."""
        n = len(b)
        # Augment matrix
        M = [A[i][:] + [b[i]] for i in range(n)]

        for i in range(n):
            # Pivot
            max_row = i
            max_val = abs(M[i][i])
            for k in range(i + 1, n):
                if abs(M[k][i]) > max_val:
                    max_val = abs(M[k][i])
                    max_row = k
            if max_val < 1e-18:
                raise ValueError("Singular matrix in hydrodynamic network solver.")

            M[i], M[max_row] = M[max_row], M[i]

            # Eliminate below
            pivot = M[i][i]
            for k in range(i + 1, n):
                factor = M[k][i] / pivot
                for j in range(i, n + 1):
                    M[k][j] -= factor * M[i][j]

        # Back-substitution
        x = [0.0 for _ in range(n)]
        for i in range(n - 1, -1, -1):
            s = M[i][n]
            for j in range(i + 1, n):
                s -= M[i][j] * x[j]
            x[i] = s / M[i][i]

        return x
