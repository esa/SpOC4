import pykep as pk
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


MU_MOON = 4.904869500000000e12
R_MOON = 1737.4e3

class TomatoProblem:
        
    """Utility class for generating tomato orbits and computing transfer cost."""

    mu_moon = MU_MOON
    r_moon = R_MOON

    def __init__(
            self, 
            file_path: str = None,
            max_revs: int = 20):
        if file_path is not None:
            self._load_problem(file_path)

        self.max_revs = max_revs



    def _load_problem(self, file_path: str):
        """Load a KTTSP instance from disk and initialize all dependent attributes."""
        file_path = self._get_file_path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"KTTSP instance file not found: {file_path}")

        header = None
        orbital_rows = []
        with file_path.open("r", encoding="utf-8") as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith("c"):
                    continue
                if line.startswith("p "):
                    parts = line.split()
                    if len(parts) != 6 or parts[1].lower() != "kttsp":
                        raise ValueError(f"Malformed KTTSP header line: {line}")
                    header = parts
                    continue

                parts = line.split()
                if len(parts) != 6:
                    raise ValueError(f"Expected 6 orbital parameters per tomato, got: {line}")
                orbital_rows.append([float(value) for value in parts])

        if header is None:
            raise ValueError(f"Missing 'p kttsp' header in file: {file_path}")
        if not orbital_rows:
            raise ValueError(f"No orbital parameter rows found in file: {file_path}")

        # Header format: p kttsp t0[d] min_tof[s] max_total_time[d] dv_threshold[m/s]
        self.t0 = pk.epoch(float(header[2])) 
        self.min_tof = float(header[3]) / pk.DAY2SEC
        self.max_time = float(header[4])
        self.dv_threshold = float(header[5])

        self.tomato_orbital_parameters = np.asarray(orbital_rows, dtype=float)
        self.n_tomatoes = self.tomato_orbital_parameters.shape[0]

        self.tomatoes = []
        for row in self.tomato_orbital_parameters:
            a, e, incl, raan, arg_periapsis, true_anomaly = row
            elements = [a, e, incl, raan, arg_periapsis, true_anomaly]
            tomato = pk.planet.keplerian(self.t0, elements, MU_MOON, 0.0, 0.0, 0.0)
            self.tomatoes.append(tomato)

    def _get_file_path(
        self,
        path: Path,
    ) -> Path:
        """
        A convenience function that searches for a graph file
        in a hierarchy of directories.

        Args:
            path (Path):
                The (relative) file path.

        Raises:
            FileNotFoundError:
                Raises an error if the file does not exist.

        Returns:
            Path:
                A filesystem path.
        """

        fpath = str(path)
        root = Path(__file__)
        for _ in range(len(root.parents)):
            root = root.parent.resolve().absolute()
            path = root / fpath

            if path.exists():
                break

        if not path.exists():
            raise FileNotFoundError(f"File '{path}' does not exist.")

        return Path(path)

    def compute_transfer(self, i_from: int, i_to: int, t_start: float, tof: float):
        """
        Compute total DeltaV for transfer between two tomatoes.

        Parameters:
            i_from (int): Index of departure tomato.
            i_to (int): Index of arrival tomato.
            t_start (float): Start time in days.
            tof (float): Time of flight in days.
        Returns:
            float: Total DeltaV required for the transfer.
        """
        tof_sec = tof * pk.DAY2SEC

        r_i, v_i = self.tomatoes[i_from].eph(t_start)
        r_j, v_j = self.tomatoes[i_to].eph(t_start + tof)

        best_dv = float('inf')
        for cw in [False, True]:
            try:
                lambert_solution = pk.lambert_problem(r_i, r_j, tof_sec, self.mu_moon, cw, self.max_revs)
            except Exception:
                continue
            for v1, v2 in zip(lambert_solution.get_v1(), lambert_solution.get_v2()):
                delta_v_i = np.linalg.norm(np.array(v1) - np.array(v_i))
                delta_v_j = np.linalg.norm(np.array(v2) - np.array(v_j))
                total_dv = delta_v_i + delta_v_j
                if total_dv < best_dv:
                    best_dv = total_dv
        return best_dv
    
    def find_transfer(self, i_from: int, i_to: int, t_start: float, 
                      dv_threshold: float, max_time: float = 5.0, n_steps: int = 1000):
        """
        Find the earliest valid transfer time between two tomatoes that meets the DeltaV threshold.

        Parameters:
            i_from (int): Index of departure tomato.
            i_to (int): Index of arrival tomato.
            t_start (float): Start time in days.
            tof (float): Time of flight in days.
            dv_threshold (float): Maximum allowed DeltaV for the transfer.
            max_time (float, optional): Maximum time to search for a valid transfer in days. 
                Default is 5.0 days.
            n_steps (int, optional): Number of time steps to evaluate within the max_time window. 
                Default is 1000.
        Returns:
            float: Time of flight for the earliest valid transfer.
        """
        time_grid = np.linspace(0, max_time, n_steps)
        for tof in time_grid[1:]:
            if self.compute_transfer(i_from, i_to, t_start, tof) < dv_threshold:
                return tof
        raise ValueError(f"Could not find valid transfer from tomato {i_from} to {i_to} within dv threshold.")
        
    def plot_orbits(self, indexes=None, t_start=None, t_end=None, n_points=200):
        """
        Plot the orbits of the tomatoes.

        Parameters:
        indexes (list, optional): List of tomato indexes to plot. If None, plots all tomatoes.
        t_start (pk.epoch or float, optional): Start epoch of the time window to plot.
            If None, defaults to self.t0. Only used when t_end is also provided.
        t_end (pk.epoch or float, optional): End epoch of the time window to plot.
            If None, the full orbit is plotted using plot_planet.
            When provided, the arc from t_start to t_end is plotted and the dot
            marks the position at t_end.
        n_points (int, optional): Number of points used to draw the arc when a time
            window is specified. Default is 200.
        """
        if not self.tomatoes:
            raise ValueError("No tomatoes available. Call create_tomatoes() before plot_orbits().")

        if indexes is None:
            indexes = range(len(self.tomatoes))

        indexes = list(indexes)
        if not indexes:
            raise ValueError("indexes must contain at least one tomato index.")

        max_index = len(self.tomatoes) - 1
        for idx in indexes:
            if not isinstance(idx, (int, np.integer)):
                raise TypeError("All indexes must be integers.")
            if idx < 0 or idx > max_index:
                raise IndexError(f"Tomato index {idx} is out of bounds. Valid range is [0, {max_index}].")

        use_window = t_end is not None
        if use_window:
            t_s = self.t0 if t_start is None else t_start
            t_s_mjd = t_s.mjd2000 if isinstance(t_s, pk.epoch) else float(t_s)
            t_e_mjd = t_end.mjd2000 if isinstance(t_end, pk.epoch) else float(t_end)
            if t_e_mjd <= t_s_mjd:
                raise ValueError("t_end must be greater than t_start.")
            times = np.linspace(t_s_mjd, t_e_mjd, n_points)

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")

        ax.scatter(0.0, 0.0, 0.0, c="gray", s=80, label="Moon")

        colors = plt.cm.tab20(np.linspace(0.0, 1.0, len(indexes)))
        for color, idx in zip(colors, indexes):
            tomato = self.tomatoes[idx]
            if use_window:
                positions = [tomato.eph(pk.epoch(t)) for t in times]
                xs = [r[0] for r, _ in positions]
                ys = [r[1] for r, _ in positions]
                zs = [r[2] for r, _ in positions]
                ax.plot(xs, ys, zs, color=color)
                r_end, _ = tomato.eph(pk.epoch(t_e_mjd))
                ax.scatter(r_end[0], r_end[1], r_end[2], c=[color], s=20)
            else:
                try:
                    pk.orbit_plots.plot_planet(
                        tomato,
                        self.t0,
                        1.0,
                        axes=ax,
                        legend=False,
                        color=color,
                    )
                except TypeError:
                    pk.orbit_plots.plot_planet(
                        tomato,
                        self.t0,
                        1.0,
                        ax=ax,
                        legend=False,
                        color=color,
                    )
                r, _ = tomato.eph(self.t0)
                ax.scatter(r[0], r[1], r[2], c=[color], s=20)

        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_title("Tomato Orbits Around the Moon")
        ax.legend(loc="best")
        plt.tight_layout()
        return fig, ax
    
    def plot_transfer(self, i_from, i_to, t_start, tof, n_points=200):
        """
        Plot the transfer trajectory between two tomatoes.

        Parameters:
            i_from (int): Index of departure tomato.
            i_to (int): Index of arrival tomato.
            t_start (float): Start time in seconds.
            tof (float): Time of flight in seconds.
            n_points (int, optional): Number of points used to draw the transfer trajectory. Default is 200.
        """

        r_i, v_i = self.tomatoes[i_from].eph(t_start)
        r_j, v_j = self.tomatoes[i_to].eph(t_start + tof)
        if not self.tomatoes:
            raise ValueError("No tomatoes available. Call create_tomatoes() before plot_transfer().")

        r_i, v_i = self.tomatoes[i_from].eph(t_start)
        r_j, v_j = self.tomatoes[i_to].eph(t_start + tof)

        tof_sec = tof * pk.DAY2SEC

        # Solve Lambert for both directions, up to 20 revolutions, and pick the best
        best_dv = float('inf')
        best_lambert = None
        best_sol_idx = 0
        for cw in [False, True]:
            for max_revs in range(self.max_revs+1):
                try:
                    lambert_solution = pk.lambert_problem(r_i, r_j, tof_sec, self.mu_moon, cw, max_revs)
                except Exception:
                    continue
                for idx, (v1, v2) in enumerate(zip(lambert_solution.get_v1(), lambert_solution.get_v2())):
                    delta_v_i = np.linalg.norm(np.array(v1) - np.array(v_i))
                    delta_v_j = np.linalg.norm(np.array(v2) - np.array(v_j))
                    total_dv = delta_v_i + delta_v_j
                    if total_dv < best_dv:
                        best_dv = total_dv
                        best_lambert = lambert_solution
                        best_sol_idx = idx

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(0.0, 0.0, 0.0, c="gray", s=80, label="Moon")

        # Draw the Lambert transfer arc using pykep's built-in plotter
        pk.orbit_plots.plot_lambert(best_lambert, N=n_points, axes=ax, color="blue", legend=False)

        ax.scatter(r_i[0], r_i[1], r_i[2], c="green", s=20, label="Departure Tomato")
        ax.scatter(r_j[0], r_j[1], r_j[2], c="red", s=20, label="Arrival Tomato")
        ax.set_xlabel("x [m]")  
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_title("Transfer Trajectory Between Tomatoes")
        ax.legend(loc="best")
        plt.tight_layout()
        return fig, ax
    
    def plot_full_trajectory(self, x, n_points=200):
            """
            Plot the full trajectory (all transfers and waiting arcs) for a given decision vector x.

            Parameters:
                x (array-like): Decision vector [t1,...,tN-1, tof1,...,tofN-1, pi1,...,piN]
                    t_i: departure time of i-th transfer (days)
                    tof_i: time of flight of i-th transfer (days)
                    pi_i: index of the tomato at step i (integer)
                n_points (int, optional): Number of points for each arc/transfer. Default is 200.
            Returns:
                fig, ax: Matplotlib figure and axes objects.
            """
            import matplotlib.pyplot as plt
            import numpy as np
            import pykep as pk

            N = (len(x) // 3) + 1  # Number of tomatoes visited
            n_legs = N - 1
            t_departures = np.array(x[:n_legs])
            tofs = np.array(x[n_legs:2*n_legs])
            pis = np.array(x[2*n_legs:]).astype(int)

            if not hasattr(self, 'tomatoes') or not self.tomatoes:
                raise ValueError("No tomatoes available. Call create_tomatoes() before plot_full_trajectory().")

            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(111, projection="3d")
            ax.scatter(0.0, 0.0, 0.0, c="gray", s=80, label="Moon")

            # Plot the full trajectory
            for k in range(n_legs):
                i_from = pis[k]
                i_to = pis[k+1]
                t_start = t_departures[k]
                tof = tofs[k]

                # Plot waiting arc on current tomato (from previous arrival to departure)
                if k == 0:
                    # First leg: wait from t0 to t_departures[0]
                    t_wait_start = self.t0.mjd2000 if hasattr(self, 't0') else 0.0
                else:
                    t_wait_start = t_departures[k-1] + tofs[k-1]
                t_wait_end = t_start
                if t_wait_end > t_wait_start:
                    times = np.linspace(float(t_wait_start), float(t_wait_end), n_points)
                    tomato = self.tomatoes[i_from]
                    positions = [tomato.eph(pk.epoch(float(t))) for t in times]
                    xs = [r[0] for r, _ in positions]
                    ys = [r[1] for r, _ in positions]
                    zs = [r[2] for r, _ in positions]
                    ax.plot(xs, ys, zs, color="orange", linestyle="dotted", label="Waiting" if k == 0 else None)

                # Plot transfer arc using Lambert
                r_i, v_i = self.tomatoes[i_from].eph(pk.epoch(float(t_start)))
                r_j, v_j = self.tomatoes[i_to].eph(pk.epoch(float(t_start + tof)))
                tof_sec = tof * pk.DAY2SEC
                best_dv = float('inf')
                best_lambert = None
                for cw in [False, True]:
                    for max_revs in range(self.max_revs+1):
                        try:
                            lambert_solution = pk.lambert_problem(r_i, r_j, tof_sec, self.mu_moon, cw, max_revs)
                        except Exception:
                            continue
                        for v1, v2 in zip(lambert_solution.get_v1(), lambert_solution.get_v2()):
                            delta_v_i = np.linalg.norm(np.array(v1) - np.array(v_i))
                            delta_v_j = np.linalg.norm(np.array(v2) - np.array(v_j))
                            total_dv = delta_v_i + delta_v_j
                            if total_dv < best_dv:
                                best_dv = total_dv
                                best_lambert = lambert_solution
                if best_lambert is not None:
                    pk.orbit_plots.plot_lambert(best_lambert, N=n_points, axes=ax, color="blue", legend=False)

                # Mark departure and arrival
                ax.scatter(r_i[0], r_i[1], r_i[2], c="green", s=20, label="Departure" if k == 0 else None)
                ax.scatter(r_j[0], r_j[1], r_j[2], c="red", s=20, label="Arrival" if k == n_legs-1 else None)

            ax.set_xlabel("x [m]")
            ax.set_ylabel("y [m]")
            ax.set_zlabel("z [m]")
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys(), loc="best")
            plt.tight_layout()
            return fig, ax
