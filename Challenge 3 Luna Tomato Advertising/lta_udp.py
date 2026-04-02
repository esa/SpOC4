import re
from pathlib import Path
import copy

import numpy as np
import heyoka as hy
import matplotlib.pyplot as plt

class celestial_morse_code():
    """
        Celestial Morse Code (UDP)

        For more information about the competition, please refer to the competition website:

        https://www.esa.int/gsp/ACT/news/spoc-2026/

        This class conforms to the pygmo UDP format.

    """

    # Initialisation
    def __init__(
            self,
            message=" Tomatoes for sale ",
        ):
        """        
        Initialize the UDP instance.

        Loads orbital databases and pre-processes their frequencies and amplitudes.
        Generates the target Morse signal for the given message.
        
        Args:
            message (str): Message to encode in Morse (default: "Tomatoes for sale")
        """
        # Extract Inputs
        self.message = message

        # Taylor adaptive integrator for CR3BP
        self.ta, self.TU, self.LU = self._create_ta()
        self.synodic_period = (self.TU / 24 / 3600 * 2*np.pi) # Days

        # Scaling morse code paramters to length of message
        dot, dash, gap, letter_gap, word_gap = 2, 5, 2, 5, 7
        message_length = self.morse_duration(self.message, dot=dot, dash=dash, gap=gap, letter_gap=letter_gap, word_gap=word_gap)
        scaling_factor = (0.9 * self.synodic_period) / message_length
        self.dot, self.dash, self.gap, self.letter_gap, self.word_gap = [x * scaling_factor for x in (dot, dash, gap, letter_gap, word_gap)]
        
        # occulatation window (box width)
        self.width = self.dot / 2 # Days
        self.dt = self.dot / 20 # Days

        # Extract signal
        self.t, self.target_signal = self.morse_to_signal(message)
        self.morse_string = self.message_to_morse_string(message)

        # Database
        self.db_dro = np.loadtxt(self._get_file_path("data/spoc4/cmc/db_dro.txt"), delimiter=",", skiprows=1)
        self.db_lyap = np.loadtxt(self._get_file_path("data/spoc4/cmc/db_lyap.txt"), delimiter=",", skiprows=1)
        self.db_axial = np.loadtxt(self._get_file_path("data/spoc4/cmc/db_axial.txt"), delimiter=",", skiprows=1)

        # Amplitude
        # DB columns: 0-5=state, 6=period, 7=freq1, [8=freq2], 8or9=amp1, [9or10=amp2]
        # DRO (1 occultation): col7=freq, col8=amp
        # Lyap/Axial (2 occultations): col7=freq1, col8=freq2, col9=amp1, col10=amp2
        # NB (simpler problem, all amplitudes = 1): We ignore the amplitude values from the DB and set all to 1, since the optimization can scale the number of spacecraft on each orbit to achieve the desired amplitude modulation.
        self.amplitudes_dro = self.db_dro[:,8] * 0 + 1
        self.amplitudes_lyap = np.concatenate([
            self.db_lyap[:,9],
            self.db_lyap[:,10]]
        ) * 0 + 1
        self.amplitudes_axial = np.concatenate([
            self.db_axial[:,9],
            self.db_axial[:,10]]
        ) * 0 + 1
        # ... and frequency
        self.frequencies_dro = self.db_dro[:,7]
        self.frequencies_lyap = np.concatenate([
            self.db_lyap[:,7],
            self.db_lyap[:,8]]
        )
        self.frequencies_axial = np.concatenate([
            self.db_axial[:,7],
            self.db_axial[:,8]]
        )

        # Number of orbits in each family
        self.num_dro = len(self.db_dro)
        self.num_lyap = len(self.db_lyap)
        self.num_axial = len(self.db_axial)
        self.num_orbits = self.num_dro + self.num_lyap + self.num_axial

        # Occulations per family
        self.occult_dro = 1
        self.occult_lyap = 2
        self.occult_axial = 2

        # Max number of spacecraft per orbit
        self.max_per_orbit = 5

        # Threshold for MSE constraint
        self.mse_thresh = 0.05

    def get_name(self):
        """
        Returns the name of the User-Defined Problem (UDP).
        
        Returns:
            str: Name of the UDP
        """
        return "Celestial Morse Code"

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
    
    # Fitness
    def fitness(self, x, postprocess=False):
        """
        Evaluates the fitness of a chromosome.
        
        Decodes the chromosome into orbit selections and phases, builds the occultation
        basis matrix, reconstructs the target morse signal, clips it (np.clip(reconstructed, 0, 1), 
        and computes a two-tiered objective: first minimizing MSE if above threshold, then minimizing number of
        orbits used if MSE is acceptable.
        
        Args:
            x (list or array): Chromosome with structure:
                - First N entries: orbit selection counts (0-5) for each orbit family
                - Remaining N*5 entries: phase values for each orbit/spacecraft combination
            postprocess (bool): If True, also return number of orbits and MSE value
        
        Returns:
            float: Objective value to minimize
            tuple (if postprocess=True): (objective, num_orbits_used, mse)
        
        Structure of x:
        x = [
                dro_orbit_1, dro_orbit_2, ..., dro_orbit_m,
                lyap_orbit_1, lyap_orbit_2, ..., lyap_orbit_n,
                axial_orbit_1, axial_orbit_2, ..., axial_orbit_p,
                dro_phase_1a, dro_phase_1b, dro_phase_1c, dro_phase_1d, dro_phase_1e, 
                    dro_phase_2a, dro_phase_2b, dro_phase_2c, dro_phase_2d, dro_phase_2e, 
                    ...
                    dro_phase_ma, dro_phase_mb, dro_phase_mc, dro_phase_md, dro_phase_me,
                lyap_phase_1a, lyap_phase_1b, lyap_phase_1c, lyap_phase_1d, lyap_phase_1e, 
                    lyap_phase_2a, lyap_phase_2b, lyap_phase_2c, lyap_phase_2d, lyap_phase_2e, 
                    ...
                    lyap_phase_ma, lyap_phase_mb, lyap_phase_mc, lyap_phase_md, lyap_phase_me,
                axial_phase_1a, axial_phase_1b, axial_phase_1c, axial_phase_1d, axial_phase_1e, 
                    axial_phase_2a, axial_phase_2b, axial_phase_2c, axial_phase_2d, axial_phase_2e, 
                    ...
                    axial_phase_ma, axial_phase_mb, axial_phase_mc, axial_phase_md, axial_phase_me
            ]
        """

        # decode x into orbits and phases
        selections, phases = self._decode_chromosome(x)

        # Reconstruct Signal
        A_matrix = self._build_occultation_matrix(selections, phases)
        reconstructed = A_matrix @ np.ones(A_matrix.shape[1])

        # clip
        reconstructed = np.clip(reconstructed, 0, 1)
        # reconstructed = np.where(reconstructed > 0.5, 1, 0) # hard clip to binary signal

        # Compute MSE
        mse = np.mean((reconstructed - self.target_signal)**2)

        # Two-tiered objective:
        num_selected = np.sum(selections)
        if mse > self.mse_thresh:
            # 1) If MSE > threshold, minimize mse
            max_num_selected = self.max_per_orbit*(self.num_orbits)
            obj = max_num_selected + 1e3*(np.sqrt(mse) - np.sqrt(self.mse_thresh))
        else:
            # 2) If MSE <= threshold, minimize number of orbits used
            obj = num_selected

        if not postprocess:
            return [obj]
        else:
            return [obj], num_selected, mse
    
    # Bounds and Constraints
    def get_bounds(self):
        """
        Returns the lower and upper bounds for the chromosome variables.
        
        Returns:
            tuple: (lower_bounds, upper_bounds)
                - Orbit selections: [0, max_per_orbit]
                - Phase values: [0, 2π]
        """
        lb = (
            [0] * self.num_dro 
            + [0] * self.num_lyap 
            + [0] * self.num_axial
            + [0] * (self.num_dro * self.max_per_orbit)
            + [0] * (self.num_lyap * self.max_per_orbit)
            + [0] * (self.num_axial * self.max_per_orbit)
        )
        ub = (
            [self.max_per_orbit] * self.num_dro 
            + [self.max_per_orbit] * self.num_lyap 
            + [self.max_per_orbit] * self.num_axial
            + [2*np.pi] * (self.num_dro * self.max_per_orbit)
            + [2*np.pi] * (self.num_lyap * self.max_per_orbit)
            + [2*np.pi] * (self.num_axial * self.max_per_orbit)
        )      
        return (lb, ub)
    
    def get_nobj(self):
        """
        Our objective is to minimize is number of orbits used
        """
        return 1

    def get_nic(self):
        """
        Returns the number of inequality constraints.
        
        Returns:
            0: No inequality constraints defined
        """
        return 0
    
    def get_nec(self):
        """
        Returns the number of equality constraints.
        
        Returns:
            0: No equality constraints defined
        """
        return 0

    def _create_ta(self):
        """
        Creates a Heyoka Taylor adaptive integrator for the Circular Restricted Three Body Problem (CR3BP).

        Returns:
            tuple: (ta, TU, LU)
                ta: configured heyoka.taylor_adaptive integrator
                TU (float): time unit in seconds
                LU (float): length unit in meters
        """
        # Parameters
        mu_earth = 3.986004418e14 # m^3 / s^2 
        mu_moon = 4.9048695e12 # m^3 / s^2

        LU = 389703e3 # m
        n = np.sqrt((mu_earth+mu_moon)/LU**3)
        TU = 1/n

        r_earth = 6378.137e3 / LU
        r_moon = 1737.4e3 / LU

        # Other parameters
        mu = mu_moon / (mu_earth + mu_moon)

        # The state
        x, y, z, vx, vy, vz, m = hy.make_vars("x", "y", "z", "vx", "vy", "vz", "m")
        # Parameters
        _mu = hy.par[0]

        # Useful expressions
        rP1_3 = ((x-(-_mu))**2 + y**2 + z**2) ** (1.5)
        rP2_3 = ((x-(1-_mu))**2 + y**2 + z**2) ** (1.5)

        # Vectors for convenience of math manipulation

        r = np.array([x, y, z])
        rP1 = np.array([-_mu, 0, 0])
        rP2 = np.array([1-_mu, 0, 0])
        r2 = np.array([x, y, z])
        v = np.array([vx, vy, vz])
        omega = np.array([0, 0, 1])

        def hy_cross3(a, b):
            # Cross product of 3D vectors
            return np.array([
                a[1]*b[2] - a[2]*b[1],
                a[2]*b[0] - a[0]*b[2],
                a[0]*b[1] - a[1]*b[0]
            ])

        # Dynamics
        fr = v
        fv = -(1-_mu)*(r-rP1)/(rP1_3) + (-_mu)*(r-rP2)/(rP2_3) - 2 * hy_cross3(omega, v) + hy_cross3(hy_cross3(omega, r), omega)

        rhs = list(fr) + list(fv)
        full_state = [x, y, z, vx, vy, vz]

        # We assemble the Taylor adaptive integrator
        dyn = [(var, dvar) for var, dvar in zip(full_state, rhs)]

        # Event
        occult_times = []
        occult_states = []
        def t_cb(ta, d_sgn): 
            x, _, _, _, _, _ = ta.state
            if x < 1-mu and x > -mu:
                occult_times.append(copy.deepcopy(ta.time))
                occult_states.append(copy.deepcopy(ta.state))
            return True

        # Occulatation: If y and z position less than radius of the Moon
        t_ev = hy.t_event(
            hy.sqrt(y**2 + z**2) - r_moon,
            callback=t_cb,
            cooldown=float(1e-6),
            direction=hy.event_direction.positive,
        )

        ta = hy.taylor_adaptive(
            dyn, 
            state=[1.0] * 6, 
            tol=1e-18, 
            t_events = [t_ev],
            high_accuracy=True, 
            compact_mode=True
        )
        ta.pars[0] = mu

        return ta, TU, LU

    def _morse_dict(self):
        # Morse code dictionary mapping characters to their Morse code representations.
        # Includes letters A-Z, digits 0-9, common punctuation, and space (as word separator).
        return {
            # Letters
            'A': '.-',     'B': '-...',   'C': '-.-.',   'D': '-..',
            'E': '.',      'F': '..-.',   'G': '--.',    'H': '....',
            'I': '..',     'J': '.---',   'K': '-.-',    'L': '.-..',
            'M': '--',     'N': '-.',     'O': '---',    'P': '.--.',
            'Q': '--.-',   'R': '.-.',    'S': '...',
            'T': '-',      'U': '..-',    'V': '...-',
            'W': '.--',    'X': '-..-',   'Y': '-.--',   'Z': '--..',

            # Numbers
            '0': '-----',  '1': '.----',  '2': '..---',
            '3': '...--',  '4': '....-',  '5': '.....',
            '6': '-....',  '7': '--...',  '8': '---..',
            '9': '----.',

            # Punctuation
            '.': '.-.-.-',
            ',': '--..--',
            '?': '..--..',
            "'": '.----.',
            '!': '-.-.--',
            '/': '-..-.',
            '(': '-.--.',
            ')': '-.--.-',
            '&': '.-...',
            ':': '---...',
            ';': '-.-.-.',
            '=': '-...-',
            '+': '.-.-.',
            '-': '-....-',
            '_': '..--.-',
            '"': '.-..-.',
            '$': '...-..-',
            '@': '.--.-.',

            # Space (word separator)
            ' ': '/'
        }

    def _decode_chromosome(self, x):
        """
        Decode chromosome into selection counts and phase vector.

        Args:
            x (Sequence[float]): Full chromosome.

        Returns:
            tuple[list[int], Sequence[float]]: (selections, phases)
        """
        selections = list(np.maximum(np.round(x[:self.num_orbits]).astype(int), 0))
        phases = x[self.num_orbits:]
        return selections, phases

    # Occultation Signal Construction
    def _heaviside_relu(self, x, k=1e6):
        """
        Soft approximation of the Heaviside step function.
        
        Implements a smooth transition from 0 to 1 with steepness controlled by k.
        Used to create sharp but differentiable window functions.
        
        Args:
            x (array): Input values
            k (float): Steepness parameter (default: 1e6, very steep)
        
        Returns:
            array: Values clipped to [0, 1] range
        """
        return np.minimum(1, np.maximum(0, k*x))

    def _repeating_box_window(self, t, center, frequency):
        """
        Creates a repeating box window signal at specified frequency.
        
        Generates a time-periodic window function that activates (=1) within +/- width/2
        around the phase center, and repeats at the given frequency. Used to modulate
        the occultation signals.
        
        Args:
            t (array): Time points
            center (float): Phase center (in time units)
            frequency (float): Repetition frequency (Hz)
        Returns:
            array: Window signal values [0, 1]
        """
        period = 1.0 / frequency

        # Convert time (days) to time (nd)
        t_nd = (t * 24 * 3600) / self.TU
        width_nd = (self.width * 24 * 3600) / self.TU

        # time relative to phase center within each period
        phase_t = (t_nd - center + period/2) % period - period/2

        left  = self._heaviside_relu(phase_t + width_nd/2)
        right = self._heaviside_relu(phase_t - width_nd/2)

        return left - right
    
    def _build_occultation_matrix(self, selections, phases):
        """
        Constructs the occultation basis matrix from orbit selections and phases.
        
        Creates a matrix where each column represents the occultation signal from one
        spacecraft on one orbit. Different orbit types have different numbers of occultations:
        - DRO: 1 occultation per orbit
        - Lyapunov: 2 occultations per orbit (π phase shift for second)
        - Axial: 2 occultations per orbit (π phase shift for second)
        
        Each column is an amplitude-modulated repeating box window at the orbit's frequency.
        
        Args:
            selections (array): Number of spacecraft per orbit (0-5)
            phases (array): Phase values for each spacecraft
        
        Returns:
            array: Basis matrix of shape (len(t), num_basis_functions)
        """
        # Create basis matrix
        A_matrix = np.zeros((len(self.t), self.max_per_orbit*(self.num_dro*self.occult_dro + self.num_lyap*self.occult_lyap + self.num_axial*self.occult_axial)))
        ai = 0 # counter

        # Construct DRO
        for orbit_idx in range(0, self.num_dro):
            if selections[orbit_idx] > 0:
                for j in range(selections[orbit_idx]):
                    ph = phases[orbit_idx * self.max_per_orbit + j]
                    A_matrix[:, ai] = self.amplitudes_dro[orbit_idx] * self._repeating_box_window(self.t, (2*np.pi - ph)/(2*np.pi)/self.frequencies_dro[orbit_idx], frequency=self.frequencies_dro[orbit_idx])
                    ai += 1
        # Construct Lyap
        for orbit_idx in range(0, self.num_lyap):
            sel = selections[self.num_dro + orbit_idx]
            if sel > 0:
                for j in range(sel):
                    ph = phases[self.num_dro * self.max_per_orbit + orbit_idx * self.max_per_orbit + j]
                    A_matrix[:, ai] = self.amplitudes_lyap[orbit_idx] * self._repeating_box_window(self.t, (2*np.pi - ph)/(2*np.pi)/self.frequencies_lyap[orbit_idx], frequency=self.frequencies_lyap[orbit_idx])
                    A_matrix[:, ai+1] = self.amplitudes_lyap[self.num_lyap + orbit_idx] * self._repeating_box_window(self.t, (ph + np.pi)/self.frequencies_lyap[self.num_lyap + orbit_idx], frequency=self.frequencies_lyap[self.num_lyap + orbit_idx])
                    ai += 2
        # Construct Axial
        for orbit_idx in range(0, self.num_axial):
            sel = selections[self.num_dro + self.num_lyap + orbit_idx]
            if sel > 0:
                for j in range(sel):
                    ph = phases[(self.num_dro + self.num_lyap) * self.max_per_orbit + orbit_idx * self.max_per_orbit + j]
                    A_matrix[:, ai] = self.amplitudes_axial[orbit_idx] * self._repeating_box_window(self.t, (2*np.pi - ph)/(2*np.pi)/self.frequencies_axial[orbit_idx], frequency=self.frequencies_axial[orbit_idx])
                    A_matrix[:, ai+1] = self.amplitudes_axial[self.num_axial + orbit_idx] * self._repeating_box_window(self.t, (ph + np.pi)/self.frequencies_axial[self.num_axial + orbit_idx], frequency=self.frequencies_axial[self.num_axial + orbit_idx])
                    ai += 2
        return A_matrix

    # Morse Code Functions
    def morse_duration(self, message, dot=2.0, dash=5.0, gap=2.0, letter_gap=5.0, word_gap=7.0):
        """
        Compute total Morse timing length for a message.

        Args:
            message (str): Plain-text message.
            dot (float): Dot duration (timing units).
            dash (float): Dash duration (timing units).
            gap (float): Intra-character gap (timing units).
            letter_gap (float): Inter-letter gap (timing units).
            word_gap (float): Inter-word gap (timing units).

        Returns:
            float: Total duration in the same timing units.
        """

        morse = self.message_to_morse_string(message)

        # Normalize word separators to "/"
        s = re.sub(r"\s*/\s*", " / ", morse.strip())
        s = re.sub(r"\s{3,}", " / ", s)

        words = [w.strip() for w in s.split("/") if w.strip()]
        total = 0

        for wi, word in enumerate(words):
            letters = [l for l in word.split(" ") if l]

            for li, letter in enumerate(letters):
                for si, sym in enumerate(letter):
                    if sym == ".":
                        total += dot
                    elif sym == "-":
                        total += dash
                    else:
                        continue  # ignore unknown chars

                    if si < len(letter) - 1:
                        total += gap  # between symbols in same letter

                if li < len(letters) - 1:
                    total += letter_gap  # between letters

            if wi < len(words) - 1:
                total += word_gap  # between words

        return total

    def morse_to_signal(self, message):
        """
        Convert a message to a binary time-series signal.

        Args:
            message (str): Plain-text message.

        Returns:
            tuple[np.ndarray, np.ndarray]:
                - time grid in days
                - binary target signal (0/1)
        """
        t = 0.0
        segments = []  # list of (t_start, t_end) for "on" regions

        for char in message.upper():
            if char == ' ':
                t += self.word_gap
                continue

            code = self._morse_dict()[char]
            for symbol in code:
                duration = self.dot if symbol == '.' else self.dash
                segments.append((t, t + duration))
                t += duration
                t += self.gap  # intra-symbol gap

            # letter gap (already added one gap above, so subtract it)
            t += (self.letter_gap - self.gap)

        total_time = t
        t_full = np.arange(0, total_time, self.dt)
        s_full = np.zeros_like(t_full)

        # Fill contiguous blocks
        for t_start, t_end in segments:
            idx_start = int(round(t_start / self.dt))
            idx_end   = int(round(t_end   / self.dt))
            s_full[idx_start:idx_end] = 1.0

        return t_full, s_full

    def message_to_morse_string(self, message):
        """
        Converts a text message to a Morse code string representation.
        
        Args:
            message (str): Message to convert
            self._morse_dict (dict): Morse code dictionary
        
        Returns:
            str: Space-separated Morse code symbols (e.g., '.- ... ---')
        """
        return ' '.join(self._morse_dict()[c] for c in message.upper())
    
    # Plotting Functions
    def plot_target(self, message=None):
        """
        Plots the target Morse code signal for a given message.
        
        Args:
            message (str): Message to plot. If None, uses self.message
        """
        if message is None:
            message = self.message
        t_signal, target_signal = self.morse_to_signal(message)
        morse_string = self.message_to_morse_string(message)

        # Create stacked subplots
        fig, ax1 = plt.subplots(1, 1, figsize=(12, 5), sharex=True)

        # Target signal
        ax1.plot(t_signal, target_signal)
        ax1.set_ylabel("Amplitude")
        ax1.grid(True)

        ax1.set_ylim([0, 1.1])

        ax1.set_title(f"Target Morse Signal\n{message}\n{morse_string}", fontsize=14)

    def plot_signal(self, x):
        """
        Plots both the target and reconstructed Morse code signals.
        
        Visualizes the target signal, the reconstructed signal from orbital basis,
        and displays the MSE and number of orbits used in the solution.
        
        Args:
            x (array): Solution chromosome to visualize
        """
        message = self.message
        t_signal, target_signal = self.morse_to_signal(message)
        morse_string = self.message_to_morse_string(message)

        # decode x into orbits and phases
        selections, phases = self._decode_chromosome(x)

        # Reconstruct Signal
        A_matrix = self._build_occultation_matrix(selections, phases)
        reconstructed = A_matrix @ np.ones(A_matrix.shape[1])

        # Clip
        reconstructed = np.clip(reconstructed, 0, 1)
        # reconstructed = np.where(reconstructed > 0.5, 1, 0) # hard clip to binary signal

        # # Normalise
        # if np.max(reconstructed) > 0:
        #     reconstructed = reconstructed / np.max(reconstructed)

        # Compute MSE
        mse = np.mean((reconstructed - self.target_signal)**2)
        num_orbits = np.sum(selections)

        # Create stacked subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

        # Target signal
        ax1.plot(t_signal, target_signal)
        ax1.set_ylabel("Amplitude")
        ax1.set_title("Target Morse Continuous Signal")
        ax1.grid(True)

        ax1.set_ylim([0, 1.1])

        # Reconstructed signal
        ax2.plot(t_signal, reconstructed)
        ax2.set_xlabel("Time (days)")
        ax2.set_ylabel("Amplitude")
        ax2.set_title(f"Reconstructed Signal from Orbital Basis \n #{num_orbits} MSE {mse:.4g}")
        ax2.grid(True)

        ax2.set_ylim([0, 1.1])

        # ax2.set_xlim(0, 10)

        # Overall title
        fig.suptitle(f"{message}\n{morse_string}", fontsize=14)

        plt.tight_layout(rect=[0, 0, 1, 0.94])
        plt.show()

    def propagate_orbit(self, orbit_idx, phase, orbit_type, nstep=1000):
        """
        Propagate a single orbit for visualization.

        Args:
            orbit_idx (int): Orbit index within the selected family database.
            phase (float): Start phase in radians.
            orbit_type (str): One of {"dro", "lyap", "axial"}.
            nstep (int): Number of propagation grid points.

        Returns:
            tuple[np.ndarray, np.ndarray, list, list]:
                (tgrid, integration, occult_times, occult_states)
        """

        if orbit_type == 'dro':
            state = self.db_dro[orbit_idx, 0:6]
            period = self.db_dro[orbit_idx, 6]
        elif orbit_type == 'lyap':
            state = self.db_lyap[orbit_idx, 0:6]
            period = self.db_lyap[orbit_idx, 6]
        elif orbit_type == 'axial':
            state = self.db_axial[orbit_idx, 0:6]
            period = self.db_axial[orbit_idx, 6]
        else:
            raise ValueError("Invalid orbit type")
        
        # Initialise Integrator
        tgrid = np.linspace(0.0, period, nstep) # + phase / 2*np.pi * period
        self.ta.time = 0.0
        self.ta.state[:6] = state

        # Propagate
        occult_times = []
        occult_states = []
        integration = self.ta.propagate_grid(tgrid)[5]

        # Reorder based on phase
        idx_phase = int(np.floor((phase / (2*np.pi)) * nstep ))

        integration_after = integration[idx_phase:, :]
        integration_before = integration[:idx_phase, :]
        tgrid_after = tgrid[idx_phase:]
        tgrid_before = tgrid[:idx_phase]

        # Concatenate to have full orbit starting from phase
        integration = np.concatenate((integration_after, integration_before), axis=0)
        tgrid = np.concatenate((tgrid_after, tgrid_before), axis=0)

        return tgrid, integration, occult_times, occult_states
    
    def plot_orbits(self, x):
        """
        Plots the orbits corresponding to the selected solution.
        
        Visualizes the trajectories of the selected orbits in the CR3BP system.
        
        Args:
            x (array): Solution chromosome to visualize
        """

        # decode x into orbits and phases
        selections, phases = self._decode_chromosome(x)
        
        integration_list = []
        occult_times_list = []
        occult_states_list = []

        print('Selected orbits:', selections)
        print('Selected Phases:', phases)
        # Construct DRO
        for orbit_idx in range(0, self.num_dro):
            if selections[orbit_idx] > 0:
                for j in range(selections[orbit_idx]):
                    ph = phases[orbit_idx * self.max_per_orbit + j]
                    _, tmp_x, tmp_occult_t, tmp_occult_x = self.propagate_orbit(orbit_idx, ph, orbit_type='dro')
                    integration_list.append(tmp_x)
                    occult_times_list.append(tmp_occult_t)
                    occult_states_list.append(tmp_occult_x)

        # Construct Lyap
        for orbit_idx in range(0, self.num_lyap):
            sel = selections[self.num_dro + orbit_idx]
            if sel > 0:
                for j in range(sel):
                    ph = phases[self.num_dro * self.max_per_orbit + orbit_idx * self.max_per_orbit + j]
                    _, tmp_x, tmp_occult_t, tmp_occult_x = self.propagate_orbit(orbit_idx, ph, orbit_type='lyap')
                    integration_list.append(tmp_x)
                    occult_times_list.append(tmp_occult_t)
                    occult_states_list.append(tmp_occult_x)
        # Construct Axial
        for orbit_idx in range(0, self.num_axial):
            sel = selections[self.num_dro + self.num_lyap + orbit_idx]
            if sel > 0:
                for j in range(sel):
                    ph = phases[(self.num_dro + self.num_lyap) * self.max_per_orbit + orbit_idx * self.max_per_orbit + j]
                    _, tmp_x, tmp_occult_t, tmp_occult_x = self.propagate_orbit(orbit_idx, ph, orbit_type='axial')
                    integration_list.append(tmp_x)
                    occult_times_list.append(tmp_occult_t)
                    occult_states_list.append(tmp_occult_x)

        
        # Figure
        fig = plt.figure(figsize=(24, 6))

        ax1 = fig.add_subplot(141, projection='3d')
        ax2 = fig.add_subplot(142)
        ax3 = fig.add_subplot(143)
        ax4 = fig.add_subplot(144)

        plt.style.use('seaborn-v0_8-colorblind')

        for k in range(len(integration_list)):

            traj_x0 = integration_list[k][:, :3]

            ax1.plot(traj_x0[:,0], traj_x0[:,1], traj_x0[:,2])
            ax1.plot(traj_x0[0,0], traj_x0[0,1], traj_x0[0,2], marker='s', color = 'red')
            for ki in range(len(occult_states_list[k])):
                ax1.scatter3D(occult_states_list[k][ki][0], occult_states_list[k][ki][1], occult_states_list[k][ki][2], marker='x', color = 'k')

            ax2.plot(traj_x0[:,0], traj_x0[:,1])
            ax2.plot(traj_x0[0,0], traj_x0[0,1], marker='s', color = 'red')
            for ki in range(len(occult_states_list[k])):
                ax2.scatter(occult_states_list[k][ki][0], occult_states_list[k][ki][1], marker='x', color = 'k')

            ax3.plot(traj_x0[:,0], traj_x0[:,2])
            ax3.plot(traj_x0[0,0], traj_x0[0,2], marker='s', color = 'red')
            for ki in range(len(occult_states_list[k])):
                ax3.scatter(occult_states_list[k][ki][0], occult_states_list[k][ki][2], marker='x', color = 'k')

            ax4.plot(traj_x0[:,1], traj_x0[:,2])
            ax4.plot(traj_x0[0,1], traj_x0[0,2], marker='s', color = 'red')
            for ki in range(len(occult_states_list[k])):
                ax4.scatter(occult_states_list[k][ki][1], occult_states_list[k][ki][2], marker='x', color = 'k')
    
        # Primaries
        ax1.scatter(-self.ta.pars[0], 0, 0, color='green', s=120, label='Earth')
        ax1.scatter(1-self.ta.pars[0], 0, 0, color='grey', s=120, label='Moon')

        ax2.scatter(-self.ta.pars[0], 0, color='green', s=120)
        ax2.scatter(1-self.ta.pars[0], 0, color='gray', s=120)

        ax3.scatter(-self.ta.pars[0], 0, color='green', s=120)
        ax3.scatter(1-self.ta.pars[0], 0, color='gray', s=120)

        ax4.scatter(0, 0, color='green', s=120)
        ax4.scatter(0, 0, color='gray', s=120)

        ax1.set_title('3D')
        ax2.set_title('XY')
        ax3.set_title('XZ')
        ax4.set_title('YZ')

        ax1.set_xlim(-1, 2)
        ax1.set_ylim(-1.5, 1.5)
        ax1.set_zlim(-0.25, 0.25)
        ax2.set_xlim(-1, 2)
        ax2.set_ylim(-1.5, 1.5)
        ax3.set_xlim(-1, 1.5)
        ax3.set_ylim(-0.25, 0.25)
        ax4.set_xlim(-1, 1)
        ax4.set_ylim(-0.25, 0.25)

        for ax in [ax2, ax3, ax4]:
            ax.grid()

        # Global figure title
        fig.suptitle('CR3BP', fontsize=18)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()

udp = celestial_morse_code()