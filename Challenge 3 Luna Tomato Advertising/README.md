<image src="../images/LTA_banner.png" align="center" alt="Challange 1: Torso Decompositions"/>

# Luna Tomato Advertising — Challenge Description

## Introduction

Luna Tomato Logistics is piloting a next-generation advertising concept: using lunar-space infrastructure to broadcast branded Morse-code messages via natural occultations.

As part of this innovation program, you have been tasked with designing a preliminary occultation-based advertising campaign. Each time a spacecraft passes in front of the Moon (from Earth’s perspective), it occults some of the moonlight, which can be interpreted as a Morse **dot** or **dash**.

Your mission is to configure and phase a fleet of spacecraft across available lunar orbit families so that the desired campaign message is reproduced accurately, while using as few assets as possible and maintaining signal quality constraints.

---

## Short Description

Given a target Morse-code signal derived from a secret message, find a set of spacecraft orbits and phases such that:

1. The combined occultation windows reproduce the Morse signal waveform.
2. The number of spacecraft used is minimised.
3. The mean-squared error (MSE) between the reconstructed and target signal remains below a feasibility threshold (`0.05`).

You are provided with orbit databases for three periodic orbit families:

- Distant Retrograde Orbits (**DRO**)
- **Lyapunov** orbits
- **Axial** orbits

Due to logistical constraints, a maximum of **5 spacecraft** can be placed in any single orbit.

To simplify the initial feasibility study of this concept, all spacecraft and orbits produce the same occultation amplitude, and the duration of the occultation is not modelled explicitly. Instead, once a spacecraft passes between the Earth and the Moon, the occultation amplitude and duration are treated as fixed values.

---

## Definition of Provided Data

The UDP (`celestial_morse_code`) exposes the following attributes after instantiation:

| Attribute | Type | Description | Default |
|---|---|---|---|
| `udp.message` | `str` | The secret message to encode | ` Tomatoes for sale ` |
| `udp.synodic_period` | `float` | Time period for the Earth-Moon system (days) | `27.8511` |
| `udp.width` | `float` | Width of one spaceraft occultation window (days) | `0.1036` |
| `udp.dot` | `float` | Duration of a Morse dot in units of `udp.width` | `2.0` |
| `udp.dash` | `float` | Duration of a Morse dash in units of `udp.width` | `5.0` |
| `udp.gap` | `float` | Intra-character gap in units of `udp.width` | `2.0` |
| `udp.letter_gap` | `float` | Inter-letter gap in units of `udp.width` | `5.0` |
| `udp.word_gap` | `float` | Inter-word gap in units of `udp.width` | `7.0` |
| `udp.num_dro` | `int` | Number of DRO orbit candidates | `188` |
| `udp.num_lyap` | `int` | Number of Lyapunov orbit candidates | `291` |
| `udp.num_axial` | `int` | Number of Axial orbit candidates | `300` |
| `udp.max_per_orbit` | `int` | Max occultations schedulable per orbit slot | `5` |
| `udp.db_dro` | `list` | DRO orbit database rows | `loaded (188 rows)` |
| `udp.db_lyap` | `list` | Lyapunov orbit database rows | `loaded (291 rows)` |
| `udp.db_axial` | `list` | Axial orbit database rows | `loaded (300 rows)` |
| `udp.t` | `array` | Time grid for signal evaluation (days) | `generated (2710 samples)` |

**Target signal helpers**

```python
message = " Tomatoes for sale "
t, target_signal = udp.morse_to_signal(message)
morse_string = udp.message_to_morse_string(message)
print(morse_string) 

```

```python
"- --- -- .- - --- . ... / ..-. --- .-. / ... .- .-.. ."
```

![Morse Code Signal](assets/morse-code-signal.png)

**Period conversion**

```python
period_days = period_rad / (2.0 * np.pi) * udp.synodic_period
```

**CR3BP dynamics** *(optional — not needed to solve the problem)*

The equations of motion for the Circular Restricted Three-Body Problem (CR3BP)
used to propagate the orbits are implemented in `_create_ta()` inside the `udp`.
These are provided in case you wish to propagate orbits forward in time (see `propagate_orbit()`) or
visualise spacecraft trajectories in the Earth-Moon rotating frame (see `plot_orbits()`).
They are **not required** to construct or evaluate a solution — the orbit databases
already contain all timing information needed.

### Orbit Database Row Format

For each family (`udp.db_dro`, `udp.db_lyap`, `udp.db_axial`), each row corresponds to a different periodic orbit, and defines the state at `t=0`, the period and the frequency and amplitude of the occultation periods.
The first occultation occurs at `t=0` and subsequent occultations occur at `t=\pi`

```python
[x, y, z, vx, vy, vz, period_rad, occultation_frequency, occultation_amplitude]
```

Example:

```python
row = udp.db_dro[i]
x, y, z, vx, vy, vz, period_rad, occ_amp, occ_freq = row[:9]
period_days = period_rad / (2.0 * np.pi) * udp.synodic_period
```

---

## Description of the Solution Format

A candidate solution is a flat Python vector `x`:

```python
N = udp.num_dro + udp.num_lyap + udp.num_axial
len(x) == N + N * udp.max_per_orbit
```

### Part 1: selection counts (`sel`)

```python
sel = x[:N]
```

- `sel[i]` is an integer in `[0, udp.max_per_orbit]`
- Orbit order is:
  1. DRO
  2. Lyapunov
  3. Axial

Interpretation:

- `sel[i] == 0`: orbit slot `i` unused
- `sel[i] == k`: `k` active occultation events for orbit slot `i`

### Part 2: phase values (`phs`)

```python
phs = x[N:]
```

For orbit slot `i`, phase block is:

```python
start = i * udp.max_per_orbit
stop = start + udp.max_per_orbit
phase_block = phs[start:stop]
active_phases = phase_block[:sel[i]]
```

- Each active phase is in `[0, 2*np.pi)`
- Inactive phase slots are ignored (can be left as `0.0`)

---

## Definition of the Objective Function

Evaluate with:

```
obj, num_selected, mse = udp.fitness(x, postprocess=True)
```

Outputs:

- `obj[0]`: scalar score used by the optimizer
- `num_selected`: total selected spacecraft
- `mse`: reconstruction mean-squared error

Implemented ranking logic:

1. Feasible (`mse <= mse_thresh`) beats infeasible.
2. Among feasible solutions, lower `num_selected` is better.
3. Among infeasible solutions, lower penalty (MSE-based) is better.

Feasibility condition:

```
mse <= mse_threshold   # typically 0.05 or 0.10
```

---

## How to Submit

```python
from _submission_helper import create_submission

create_submission(
    "spoc-4-luna-tomato-advertising",
    "tie-breaker",
    best_x,
    "solutions/my_submission.json",
)
```

This writes a submission JSON and uploads if your credentials are configured.

---

## Hints, Tips, First Steps

### 1) Build a baseline solution

Start feasible first (even if many spacecraft), then improve.

### 2) Tune phases

With `sel` fixed, optimize `phs` to reduce MSE.

### 3) Validate often

```python
obj, num_selected, mse = udp.fitness(x, postprocess=True)
udp.plot_signal(x)
udp.plot_orbits(x)
```

### 4) Keep feasibility margin

During iterative edits, keep:

```
mse < 0.05
```

so intermediate solutions stay submittable.
