<image src="../images/ltl.png" align="center" alt="Challange
1: Luna Tomato Logistics"/>

# The Luna Tomato Logistics Challenge

In the near future, as humanity’s first permanent lunar settlements
stretched their silver domes across the dusty lunar surface, a
surprisingly earthly problem emerged: the desire for a good tomato
sauce. While oxygen, fuel, and machinery could be synthesized or
recycled, nothing replaced the rich taste of tomatoes grown under the
blue skies of Earth. What began as a small luxury quickly became a
logistical obsession — how do you reliably transport fragile,
perishable cargo across 384,400 kilometers of space? 

## Description

A logistic transport network has to be designed and planned between
nodes representing orbits in the cislunar space. Some orbits around
the Earth (E) need to be connected to orbits around the Moon (L for
Luna) which are then connected to final destinations (D for
Destination or Depot). An overall network must be designed that
maximizes the amount of material that can be transferred from Earth
orbits to final destinations via Moon orbits.

The Luna Tomato Logistics GmbH (LTL, not to be confused with linear
temporal logic) provides a well-developed infrastructure from Earth
orbits (E) to destinations (D) via Moon orbits (L). Each connection in
this network allows you to transport a fixed amount of tomatoes to a
destination on the Moon. However, each orbit in E and L, and every
destination in D can be only part of one transfer due to licensing
traffic regulations. 

## Warm-Up Problems

Your first task (problems marked as **Beginner**) in this challenge is
to figure out the best selection of transfers in order to *maximize*
the number of tomatoes delivered to Moon destinations while being
consistent with the licensing regulations. 

### Formal Specification 

Three sets $E, L$, and $D$ are considered. They contain nodes (orbits)
of possible futuristic infrastructures around the Earth, around the
Moon, as well as destinations on the Moon surface. A *transfer* is a
connection from an orbit in $E$ to a destination from $D$ via an orbit
in $L$, i.e., a triple $(e,l,d)$ with $e\in E$, $l\in L$, and $d\in
D$. Given a set $T\subseteq E\times L\times D$ of transfers, 
a *solution* is a set $S \subseteq T$ such that for any two distinct
triples $(e_1, l_1, d_1)$ and $(e_2, l_2, d_2)$ one has $e_1\ne e_2$,
$l_1\ne l_2$, and $d_1\ne d_2$. In words, each orbit and every destination can only be used once. The *fitness function* for $S$
is computed by summing up all the mass delivered via the selected
transfers. 

The input for the **Beginner** problems is a set $T$ together with
a weight function $w\colon T\rightarrow \mathbb{R}$ that indicates the
mass of tomatoes that can be moved to the Moon using the given
transfer. We encode the Earth and Moon orbits as well as the
destinations by simple integers, i.e., $E = \{1, 2, \dots\}$, $L =
\{1, 2, \dots\}$, $D = \{1, 2, \dots\}$. The input is specified by a
text file that contains one line per transfer in $T$, specifying the
corresponding element in $E$, $L$, $D$, and its weight. For instance,
the following file specifies an instance with 10 transfers.

```
1 1 1 6.202
2 2 1 6.386
3 2 2 2.641
4 3 2 0.577
5 3 3 0.254
6 4 3 8.042
7 4 4 8.221
8 5 4 8.800
9 5 5 0.213
10 6 5 4.275
```

A solution is represented by a binary array of length $|T|$, where the
$i$th index is $1$ if the $i$th transfer should be part of the
solution (and is $0$ otherwise). Some sample solutions for the above
instance are:

- [0,0,0,0,0,0,0,0,0,0] of weight $0.0$ (no transfer selected)
- [1,0,0,0,0,0,0,0,0,0] of weight $6.202$ (only first transfer selected)
- [1,0,1,0,0,0,0,0,0,0] of weight $8.843$ (first and third transfer selected)
  
The following solutions are *invalid* since they conflict with
the licensing regulation:

- [1,1,0,0,0,0,0,0,0,0] is invalid since destination $1$ is used twice
- [0,1,1,0,0,0,0,0,0,0] is invalid since Moon orbit $2$ is used twice

Invalid solutions will be scored with $0$.

**Note:** Since PyGMO internally always assumes a minimization problem, the fitness function will be shown negated on the score board.

## The Real Problem: Computing the Transfers

Unfortunately, the previous problem was drastically simplified. For
the real problem (marked as **Advanced**), the shape of a transfer
first needs to be computed, which will ultimatively define the
amount of tomatoes the transfer can deliver to the Moon. While LTL
provides a transfer network from Moon orbits to destinations, you will
need to design impulsive transfers (max. 3 impulses) from Earth
orbits to Moon orbits in order to maximise the overall amount of
transferrable tomatoes. As before, every Earth and Moon orbit, as well
as every destination is allowed to be used only once.

### Formal Specification

For this problem, the goal is to transfer as many tomatoes as possible
to Moon destinations **within 200 days**. LTL provides a network that
allows to transport a fixed amount of tomatoes per day from Moon
orbits to destinations (these rates are defined in `LTL.txt`). The
orbits around the Earth and the Moon are defined by triplets $(a,e,i)$
which represent their Keplerian orbital elements with respect to a
reference frame placed on the corresponding body (in
`Earth_orbits.txt` and `Moon_orbits.txt`). These orbits are
controlled rigorously by LTL, hence their important parameters are not
changing with time. 

A transfer from an Earth orbit $e$ to a Moon orbit $l$ is computed as
a multiple impulsive rendevouz in the BCP problem. The spacecraft wet
mass is $m_w=5000$ kg, the system mass $m_{dry}=500$ kg and the mass
of material delivered along a designed transfer is thus computed as
$m_l = m_w \exp^{-\frac{\Delta V_{tot}}{I_{sp} g_0}} -
m_{dry}$. $I_{sp} = 311$ (s), $g_0 = 9.80665$ (m/s$^2$). 

The dynamics to be used is the BCP (bi-circular problem) as stated by Simo' in his 1995 paper and including Earth, Moon and the Sun:

Simó, Carles, et al. "The bicircular model near the triangular libration points of the RTBP." From Newton to Chaos: modern techniques for understanding and coping with chaos in n-body dynamical systems. Boston, MA: Springer US, 1995. 343-370.

The equations of motion are:

$$
\left\{
   \begin{array}{l}
       \dot{\mathbf r} = \mathbf v \\
       \dot v_x = 2v_y + x - (1 - \mu) \frac{x + \mu}{r_1^3} - \mu \frac{x + \mu - 1}{r_2^3} - \frac{\mu_s}{r_s^3} (x- \rho_s \cos(\omega_s t)) - \frac{\mu_s}{\rho_s^2} \cos(\omega_s t)\\
       \dot v_y = -2 v_x + y - (1 - \mu) \frac{y}{r_1^3} - \mu \frac{y}{r_2^3} - \frac{\mu_s}{r_s^3} (y - \rho_s \sin(\omega_s t)) - \frac{\mu_s}{\rho_s^2} \sin(\omega_s t)\\
       \dot v_z = -(1 - \mu) \frac{z}{r_1^3} - \mu \frac{z}{r_2^3} - \mu_s \frac{z}{r_s^3} 
\end{array}\right.
$$

The values of the various parameters in the above equations are given in the following table:


| parameter  | value (non-dimensional units) |
| ------     |                        ------ |
| $\mu$      |           0.01215058439470971 |
| $\mu_s$    |          3.3294604877306713E5 |
| $\rho_s$   |                  3.88811143E2 |
| $\omega_s$ |               -9.25195985E-01 |

The actual SI units used to interpret the non dimensional values are:

| Quantity | value              |
| ------   | ------             |
| $L$      | 3.84405000E8 (m)   |
| $T$      | 3.7567696752E5 (s) |
| $V$      | $\frac LT$ (m/s)   |

#### The Earth Orbits (E)

The allowed starting orbits around the Earth are defined via their
Keplerian elements $a,e,i$ as measured in $t_0$ in an Earth centered
reference frame having the $xy$ plane overlapped to the Earth-Moon
orbit plane and the $x$ axis aligned at $t_0$ with the Earth-Moon
direction. The remaining elements are considered to be not important
and thus need not to be targeted. Note that this defnition is very
unusual as it changes the meaning of "inclination" we are familiar
with, but LTL has its own proprietary ways and we must stick with
them. 

#### The Moon Orbits (L)

The allowed transit orbits around the Moon are defined via their
Keplerian elements $a,e,i$ as measured in $t_0$ in a Moon centered
reference frame having the $xy$ plane overlapped to the Earth-Moon
orbit plane and the $x$ axis aligned at $t_0$ with the Earth-Moon
direction. The remaining elements are considered to be not important
and thus need not to be targeted. Note that this defnition is very
unusual as it changes the meaning of "inclination" we are familiar
with, but LTL has its own proprietary ways and we must stick with
them. 

#### The Solution Format

A solution is a subset $S \subseteq E\times L\times D$ as well as the
information to compute the Earth-Moon trajectories. Therefore, every
element in $S$ is specified by 21 elements:

```
e_id, l_id, d_id, t0, x0 , y0, z0, vx0, vy0, vz0, DVx0, DVy0, DVz0, DVx1, DVy1, DVz1, DVx2, DVy2, DVz2, T1, T2
```

two impulse trajectories can be included by nullifying the mid-course
impulse (DVx1, DVy1, DVz1). The solution vector (chromosome) for this
problem allows you to select up to 400 transfers, so the array should
have length $400\cdot 21=8400$, which will be interpreted in chunks of
size 21. Solutions that uses less than 400 transfers can deactivate
unused entries by setting $e_id$ to $-1$ (the remaining 20 elements of
this transfer than have no meaning).

#### The Fitness Function

The fitness function for $S$ is computed by summing up all the mass
delivered via the selected transfers, each computed first by
evaluating the mass $m_l$ delivered from the Earth orbit to the Moon
orbit and then discounting it via the formula $m_d = \min(m_l, (200 -
\Delta T) c_{ld})$. The coefficients $c_{ld}$ represent the capacity
(kg per day) provided by LTL to transport mass from orbit $l\in L$ to
destination $d \in D$. $\Delta T$ represent the Earth-Moon transfer
duration in days. 

**Note:** Since PyGMO internally always assumes a minimization problem, the fitness function will be shown negated on the score board.

#### Validation Tolerances

Some generous tolerances are used when verifying that the orbits
acquired are actually having the correct semi-major axis, eccentricity
and inclination. These tolerances are 1e-6 on all values. The
semi-major axes are checked on their non dimensional values. 

### Submitting

To submit a solution,  prepare a simple text / JSON file with the following content:

```
  [
    {
      "decisionVector": [ x1, x2, ... , ... ],
      "problem": "<problem>",
      "challenge": "spoc-4-luna-tomato-logistics",
    }
  ] 
```

The decision vector is defined above (note that it is different from
 **Beginner** and **Advanced** problems). Mind the `<problem>` that needs to be replaced by the actual problem id. Once the file is ready, you can [submit it](https://optimize.esa.int/submit). In case you use Python, you can also use
 the [submission helper](https://api.optimize.esa.int/data/tools/submission_helper.py) via

```python
from submisson_helper import create_submission
create_submission("spoc-4-luna-tomato-logistics","<problem>",x,"submission_file.json","submission_name","submission_description")
```

### Utilities / Hints

* We will be using [GitHub](https://github.com/esa/SpOC4) as our hub for communication with competitors. Our primary means of communication will be the [Discussions](https://github.com/esa/SpOC4/discussions) feature on our repository.

### Instance Files

- [Matching I](https://api.optimize.esa.int/data/spoc4/ltl/random_harder_5000.txt)
- [Matching II](https://api.optimize.esa.int/data/spoc4/ltl/random_hard_10000.txt)
- Trajectory Matching
  - [Earth Orbits](https://api.optimize.esa.int/data/spoc4/ltl/Earth_orbits.txt)
  - [Moon Orbits](https://api.optimize.esa.int/data/spoc4/ltl/Moon_orbits.txt)
  - [LTL Rates](https://api.optimize.esa.int/data/spoc4/ltl/LTL.txt)
