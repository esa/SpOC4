For the fourth time, the European Space Agency's Advanced Concepts Team
(ACT) presents in cooperation with the [The Genetic and Evolutionary
Computation Conference
(GECCO)](https://gecco-2026.sigevo.org/HomePage) the *Space
Optimization Competition (SpOC)*. Look forward to challenging problems set in a futuristic space mission scenario.

<Image src="images/spoc-banner.jpg" align="center" alt="SpOC 2024"/>

## SpOC 4: Space Logistics

SpOC 4 features three distinct challenges, each with its own flavor:

1. [Luna Tomato Logistics](https://optimise.esa.int/challenge/spoc-4-luna-tomato-logistics/About)

Design a transport network to maximize the delivery of precious
tomatoes from Earth orbits to lunar destinations. Each transfer must
respect strict licensing regulations, but be careful: no orbit or
destination can be used more than once! For advanced competitors,
you’ll need to compute the actual transfer trajectories and maximize
the total mass delivered, taking into account the complexities of
orbital mechanics and transfer capacities. 

[To the challenge.](https://optimise.esa.int/challenge/spoc-4-luna-tomato-logistics/About)

2. [Keplerian Tomato Traveling Salesperson](https://optimise.esa.int/challenge/spoc-4-keplerian-tomato-traveling-salesperson/About)

Disaster! Tomatoes are stranded in lunar orbit after a cargo
mishap. As a logistics expert, you must collect all the tomatoes in
the shortest possible time, planning a route that respects the limits
of your spacecraft’s maneuvering capability. This is a traveling
salesperson problem with a twist: the cost of moving between tomatoes
depends on Keplerian dynamics, and not all transfers are
feasible. Optimize your route, timing, and maneuvers to save the Luna
Tomato Industry! 

[To the challenge.](https://optimise.esa.int/challenge/spoc-4-keplerian-tomato-traveling-salesperson/About)

3. [Luna Tomato Advertising](https://optimise.esa.int/challenge/spoc-4-luna-tomato-advertising/About)

Can you broadcast advertisement across the lunar sky? Configure a
fleet of spacecraft in lunar orbits to reproduce a Morse-code
advertising campaign, using occultations as dots and dashes. Your
goal: use as few spacecraft as possible while keeping the signal crisp
and clear. Choose from Distant Retrograde, Lyapunov, and Axial orbits,
and optimize the phases to match the target waveform. Everybody must
know that tomatoes are for sale. 

[To the challenge.](https://optimise.esa.int/challenge/spoc-4-luna-tomato-advertising/About)


## Competition Structure

SpOC 4.0 contains three distinct problems centered around a futuristic
space mission. Starting from **1 April 2026, End of Day, Anywhere on Earth (AoE)** you have **three** months
to tackle these challenges to secure a spot on the leaderboard, i.e., until
**30 June 2026, End of Day, Anywhere on Earth (AoE)**.

**Detailed technical descriptions for the three challenges to be
solved will be made available on the [Optimise
platform](https://optimise.esa.int) from the same date.** 

## Guidelines and Rules

The competition will be hosted on the
  [Optimise](https://optimise.esa.int/) platform developed by the
  Advanced Concepts Team. Participants will need to register online on the
  platform, and solution entries will need to be submitted via
  Optimise for validation. While SpOC is organized in cooperation with
  GECCO 2026, it is *not* required to attend GECCO 2026 in order to
  participate in SpOC.

- Your objective is to propose and implement metaheuristic algorithms
  to solve the proposed optimisation challenges. 
- In order to validate your solutions, we will provide you with Python
  validation code for each of the challenges. This code includes
  problem definitions in the [Pygmo](https://esa.github.io/pygmo2/#)
  user-defined problem (UDP) format, examples of solutions, and
  visualisation tools. 
- You have until **30 June 2026** to submit your entries via the dedicated portal [Optimise](https://optimise.esa.int/).
- Please comply with our [basic code of
  honour](https://optimise.esa.int/terms). The ACT reserves the right
  to exclude users from the competition if they abuse the evaluation
  system.

## Scoring and Winner Selection 

This year, SpOC contains three challenges, of which two are mandatory
and one is a tie-breaker. In detail, you will obtain a *local score*
$s_i$ for the [Luna Tomato Logistics](https://optimise.esa.int/challenge/spoc-4-luna-tomato-logistics/About) challenge and for the [Keplerian Tomato Traveling Salesperson](https://optimise.esa.int/challenge/spoc-4-keplerian-tomato-traveling-salesperson/About) challenge,
each computed via the the rules below,
then your *global score* is the sum of the two. The global score defines your place on the overall SpOC
leaderboard. Ties are broken with the score of the [Luna Tomato Advertising](https://optimise.esa.int/challenge/spoc-4-luna-tomato-advertising/About) challenge. 

The two main challenges contain *three problems* of different
difficulty. Every problem has its own leaderboard (visible on
[Optimise](https://optimise.esa.int/)) that ranks participants
according to the objective of the challenge. The top ten ranks on the
leaderboard of an easy instances get $e_i=10,9,8,\dots,1$ points, similarly
$m_i=\frac{4}{3}10,\frac{4}{3}9,\dots,\frac{4}{3}$ points for the top
ten teams on medium problems, and 
$h_i=(\frac{4}{3})^2 10,(\frac{4}{3})^2 9,\dots,(\frac{4}{3})^2$ points for hard
problems. The score of the challenge is determined by the sum of the
scores obtained in the individual problems.

We wish all participants the best of luck and are excited to see what you accomplish!


## Timeline

The Space Optimization Competition starts at the **first of April, End
of Day, Anywhere on Earth (AoE):** 

<table width="100%" border='0' cellspacing="0" cellpadding="0"><tr><td align="center"><img src="https://i.countdownmail.com/4xe8n9.gif" style="display:inline-block!important;width:100%!important;max-width:586px!important;" border="0" alt="countdownmail.com"/></td></tr></table>

- 1 April, End of Day, AoE Submissions open
- 30 June, End of Day, AoE Submissions closes
- GECCO 2026, Winner Announcement

**NOTE**: The submission portal remains open after **30 June
2026, End of Day, Anywhere on Earth (AoE)**. Submissions received after that date will not be taken into
consideration for the competition, but still appear on the
leaderboard. 

## Contact

Our primary means of communication with competitors will be the
*Discussions* feature on this repository. Please use it to ask
any questions you may have about the challenges or to exchange
information with us. We will do our best to respond to your questions
in a timely manner. 

If you encounter a bug in the code, please use this repository's
*Issues* feature to report it. 
