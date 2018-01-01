# DecisionMaking_models
Attempts to model and simulate behavioral & neural data during perceptual decision making. Particular interest is given to how the amplitude of psychophysical kernel, which quantifies how strongly an observer relies on the stimulus to guide a behavior (Nienborg & Cumming, Nature 2009), changes over time within the trial in the simulated model.

Note that these model implementations are generally simplified and not strictly based on the original studies. Basically this repository is to understand how these decision-making models work.

- Signal detection theory
  - Integrated-to-bound model (Kiani et al., Journal of Neuroscience 2008)
  - SDT with weighted-integration model
  - SDT with pooling noise
  - Leaky-integrator model (Ossmy et al., Current Biology 2013) ...all implemented in 'SDTvariants.m'
- Bayesian perception
  - Case of 2AFC (Vincent, Journal of Mathematical Psychology 2015) ...in the forked repository 'bayesian2afc'
  - Probablistic inference by neural sampling (Ralf et al., Neuron 2016) ... in the forked repository 'sampling_decision'
  - Bayesian sampling (Sanborn & Chater, Trends in Cognitive Sciences 2016) ...upcoming
- Stimulus-to-MT model (Yates et al., Nature Neuroscience 2017) ...implemented in 'Yates_simple_example.m' & 'Yates_history.m'
- Attractor model
  - nonlinear integration (general attractor behavior) ...implemented in 'Evidence_Accumulation2.m' & 'plot_EvidenceAccumulation.m'
  - Energy profile with 'valleys' (Wong & Wang et al., Journal of Neuroscience 2006; Wimmer et al., Nature Communications, 2015) ...in the forked repository 'Hierarchical_network'
- Reinforcement learning
  - Partially observable Marcov decision process (Lak et al., Current Biology 2017) ...in the forked repository 'POMDP'
