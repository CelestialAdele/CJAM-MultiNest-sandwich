# CJAM+MultiNest -- estimation of orbital dynamics of Milky Way globular clusters with Bayesian inference

Utlizing CJAM (https://github.com/lauralwatkins/cjam), MultiNest (https://github.com/JohannesBuchner/PyMultiNest), and globular cluster data from Baumgardt's online Globular Cluster Database (https://people.smp.uq.edu.au/HolgerBaumgardt/globular/), estimate dynamical parameters of MW globular clusters.

Launch Gooey2.py to select a globular cluster, dataset, MGE file, and more. From there, MultiNest is called in CJAM_nestsamp2.py, which calls CJAM from nestsamp_watkins_script2.py.
