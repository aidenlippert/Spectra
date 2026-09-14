# Controlled asymmetric quartic reference

A full quartic reference for epsilon=1/1000 **completed** inside a measured 120-second process deadline. It did not produce an accurate numerical optimum. Its exact exported certificate remains valid.

The Hamiltonian has five distinct pair hoppings and conserves pair charges plus simultaneous left/right exchange. The reference uses those symmetries only: 581 original dictionaries of maximum size 186 split into 1,122 exchange eigenspaces of maximum size 96, with 59,181 PSD scalar variables. It does not use flavor permutation symmetry.

The controlled parent process recorded the worker as live at ten-second intervals, then recorded normal exit at **106.58 seconds**. Assembly took 8.91 seconds and the numerical solve took 74.99 seconds; remaining time includes factor export and exact checks.

The solver returned `optimal_inaccurate` with objective **3.281055788864656**. The independently checked upper bound is **3.2810047778834033**. The numerical objective exceeds that upper bound by **5.1011e-5**, so it cannot be interpreted as a valid lower bound or accurate relaxation optimum.

After clipping numerical Gram eigenvalues, rational export, and exact coefficient replay, the accepted interval is

```text
3.2761969937601374 <= E0 <= 3.2810047778834033
width = 0.004807784123265976
```

This interval is weaker than the smaller adapted proof. That comparison diagnoses numerical accuracy limitations of this run; it does not prove the full quartic cone is weak or that the asymmetric quartic optimum differs from the physical optimum.

Accepted artifacts are under `results/marginal_asymmetric_reference/controlled_1_1000/`, including process observations, worker log, exact certificate, and independent replay. Reproduction uses `python -m experiments.marginal_asymmetric_reference --epsilon 1/1000 --budget 120`; the controller refuses to overwrite an existing controlled run. Exact replay uses `python3 -S -m experiments.marginal_transfer_verify results/marginal_asymmetric_reference/controlled_1_1000/certificate.json`.

Earlier uncontrolled observations at the parent results directory are superseded. An initial receipt recorded an 11.08-second interrupted observation; a same-command process was later observed alive for more than three minutes and explicitly terminated by the root agent. Those observations do not form a reliable solver timing experiment and do not establish a timeout or scalability conclusion. The controlled run above replaces those claims.
