# H6 factor-row compression audit

This is a post-discovery audit of the downloaded cubic-precision H6 certificate
at [certificate.json](/Users/aidenlippert/Documents/Spectra/results/lambda_runs/cubic_precision/downloaded/precision_scs/h6_none/certificate.json).
It does not rerun discovery or use the physical upper bound to choose factors.

The source has 50 blocks, 1,130 encoded factor rows, 126,278 nonzero integer
entries, and an encoded size of about 949 KiB. Complete rows were removed by
norm threshold and each modified certificate was sent through the existing
exact verifier. Thresholds 0, 10,000, 100,000, and 1,000,000 retained 1,130,
1,118, 717, and 238 rows respectively; all four exact replays were accepted.
Encoded JSON sizes were 1,128,519; 1,120,750; 873,625; and 517,312 bytes.

The reported row-loss quantity is the sum of squared Euclidean norms of
discarded integer rows. It is a bookkeeping bound on removed factor weight,
not an operator-norm error bound or a guarantee that the lower bound remains
unchanged. Exact replay is the acceptance gate. The compressed certificates
therefore demonstrate storage reduction for this already-discovered proof,
not a discovery or scaling theorem.
