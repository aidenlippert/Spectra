# Lambda GPU experiment

Authorized spending ceiling: $2,600 total. Start with one A100 40GB SXM4 at the displayed $1.99/hour, Arizona (us-west-2), using the registered aiden-mac SSH key. The user confirmed the final license and billing agreement.

The source archive is a frozen snapshot, with per-file SHA-256 hashes in launch_receipt.json. It contains only Python experiment code and the two required numerical inputs. No credentials are in the archive.

Mac baseline: 64 evaluations, batch 16, ratio 0.5. See mac_baseline/receipt.json. This is a numerical equivalence and performance benchmark, not an exact certificate.

Operational policy: no additional instance until measured throughput warrants it; retrieve results before termination; terminate through Lambda (guest shutdown alone is not treated as stopping billing). Do not leave an idle instance after the experiment unless the active continuation task receives the explicitly authorized final ownership handoff.
