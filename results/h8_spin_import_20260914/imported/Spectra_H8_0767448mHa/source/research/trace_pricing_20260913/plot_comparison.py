"""Export a scientific comparison of actual incumbent certificates."""
from pathlib import Path
import argparse
import json


def incumbent_curve(points, key, end):
    grouped = {}
    for p in points:
        grouped[p[key]] = min(grouped.get(p[key], float('inf')), p['width_mHa'])
    xs = []; ys = []; best = float('inf')
    for x in sorted(grouped):
        best = min(best, grouped[x]); xs.append(x); ys.append(best)
    if xs[-1] < end:
        xs.append(end); ys.append(best)
    return xs, ys


def run(results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data = json.loads((results/'fresh_replay.json').read_text())
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), sharey=True)
    colors = {'coefficient': '#69798c', 'trace': '#126caa'}
    labels = {'coefficient': 'Coefficient-ranked control', 'trace': 'Trace-normalized ranking'}
    for ax, key, end, xlabel in zip(axes, ('directions', 'complete_at_seconds'), (256, 360),
            ('Combination-count ceiling', 'Incremental time ceiling (seconds)')):
        for name, points in data['comparison_points'].items():
            x, y = incumbent_curve(points, key, end)
            ax.step(x, y, where='post', color=colors[name], label=labels[name], linewidth=2)
            observed_x = [p[key] for p in points]; observed_y = [p['width_mHa'] for p in points]
            ax.scatter(observed_x, observed_y, s=11, color=colors[name], zorder=3)
        ax.axhline(1.6, color='#a8503d', linestyle='--', linewidth=1.3, label='1.6 mHa target')
        ax.set_xlabel(xlabel); ax.grid(alpha=.18); ax.spines[['top', 'right']].set_visible(False)
        ax.set_ylim(0, 8); ax.set_xlim(64 if key == 'directions' else 0, end)
    axes[0].set_ylabel('Best certified interval found (mHa)')
    axes[1].legend(fontsize=8, frameon=False, loc='upper right')
    fig.suptitle('H6: same 64-combination seed, same six-minute budget', fontsize=13, x=.06, ha='left')
    fig.text(.06, .01, 'Steps show completed certificates under each ceiling. One historical comparison; frozen ten-pattern Hamiltonian and tail.', fontsize=8, color='#526172')
    fig.tight_layout(rect=(0, .06, 1, .92))
    fig.savefig(results/'comparison.png', dpi=180); fig.savefig(results/'comparison.svg'); plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--results', type=Path, required=True)
    args = parser.parse_args(); run(args.results)
