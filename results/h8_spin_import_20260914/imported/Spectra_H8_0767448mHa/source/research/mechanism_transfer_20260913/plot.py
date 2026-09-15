"""Export the held-out interval comparison, including the fitted H6 reference."""
from pathlib import Path
import json


def run():
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=Path(__file__).resolve().parents[2];out=root/'results/mechanism_transfer_20260913'
    data=json.loads((out/'fresh_replay.json').read_text())
    rows={(r['case'],r['arm']):r for r in data['rows']}
    cases=['h6_train','h6_1p2','h6_1p8','h8_heldout'];x=np.arange(4);w=.23
    fig,ax=plt.subplots(figsize=(10.2,4.8))
    for offset,arm,label,color in ((-w,'quadratic','Quadratic baseline','#8994a3'),
        (0,'coupled','Coupled contraction rule','#1676a5'),(w,'separate','Same rule, separate contributions','#bcc7ca')):
        ax.bar(x+offset,[rows[c,arm]['width_mHa'] for c in cases],w,label=label,color=color)
    for i,c in enumerate(cases):
        gain=rows[c,'quadratic']['width_mHa']-rows[c,'coupled']['width_mHa']
        top=max(rows[c,a]['width_mHa'] for a in ('quadratic','coupled','separate'))
        ax.text(i,top+.5,f'Gain: {gain:.3f} mHa',ha='center',fontsize=9)
    ax.scatter([-.4],[3.415421065516425],marker='D',color='#874d32',s=35,zorder=5,label='Fitted 216-combination H6 proof')
    ax.axhline(1.6,color='#874d32',linestyle='--',linewidth=1,label='1.6 mHa target')
    ax.set_xticks(x,['H6 · 1.4 Å\nTraining','H6 · 1.2 Å\nNew geometry','H6 · 1.8 Å\nNew geometry','H8 · 1.4 Å\nHeld out from rule selection'])
    ax.set_ylabel('Certified interval width (mHa)');ax.set_ylim(0,31.5)
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
    ax.legend(loc='upper left',fontsize=8,frameon=False)
    fig.suptitle('The portable rule retains only a small part of the fitted H6 gain',x=.07,ha='left',fontsize=13)
    fig.text(.07,.01,'108 generated combinations for H6; 144 for H8. Fixed rule, with the same Hamiltonian and upper within each comparison.',fontsize=8,color='#526172')
    fig.tight_layout(rect=(0,.04,1,.93));fig.savefig(out/'transfer.png',dpi=180);fig.savefig(out/'transfer.svg');plt.close(fig)


if __name__=='__main__':run()
