import os

def savefig(fig, name, figdir='figs'):
    fig.savefig(os.path.join(figdir, name + '.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(figdir, name + '.pdf'))