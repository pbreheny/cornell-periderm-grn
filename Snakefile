# --- Global configuration -----------------------

from itertools import product

strength = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
alpha = [0.01, 0.25, 0.50, 0.75, 1]

grid = list(product(strength, alpha))

# --- Rules --------------------------------------

rule all:
    input:
        "figs/fig3.pdf",
        "figs/fig5.pdf",
        "figs/fig6.pdf"

rule fig3:
    input:
        ".r-installed",
        "data/network/gl/.gl-complete",
        "data/network/inferelator/.bbsr-complete",
        "data/network/inferelator/.enet-complete"
    output: "figs/fig3.pdf"
    script: "scripts/fig3.R"

rule fig5:
    input:
        ".r-installed"
    output: "figs/fig5.pdf"
    script: "scripts/fig5.R"

rule fig6:
    input:
        ".r-installed"
    output: "figs/fig6.pdf"
    script: "scripts/fig6.R"

rule glasso:
    input:
        ".r-installed"
    output:
        touch("data/network/gl/.gl-complete")
    run:
        for strength, alpha in grid:
            shell(f"bin/lasso gl kna10 mc2 {strength} {alpha}")

rule inferelator_bbsr:
    output: "data/network/inferelator/.bbsr-complete"
    shell: "bin/batch-inferelator --ppj=2 bbsr kna10 mc2 0,0.5,1,1.5,2,2.5,3,3.5,4"

rule inferelator_enet:
    output: "data/network/inferelator/.enet-complete"
    shell: "bin/batch-inferelator --ppj=2 enet kna10 mc2 0.01,0.25,0.50,0.75,1"

rule setup_r:
    output: ".r-installed"
    script: "scripts/setup-r.R"

