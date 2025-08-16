# --- Global configuration -----------------------

from itertools import product

# strength and alpha are multiplied by 100 here to avoid . in filenames
strength = ["000", "050", "100", "150", "200", "250", "300", "350", "400"]
alpha = ["001", "025", "050", "075", "100"]
gl_grid = list(product(strength, alpha))

# --- Rules --------------------------------------

rule all:
    input:
        "figs/fig3.pdf",
        "figs/fig5.pdf",
        "figs/fig6.pdf"

rule fig3:
    input:
        ".r-installed",
        expand(
            "data/network/gl/kna10-mc2-s{gls}_a{gla}.rds",
            gls = [s for s, _ in gl_grid],
            gla = [a for _, a in gl_grid]
        ),
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
        "data/network/gl/kna10-mc2-s{gls}_a{gla}.rds"
    params:
        gls_dec = lambda wildcards: float(wildcards.gls) / 100,
        gla_dec = lambda wildcards: float(wildcards.gla) / 100
    shell:
        "bin/lasso gl kna10 mc2 {params.gls_dec} {params.gla_dec}"

rule inferelator_bbsr:
    output: "data/network/inferelator/.bbsr-complete"
    shell: "bin/batch-inferelator --ppj=2 bbsr kna10 mc2 0,0.5,1,1.5,2,2.5,3,3.5,4"

rule inferelator_enet:
    output: "data/network/inferelator/.enet-complete"
    shell: "bin/batch-inferelator --ppj=2 enet kna10 mc2 0.01,0.25,0.50,0.75,1"

rule setup_r:
    output: ".r-installed"
    script: "scripts/setup-r.R"

