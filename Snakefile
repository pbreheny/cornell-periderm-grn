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
            "data/network/gl/kna10-mc2-s{gls}-a{gla}.rds",
            gls = [s for s, _ in gl_grid],
            gla = [a for _, a in gl_grid]
        ),
        expand(
            "data/network/inferelator/kna10-mc2-bbsr-s{bbsrs}.rds",
            bbsrs = [s for s in strength]
        ),
        expand(
            "data/network/inferelator/kna10-mc2-enet-a{eneta}.rds",
            eneta = [a for a in alpha]
        )
    output:
        "figs/fig3.pdf"
    script:
        "scripts/fig3.R"

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
        "data/network/gl/kna10-mc2-s{gls}-a{gla}.rds"
    params:
        gls_dec = lambda wildcards: float(wildcards.gls) / 100,
        gla_dec = lambda wildcards: float(wildcards.gla) / 100
    shell:
        "bin/lasso gl kna10 mc2 {params.gls_dec} {params.gla_dec}"

rule inferelator_bbsr:
    input:
        ".inferelator-installed"
    output:
        "data/network/inferelator/kna10-mc2-bbsr-s{bbsrs}.rds"
    params:
        bbsrs_dec = lambda wildcards: float(wildcards.bbsrs) / 100,
    shell:
        "bin/fit-inferelator bbsr kna10 mc2 {params.bbsrs_dec}"

rule inferelator_enet:
    input:
        ".inferelator-installed"
    output:
        "data/network/inferelator/kna10-mc2-enet-a{eneta}.rds"
    params:
        eneta_dec = lambda wildcards: float(wildcards.eneta) / 100,
    shell:
        "bin/fit-inferelator enet kna10 mc2 {params.eneta_dec}"

rule setup_r:
    output: ".r-installed"
    script: "scripts/setup-r.R"

rule setup_inferelator:
    output: ".inferelator-installed"
    shell: "python -m pip install inferelator"

