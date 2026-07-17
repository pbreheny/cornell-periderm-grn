container: "container/periderm-grn.sif"
from itertools import product


# --- Global configuration -----------------------

# strength and alpha are multiplied by 100 here to avoid . in filenames
strength = ["000", "050", "100", "150", "200", "250", "300", "350", "400"]
alpha = ["001", "025", "050", "075", "100"]
gl_grid = list(product(strength, alpha))


# --- Rules --------------------------------------

rule all:
    input:
        "figs/fig3.pdf",
        "figs/fig5a.pdf",
        "figs/fig5b.pdf",
        "figs/fig5c.csv",
        "figs/fig6.pdf"

rule fig3:
    input:
        expand(
            "network/gl/kna10-mc2-s{gls}-a{gla}.rds",
            gls = [s for s, _ in gl_grid],
            gla = [a for _, a in gl_grid]
        ),
        expand(
            "network/inferelator/kna10-mc2-bbsr-s{bbsrs}.rds",
            bbsrs = [s for s in strength]
        ),
        expand(
            "network/genie/kna10-mc2-gne-s{gs}.rds",
            gs = [s for s in strength]
        )
    output:
        "figs/fig3.pdf"
    script:
        "scripts/fig3.R"

rule fig5:
    input:
        "network/wgcna/kna10-mc2-cor-s100.rds",
        "data/go/cluster-a.json",
        "data/go/cluster-b.json"
    output:
        "figs/fig5a.pdf",
        "figs/fig5b.pdf",
        "figs/fig5c.csv"
    script: "scripts/fig5.R"

rule fig6:
    input:
        "network/wgcna/kna10-mc2-cor-s100.rds"
    output: "figs/fig6.pdf"
    script: "scripts/fig6.R"

rule glasso:
    output:
        "network/gl/kna10-mc2-s{gls}-a{gla}.rds"
    params:
        gls_dec = lambda wildcards: float(wildcards.gls) / 100,
        gla_dec = lambda wildcards: float(wildcards.gla) / 100
    shell:
        "scripts/lasso gl kna10 mc2 {params.gls_dec} {params.gla_dec}"

rule inferelator_bbsr:
    threads: 10
    output:
        "network/inferelator/kna10-mc2-bbsr-s{bbsrs}.rds"
    params:
        bbsrs_dec = lambda wildcards: float(wildcards.bbsrs) / 100,
    shell:
        "scripts/fit-inferelator bbsr kna10 mc2 {params.bbsrs_dec}"

rule inferelator_enet:
    threads: 10
    output:
        "network/inferelator/kna10-mc2-enet-a{eneta}.rds"
    params:
        eneta_dec = lambda wildcards: float(wildcards.eneta) / 100,
    shell:
        "scripts/fit-inferelator enet kna10 mc2 {params.eneta_dec}"

rule wgcna:
    output:
        expand(
            "network/wgcna/kna10-mc2-cor-s{strength}.rds",
            strength = strength
        )
    script: "scripts/wgcna.R"

rule genie:
    threads: 10
    output:
        expand(
            "network/genie/kna10-mc2-gne-s{strength}.rds",
            strength = strength
        )
    script: "scripts/genie.R"
