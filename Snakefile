rule glasso:
    output: "data/network/gl/.gl-complete"
    shell: "bin/batch-lasso gl kna10 mc2 0,0.5,1,1.5,2,2.5,3,3.5,4 1,25,50,75,100"

rule inferelator_bbsr:
    output: "data/network/inferelator/.bbsr-complete"
    shell: "bin/batch-inferelator --ppj=2 bbsr kna10 mc2 0,0.5,1,1.5,2,2.5,3,3.5,4"

rule inferelator_enet:
    output: "data/network/inferelator/.enet-complete"
    shell: "bin/batch-inferelator --ppj=2 enet kna10 mc2 0.01,0.25,0.50,0.75,1"

rule setup_r:
    output: ".r-installed"
    script: "scripts/setup-r.R"

rule fig4:
    output: "figs/fig4.pdf"
    script: "scripts/fig4.r"

rule fig6:
    output: "figs/fig4.pdf"
    script: "scripts/fig4.r"
