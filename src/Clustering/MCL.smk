

rule MCL:
    """
    Note: hardcoded to keeping 80 %
    """
    input:
        edges="work/edgelists/joined/{joinedname}_0.8.csv"
    output:
        table="work/edgelists/clustering/table/{joinedname}_{inflation}_MCL.clusters"
    shell:
        """
        mcl {input.edges} --abc -I {wildcards.inflation} -o {output}
        """


rule format_output:
    input:
        clusters = "work/edgelists/clustering/table/{joinedname}_{inflation}_MCL.clusters"
    output:
        csv ="work/edgelists/clustering/table/{joinedname}_{inflation}_MCL.csv"
    run:
        with open(output.csv, "w") as w:
            w.write("term\tcluster\n")
            with open(input.clusters, "r") as f:
                cluster = 0
                for l in f:
                    terms = l.strip().split("\t")
                    for term in terms:
                        w.write(f"{term}\t{cluster}\n")
                    cluster +=1

rule MCL_cutoff:
    """
    Note: found 1.6 to be a good value 
    """
    input:
        edges="work/edgelists/joined/{joinedname}_{cutoff}.csv"
    output:
        table="work/edgelists/clustering/table/{joinedname}_{cutoff}_cutoffMCL.clusters"
    shell:
        """
        mcl {input.edges} --abc -I 1.6 -o {output}
        """

rule format_output_cutoff:
    input:
        clusters = "work/edgelists/clustering/table/{joinedname}_{cutoff}_cutoffMCL.clusters"
    output:
        csv ="work/edgelists/clustering/table/{joinedname}_{cutoff}_cutoffMCL.csv"
    run:
        with open(output.csv, "w") as w:
            w.write("term\tcluster\n")
            with open(input.clusters, "r") as f:
                cluster = 0
                for l in f:
                    terms = l.strip().split("\t")
                    for term in terms:
                        w.write(f"{term}\t{cluster}\n")
                    cluster +=1

