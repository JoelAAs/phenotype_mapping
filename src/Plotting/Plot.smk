import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

rule seaborn_range_plot:
    input:
        data_range = "work/edgelists/clustering/results/{joinedname}_{method}.csv"
    output:
        range_plot = "work/edgelists/plots/{joinedname}_{method}.png"
    run:
        data_df = pd.read_csv(input.data_range, sep="\t")
        plot_df = pd.melt(
            data_df,
            value_vars=["med_score", "hpo_score"],
            id_vars="param", value_name="score")
        sns.set_theme(style="whitegrid")
        score_plt = sns.relplot(
            data=plot_df,
            x ="param",
            y = "score",
            hue="variable",
            kind="line")
        score_plt.legend.remove()
        if wildcards.method == "MCL":
            score_plt.set_xlabels("Inflation coefficient")
        else:
            score_plt.set_xlabels("% lowest weights removed")
        score_plt.set(ylim=(0, 1.05))
        plt.legend(loc="lower right")
        plt.savefig(output.range_plot, dpi=500)
