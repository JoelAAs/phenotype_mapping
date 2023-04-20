include: "GeneCombinations/GetTermPairs.smk"
include: "GetGeneDepth/GetGeneDepths.smk"
include: "GetGeneDepth/GetUniqueGeneDepthPairs.smk"
include: "GetPossiblePaths/GetPathsAtDepth.smk"
include: "AggregateDistances/CalculateTermDistance.smk"
include: "AggregateDistances/CalculateGeneDistance.smk"
include: "AggregateDistances/AggregateDistanceMartix.smk"



project = config["project"]

### Input all
rule all:
    input:
        "work/{project}/term_distance/all_term_distances.csv".format(project=project)
