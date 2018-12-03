if __name__ == "__main__":
    import pandas as pd
    import json

    meta = pd.read_csv("papers/carbapenem_manuscript/supplements/samples_meta.original.tsv",sep="\t")

    types = pd.read_csv("mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/ST_types.tsv",sep="\t")
    types = {row["isolate"]:row["ST"] for index,row in types.iterrows()}

    plasmid_incs = json.load(open("mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/plasmid_types.json"))

    ST = []
    INCs = []
    toDrop = []
    for index,row in meta.iterrows():
        isolateName = "carb%03d"%int(row["Internal #"])
        if (isolateName in types):
            ST.append( types[isolateName] )
            if isolateName in plasmid_incs:
                INCs.append( [ ",".join( str(elem) for elem in set(x) ) for x in plasmid_incs[isolateName].values() ] )
            else:
                INCs.append( "" )
        else: 
            ST.append("-")
            INCs.append( "" )
    
    meta["ST"] = ST
    meta["inc"]= INCs

    meta.to_csv("papers/carbapenem_manuscript/supplements/samples_meta.tsv",sep="\t")


