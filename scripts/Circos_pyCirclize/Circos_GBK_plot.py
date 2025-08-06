from pycirclize import Circos
from pycirclize.parser import Genbank
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import pandas as pd

# Load Genbank file
gbk_file = "PATH/TO/GBK/File"
telomere_tsv_file = "PATH/To/telomeric/repeat/windows.tsv"
gbk = Genbank(gbk_file)

Output_name = "Output_name"

# Initialize circos instance
seqid2size = gbk.get_seqid2size()
space = 0 if len(seqid2size) == 1 else 2
circos = Circos(sectors=seqid2size, space=space)
#circos.text("Aspergillus oryzae STRAINNR", size=12, r=20)
circos.text("NAME ON PLOT", size=12, r=20)
seqid2features = gbk.get_seqid2features(feature_type=None)
seqid2seq = gbk.get_seqid2seq()

# Add tracks for Genbank features
for sector in circos.sectors:
    # Plot outer track
    outer_track = sector.add_track((98, 100))
    outer_track.axis(fc="lightgrey")
    # Add start and end ticks manually
    start_position = 0
    end_position = sector.size
    start_label = "0 Mb"
    end_label = f"{sector.size / 10**6:.1f} Mb"
    outer_track.xticks(
        [start_position, end_position], 
        [start_label, end_label],        
        label_orientation="vertical",
        tick_length=1,
    )
    # Plot chromosome/contig name
    sector.text(sector.name, size=10)
    # Add tracks for features
    f_cds_track = sector.add_track((90, 97), r_pad_ratio=0.1)
    r_cds_track = sector.add_track((83, 90), r_pad_ratio=0.1)
    trna_track = sector.add_track((76, 83), r_pad_ratio=0.1)

    # Plot Forward CDS, Reverse CDS, tRNA, BGC
    features = seqid2features[sector.name]
    for feature in features:
        if feature.type == "CDS" and feature.location.strand == 1:
            f_cds_track.genomic_features(feature, fc="red")
        elif feature.type == "CDS" and feature.location.strand == -1:
            r_cds_track.genomic_features(feature, fc="blue")
        elif feature.type == "tRNA":
            trna_track.genomic_features(feature, color="magenta", lw=0.1)


###antiSMASH###
# Extract antiSMASH clusters from features
def extract_antismash_clusters(features):
    clusters = {}
    for feature in features:
        if feature.type == "CDS" and "note" in feature.qualifiers:
            notes = feature.qualifiers["note"]
            for note in notes:
                if "antiSMASH:Cluster_" in note:
                    cluster_id = note.split("antiSMASH:")[1].split()[0]  # Extract Cluster ID
                    if cluster_id not in clusters:
                        clusters[cluster_id] = []
                    clusters[cluster_id].append(feature)
    return clusters


# Add a new track for antiSMASH clusters
for sector in circos.sectors:
    antismash_track = sector.add_track((69, 76), r_pad_ratio=0.1)
    features = seqid2features[sector.name]
    clusters = extract_antismash_clusters(features)

    # Plot each cluster
    for cluster_id, cluster_features in clusters.items():
        for feature in cluster_features:
            start = feature.location.start
            end = feature.location.end
            antismash_track.rect(start, end, fc="green", ec="black", label=cluster_id)
###antiSMASH###





###CAZy###
# Extract CAZy enzymes from features
def extract_cazy_enzymes(features):
    cazy_enzymes = {}
    for feature in features:
        if feature.type == "CDS" and "note" in feature.qualifiers:
            notes = feature.qualifiers["note"]
            for note in notes:
                if "CAZy:" in note:
                    cazy_family = note.split("CAZy:")[1].split()[0]  # Extract CAZy family (e.g., GT8, GH78)
                    if cazy_family not in cazy_enzymes:
                        cazy_enzymes[cazy_family] = []
                    cazy_enzymes[cazy_family].append(feature)
    return cazy_enzymes


# Add a new track for CAZy enzymes
for sector in circos.sectors:
    cazy_track = sector.add_track((63, 69), r_pad_ratio=0.1)
    features = seqid2features[sector.name]
    cazy_enzymes = extract_cazy_enzymes(features)

    # Plot each CAZy family
    for cazy_family, cazy_features in cazy_enzymes.items():
        for feature in cazy_features:
            start = feature.location.start
            end = feature.location.end
            cazy_track.rect(start, end, fc="purple", ec="black", label=cazy_family)
###CAZy###




### Telomere ###
# Function to extract telomere data from TSV
def extract_telomere_data(tsv_file):
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(tsv_file, sep="\t")
    ###
    #The telomere_data number i.e. 5 for aspergillus and trichoderma and 20 for neurospora is highly dependant on how you use TIDK, 
    # but it should just be clear that there are alot more telomeric repeats at the ends of the contigs so set your number after that.  

    ###
    # Filter rows where forward or reverse repeat number is greater than 0
    #Aspergillus/Trichoderma
    telomere_data = df[(df["forward_repeat_number"] > 5) | (df["reverse_repeat_number"] > 5)]
    #Neurospora
    #telomere_data = df[(df["forward_repeat_number"] > 20) | (df["reverse_repeat_number"] > 20)]
    # Extract relevant columns: id, window, forward_repeat_number, reverse_repeat_number
    telomere_data = telomere_data[["id", "window", "forward_repeat_number", "reverse_repeat_number"]]
    return telomere_data


# Add a new track for telomeres

telomere_data = extract_telomere_data(telomere_tsv_file)
# Save the extracted telomere data to a new file
output_file = "telomeric_repeat_filtered.tsv"
telomere_data.to_csv(output_file, sep="\t", index=False)
print(f"Telomere data saved to {output_file}")


for sector in circos.sectors:
    # Plot outer track
    outer_track = sector.add_track((98, 100))
    outer_track.axis(fc="lightgrey", ec="black", lw=1)  # Default axis for the track
    telomeres = telomere_data[telomere_data["id"] == sector.name]

    # Add bold orange borders at the start and end of the chromosome if telomeres are present
    for _, row in telomeres.iterrows():
        window_start = int(row["window"])
        if window_start == 10000:  # Forward telomere (start of the chromosome)
            # Add a bold orange rectangle at the start
            outer_track.rect(0, 10000, fc="none", ec="orange", lw=4)
        else:  # Reverse telomere (end of the chromosome)
            # Add a bold orange rectangle at the end
            outer_track.rect(sector.size - 10000, sector.size, fc="none", ec="orange", lw=4)
### Telomere ###



# Plot the circos figure
fig = circos.plotfig()

# Add legend
handles = [
    Patch(color="red", label="Forward CDS"),
    Patch(color="blue", label="Reverse CDS"),
    Patch(color="magenta", label="tRNA"),
    Patch(color="green", label="BGC"),
    Patch(color="Purple", label="CAZy Enzymes"),
    Patch(color="orange", label="Telomere"),
]
_ = circos.ax.legend(handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=8)

# Save the figure
fig.savefig(Output_name, dpi=300)