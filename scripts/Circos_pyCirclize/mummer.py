from pycirclize import Circos
from pygenomeviz.align import MUMmer
from pygenomeviz.parser import Fasta
import os
TICKS_INTERVAL = 1000000

# Load query & reference 
#ref = Fasta(r"Path/to/reference/genome.fasta")
#query = Fasta(r"Path/to/query/genome.fasta")
# Naming
#ref_name = "name on plot"
#query_name = "name on plot"


ref = Fasta(r"Path/to/reference/genome.fasta")
query = Fasta(r"Path/to/query/genome.fasta")
ref_name1 = "name on plot"
query_name1 = "name on plot"


# Initialize circos instance with proper ordering by size
ref_sectors_by_size = dict(sorted(ref.get_seqid2size().items(), key=lambda x: x[1], reverse=True))
query_sectors_by_size = dict(sorted(query.get_seqid2size().items(), key=lambda x: x[1], reverse=True))

circos = Circos(
    sectors=dict(**ref_sectors_by_size, **dict(reversed(list(query_sectors_by_size.items())))),
    start=-358,
    end=2,
    space=4,
    sector2clockwise={seqid: False for seqid in query_sectors_by_size.keys()},
)

# Names on the plot
#ref
circos.text(f"{ref_name1}\n({ref.full_genome_length/1000000:.1f} Mb)", r=140, deg=35, size=10)
#query
circos.text(f"{query_name1}\n({query.full_genome_length/1000000:.1f} Mb)", r=140, deg=-35, size=10)


# Create mapping from contig name to number based on circos sector order
ref_contig_names = list(ref.get_seqid2size().keys())
query_contig_names = list(query.get_seqid2size().keys())

# Create numbering based on the order they appear in circos.sectors
ref_numbering = {}
query_numbering = {}

ref_counter = 1
query_counter = len(query_contig_names)  # Start from the total number of query contigs

for sector in circos.sectors:
    if sector.name in ref_contig_names:
        ref_numbering[sector.name] = ref_counter
        ref_counter += 1
    elif sector.name in query_contig_names:
        query_numbering[sector.name] = query_counter
        query_counter -= 1  # Decrement for reverse numbering



# Plot genomic sector axis & xticks
for sector in circos.sectors:
    track = sector.add_track((99.8, 100))
    track.axis(fc="black")
    if sector.size >= TICKS_INTERVAL:
        track.xticks_by_interval(
            TICKS_INTERVAL,
            label_formatter=lambda v: f"{v/1000000:.1f} Mb",
            label_orientation="vertical",
        )

     # Add numbered contig labels
    if sector.name in ref_numbering:
        # Reference genome contigs
        sector.text(f"{ref_numbering[sector.name]}", size=8, r=120, orientation="horizontal")
    elif sector.name in query_numbering:
        # Query genome contigs
        sector.text(f"{query_numbering[sector.name]}", size=8, r=120, orientation="horizontal")


# MUMmer genome comparison & plot links
align_coords = MUMmer([query, ref]).run()
for ac in align_coords:
    region1 = (ac.query_name, ac.query_start, ac.query_end)
    region2 = (ac.ref_name, ac.ref_start, ac.ref_end)
    color = "red" if ac.is_inverted else "grey"
    circos.link(region1, region2, color=color, r1=98, r2=98)

#Output file name 
ref_name_file = os.path.splitext(os.path.basename(ref.name))[0]
query_name_file = os.path.splitext(os.path.basename(query.name))[0]
output_filename = f"circos_plot_{ref_name_file}_vs_{query_name_file}.png"
fig = circos.plotfig()
fig.savefig(output_filename, dpi=300)  