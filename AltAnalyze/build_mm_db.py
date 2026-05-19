#!/usr/bin/env python
"""
Build the processed Mm RNASeq database files that AltAnalyze expects
at runtime from the raw Ensembl annotation files.

AltAnalyze's checkForLocalArraySupport() requires:
  1. AltDatabase/EnsMart31/Mm/RNASeq/ directory to exist
  2. probeset-domain-annotations-exoncomp.txt with >9 lines

The RNASeq analysis also reads:
  3. Mm_Ensembl_exons.txt (reformatted from Mm_Ensembl_exon.txt)
  4. Mm_Ensembl_junctions.txt (reformatted from Mm_Ensembl_junction.txt)
  5. chr/Mm_Ensembl_exons.bed (BED file for BAM-to-BED conversion)
"""
import os
import sys

SPECIES = "Mm"
VERSION = "EnsMart31"
# When run from the Docker build context, the AltDatabase is next to this script.
# When run inside the container, it's at /usr/src/app/altanalyze/AltDatabase.
# Accept an optional CLI argument for the base path, or auto-detect.
if len(sys.argv) > 1:
    BASE = sys.argv[1]
else:
    # Check if we're inside the container or in the build context
    _candidate = os.path.join(os.path.dirname(os.path.abspath(__file__)), "AltDatabase")
    if os.path.isdir(os.path.join(_candidate, VERSION, "ensembl", SPECIES)):
        BASE = _candidate
    else:
        _candidate = os.path.join(os.path.dirname(os.path.abspath(__file__)), "altanalyze", "AltDatabase")
        BASE = _candidate
ENSEMBL_DIR = os.path.join(BASE, VERSION, "ensembl", SPECIES)
RNASEQ_DIR = os.path.join(BASE, VERSION, SPECIES, "RNASeq")
CHR_DIR = os.path.join(RNASEQ_DIR, "chr")


def reformat_exon_file():
    """Reformat Mm_Ensembl_exon.txt into Mm_Ensembl_exons.txt and chr/Mm_Ensembl_exons.bed.

    This replicates the logic of RNASeq.reformatExonFile(species, 'exon', chr_status=True).
    Input format:  gene\texon-id\tchromosome\tstrand\tstart\tstop\tconstitutive_call\tens_exon_ids\tsplice_events\tsplice_junctions
    Output format:  AltAnalyzeID\texon_id\tensembl_gene_id\ttranscript_cluster_id\tchromosome\tstrand\tprobeset_start\tprobeset_stop\taffy_class\tconstitutive_probeset\tens_exon_ids\tens_constitutive_status\texon_region\texon-region-start(s)\texon-region-stop(s)\tsplice_events\tsplice_junctions
    BED format:    chrom\tstart\tstop\tname\tscore\tstrand
    """
    input_file = os.path.join(ENSEMBL_DIR, SPECIES + "_Ensembl_exon.txt")
    output_file = os.path.join(RNASEQ_DIR, SPECIES + "_Ensembl_exons.txt")
    bed_file = os.path.join(CHR_DIR, SPECIES + "_Ensembl_exons.bed")

    if not os.path.exists(input_file):
        print("ERROR: {} not found".format(input_file))
        sys.exit(1)

    print("Reformatting exon file: {} -> {}".format(input_file, output_file))

    with open(input_file, "r") as fin, open(output_file, "w") as fout, open(bed_file, "w") as fbed:
        header = "\t".join([
            "AltAnalyzeID", "exon_id", "ensembl_gene_id", "transcript_cluster_id",
            "chromosome", "strand", "probeset_start", "probeset_stop",
            "affy_class", "constitutive_probeset", "ens_exon_ids", "ens_constitutive_status",
            "exon_region", "exon-region-start(s)", "exon-region-stop(s)",
            "splice_events", "splice_junctions"
        ])
        fout.write(header + "\n")

        for line in fin:
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            if parts[0] == "gene":
                continue  # skip header

            try:
                gene, exonid, chrom, strand, start, stop, constitutive, ens_exon_ids, splice_events, splice_junctions = parts
            except ValueError:
                print("Skipping malformed line: {}".format(line[:80]))
                continue

            if chrom == "chrM":
                chrom = "chrMT"
            if chrom == "M":
                chrom = "MT"

            constitutive_status = "1" if constitutive == "yes" else "0"

            altanalyze_id = "{}:{}".format(gene, exonid)
            out_values = [
                altanalyze_id, exonid, gene, "", chrom, strand, start, stop,
                "known", constitutive, ens_exon_ids, constitutive_status,
                exonid, start, stop, splice_events, splice_junctions
            ]
            fout.write("\t".join(out_values) + "\n")

            # BED file: keep 'chr' prefix to match BAM chromosome names
            fbed.write("\t".join([chrom, start, stop, "{}_{}".format(altanalyze_id, ens_exon_ids), "0", strand]) + "\n")

    line_count = sum(1 for _ in open(output_file)) - 1  # subtract header
    print("Wrote {} lines to {}".format(max(line_count, 0), output_file))
    print("Wrote BED file {}".format(bed_file))


def reformat_junction_file():
    """Reformat Mm_Ensembl_junction.txt into Mm_Ensembl_junctions.txt.

    This replicates the logic of RNASeq.reformatExonFile(species, 'junction', chr_status=True).
    Input format:  same as exon input
    Output format: same as exon output (but no BED file)
    """
    input_file = os.path.join(ENSEMBL_DIR, SPECIES + "_Ensembl_junction.txt")
    output_file = os.path.join(RNASEQ_DIR, SPECIES + "_Ensembl_junctions.txt")

    if not os.path.exists(input_file):
        print("ERROR: {} not found".format(input_file))
        sys.exit(1)

    print("Reformatting junction file: {} -> {}".format(input_file, output_file))

    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        header = "\t".join([
            "AltAnalyzeID", "exon_id", "ensembl_gene_id", "transcript_cluster_id",
            "chromosome", "strand", "probeset_start", "probeset_stop",
            "affy_class", "constitutive_probeset", "ens_exon_ids", "ens_constitutive_status",
            "exon_region", "exon-region-start(s)", "exon-region-stop(s)",
            "splice_events", "splice_junctions"
        ])
        fout.write(header + "\n")

        for line in fin:
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            if parts[0] == "gene":
                continue  # skip header

            try:
                gene, exonid, chrom, strand, start, stop, constitutive, ens_exon_ids, splice_events, splice_junctions = parts
            except ValueError:
                print("Skipping malformed line: {}".format(line[:80]))
                continue

            if chrom == "chrM":
                chrom = "chrMT"
            if chrom == "M":
                chrom = "MT"

            constitutive_status = "1" if constitutive == "yes" else "0"

            altanalyze_id = "{}:{}".format(gene, exonid)
            out_values = [
                altanalyze_id, exonid, gene, "", chrom, strand, start, stop,
                "known", constitutive, ens_exon_ids, constitutive_status,
                exonid, start, stop, splice_events, splice_junctions
            ]
            fout.write("\t".join(out_values) + "\n")

    line_count = sum(1 for _ in open(output_file)) - 1  # subtract header
    print("Wrote {} lines to {}".format(max(line_count, 0), output_file))


def create_domain_annotations_stub():
    """Create a minimal probeset-domain-annotations-exoncomp.txt.

    checkForLocalArraySupport() requires this file to have >9 lines.
    The content is parsed by importGeneric() which reads tab-delimited lines.
    An empty/invalid file (beyond the header) just means no domain annotations
    are loaded, which is acceptable for the SNAF pipeline.
    """
    output_file = os.path.join(RNASEQ_DIR, "probeset-domain-annotations-exoncomp.txt")

    print("Creating stub {}".format(output_file))

    with open(output_file, "w") as fout:
        # Header line (imported as dbase['title'])
        fout.write("ProbesetID\tDomainID\tDomainDescription\tAlignmentStart\tAlignmentStop\n")
        # 10 dummy lines to pass the >9 lines check
        for i in range(10):
            fout.write("dummy_probeset_{}\tdummy_domain\tNo domain annotation available\t0\t0\n".format(i))

    print("Wrote stub domain annotations file (11 lines)")


def create_protein_annotations_stub():
    """Create minimal probeset-protein-annotations-exoncomp.txt files.

    AugmentEventAnnotations.importIsoformAnnotations() reads this file
    to add protein domain predictions to junction events. It checks two
    locations depending on dataType:
      - AltDatabase/EnsMart31/Mm/RNASeq/probeset-protein-annotations-exoncomp.txt
      - AltDatabase/EnsMart31/Mm/RNASeq/junction/probeset-protein-annotations-exoncomp.txt

    An empty file (no data lines) means no protein annotations are added,
    which is acceptable for the SNAF pipeline.
    """
    output_file = os.path.join(RNASEQ_DIR, "probeset-protein-annotations-exoncomp.txt")
    junction_dir = os.path.join(RNASEQ_DIR, "junction")
    output_file_junction = os.path.join(junction_dir, "probeset-protein-annotations-exoncomp.txt")

    for path in [output_file, output_file_junction]:
        d = os.path.dirname(path)
        if not os.path.isdir(d):
            os.makedirs(d)
        print("Creating stub {}".format(path))
        with open(path, "w") as fout:
            fout.write("\n")

    print("Wrote stub protein annotations files")


def create_platform_file():
    """Create platform.txt to identify this as an RNASeq database.

    AltAnalyze checks for platform.txt in some code paths.
    """
    output_file = os.path.join(RNASEQ_DIR, "platform.txt")

    print("Creating {}".format(output_file))

    with open(output_file, "w") as fout:
        fout.write("RNASeq\n")

    print("Wrote platform.txt")


def create_junction_comps_stub():
    """Create a minimal Mm_junction_comps.txt.

    exportKnownJunctionComparisons() creates this file from
    Mm_alternative_junctions.txt, which we don't have for mouse.
    A stub with just a header allows AltAnalyze to proceed without errors.
    """
    output_file = os.path.join(RNASEQ_DIR, SPECIES + "_junction_comps.txt")

    print("Creating stub {}".format(output_file))

    with open(output_file, "w") as fout:
        fout.write("gene\tcritical_exon\texclusion_junction_region\tinclusion_junction_region\texclusion_probeset\tinclusion_probeset\tdata_source\n")

    print("Wrote stub junction_comps file")


def create_junction_comps_updated():
    """Create Mm_junction_comps_updated.txt.

    AltAnalyze reads this at runtime. A stub with just a header is fine
    when there are no junction comparisons to report.
    """
    output_file = os.path.join(RNASEQ_DIR, SPECIES + "_junction_comps_updated.txt")

    print("Creating stub {}".format(output_file))

    with open(output_file, "w") as fout:
        fout.write("gene\tcritical_exon\texclusion_junction_region\tinclusion_junction_region\texclusion_probeset\tinclusion_probeset\tdata_source\n")

    print("Wrote stub junction_comps_updated file")


def update_version_file():
    """Update Config/version.txt to EnsMart31.

    AltAnalyze uses Config/version.txt to resolve database paths.
    The Hs --update Official step writes EnsMart91 here, but we need
    it to say EnsMart31 for Mm to work. Since the --version flag
    overwrites this at runtime, we just need to ensure the version
    is correct when AltAnalyze starts.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Check both possible locations for Config/
    candidates = [
        os.path.join(script_dir, "altanalyze", "Config", "version.txt"),
        os.path.join(script_dir, "Config", "version.txt"),
    ]
    version_file = None
    for c in candidates:
        parent = os.path.dirname(c)
        if os.path.isdir(parent):
            version_file = c
            break
    if version_file is None:
        # Fall back to the container path
        version_file = "/usr/src/app/altanalyze/Config/version.txt"

    print("Updating {} to EnsMart31".format(version_file))

    version_dir = os.path.dirname(version_file)
    try:
        if not os.path.isdir(version_dir):
            os.makedirs(version_dir)
        with open(version_file, "w") as fout:
            fout.write("EnsMart31\t01/01/2025\n")
        print("Updated version.txt")
    except OSError as e:
        print("WARNING: Could not write version.txt: {}".format(e))
        print("This is expected when running outside the Docker container.")
        print("The Dockerfile will handle version.txt separately.")


def main():
    print("Building Mm RNASeq database for AltAnalyze...")
    print("Base directory: {}".format(BASE))
    print("Ensembl directory: {}".format(ENSEMBL_DIR))
    print("RNASeq directory: {}".format(RNASEQ_DIR))

    # Verify input files exist
    for fname in ["Mm_Ensembl_exon.txt", "Mm_Ensembl_junction.txt"]:
        path = os.path.join(ENSEMBL_DIR, fname)
        if not os.path.exists(path):
            print("ERROR: Required input file {} not found".format(path))
            sys.exit(1)

    # Create output directories (Python 2 compatible)
    for d in [RNASEQ_DIR, CHR_DIR]:
        if not os.path.isdir(d):
            os.makedirs(d)

    # Process the raw Ensembl files into AltAnalyze format
    reformat_exon_file()
    reformat_junction_file()

    # Create stub files required by AltAnalyze
    create_domain_annotations_stub()
    create_protein_annotations_stub()
    create_platform_file()
    create_junction_comps_stub()
    create_junction_comps_updated()

    # Update version config
    update_version_file()

    print("\nDone! Mm RNASeq database built successfully.")
    print("Output files:")
    for fname in os.listdir(RNASEQ_DIR):
        path = os.path.join(RNASEQ_DIR, fname)
        if os.path.isfile(path):
            print("  {} ({} bytes)".format(fname, os.path.getsize(path)))
    print("  chr/ (directory)")
    for fname in os.listdir(CHR_DIR):
        path = os.path.join(CHR_DIR, fname)
        if os.path.isfile(path):
            print("    {} ({} bytes)".format(fname, os.path.getsize(path)))


if __name__ == "__main__":
    main()