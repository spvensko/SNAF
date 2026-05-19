#!/usr/bin/env python3
"""
Build mouse (GRCm39/mm39) reference databases for mSNAF and AltAnalyze from GENCODE vM31.

Key differences from the GRCm38/Ensembl build:
  - GTF source: GENCODE vM31 (not Ensembl release 100)
  - GENCODE GTF uses chr-prefixed chromosomes (chr1, chrX, etc.)
  - GENCODE gene/transcript IDs have version suffixes (ENSMUSG00000064842.3)
  - Assembly: GRCm39 (not GRCm38)

Output directory structure:
  Mm_Alt31_db/
    Mm_Ensembl_exon_add_col.txt
    Mm_mRNA-ExonIDs.txt
    Mm_gene-seq-2000_flank.fa
    Mus_musculus.GRCm39.gtf
    Mm_df_start_codon.txt
  AltDatabase/EnsMart31/ensembl/Mm/
    Mm_Ensembl_exon.txt
    Mm_Ensembl_junction.txt
    Mm_Ensembl-annotations.txt
    Mm_Ensembl_transcript-annotations.txt
    Mm_Ensembl_transcript-biotypes.txt
    Mm.bed

Usage:
    python3 build_mm39_reference_db.py --output-dir /path/to/output [--skip-fasta]
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
import urllib.request
from collections import defaultdict

# GENCODE vM31 configuration
GENCODE_RELEASE = "M31"
ENSEMBL_RELEASE = 108  # Corresponding Ensembl release for EnsMart version
ASSEMBLY = "GRCm39"
GENE_PREFIX = "ENSMUSG"
TRANSCRIPT_PREFIX = "ENSMUST"
PROTEIN_PREFIX = "ENSMUSP"
EXON_PREFIX = "ENSE"
DB_DIR = "Mm_Alt31_db"
ENSMART_VERSION = "EnsMart31"

GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.primary_assembly.annotation.gtf.gz"
FASTA_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/GRCm39.primary_assembly.genome.fa.gz"

# AltAnalyze expects chr prefix in chromosome names to match BAM files
STRIP_CHR_PREFIX = False


def download_file(url: str, dest: str) -> None:
    if os.path.exists(dest):
        print(f"  [skip] {dest} already exists")
        return
    print(f"  [download] {url}")
    urllib.request.urlretrieve(url, dest)


def strip_id_version(gene_id: str) -> str:
    """Strip version suffix from GENCODE gene/transcript IDs.
    ENSMUSG00000064842.3 -> ENSMUSG00000064842
    """
    if "." in gene_id:
        return gene_id.split(".")[0]
    return gene_id


def normalize_chrom(chrom: str) -> str:
    """Strip chr prefix if configured.
    chr1 -> 1, chrX -> X, chrMT -> MT
    """
    if STRIP_CHR_PREFIX and chrom.startswith("chr"):
        return chrom[3:]
    return chrom


def parse_gencode_gtf(gtf_path: str) -> tuple[dict, dict, list, dict]:
    """Parse GENCODE GTF and return gene, transcript, exon, and start_codon records.

    Handles GENCODE-specific format:
    - Strips version suffixes from gene/transcript IDs
    - Normalizes chr-prefixed chromosome names
    - Uses gene_type/transcript_type instead of gene_biotype/transcript_biotype
    """
    genes = {}
    transcripts = {}
    exons = []
    start_codons = defaultdict(list)

    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            start, end = int(start), int(end)
            chrom = normalize_chrom(chrom)

            attr_dict = {}
            for attr in attrs.split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                key, _, value = attr.partition(' ')
                value = value.strip('"')
                attr_dict[key] = value

            if feature == 'gene':
                ensg = strip_id_version(attr_dict.get('gene_id', ''))
                gene_name = attr_dict.get('gene_name', '')
                genes[ensg] = (chrom, start, end, strand, gene_name)

            elif feature == 'transcript':
                ensg = strip_id_version(attr_dict.get('gene_id', ''))
                enst = strip_id_version(attr_dict.get('transcript_id', ''))
                biotype = attr_dict.get('transcript_type', attr_dict.get('gene_type', ''))
                transcripts[enst] = (ensg, chrom, start, end, strand, biotype)

            elif feature == 'exon':
                ensg = strip_id_version(attr_dict.get('gene_id', ''))
                enst = strip_id_version(attr_dict.get('transcript_id', ''))
                exon_number = attr_dict.get('exon_number', '')
                exon_id = strip_id_version(attr_dict.get('exon_id', ''))
                try:
                    exon_number = int(exon_number)
                except (ValueError, TypeError):
                    exon_number = 0
                exons.append((ensg, enst, chrom, start, end, strand, exon_number, exon_id))

            elif feature == 'start_codon':
                ensg = strip_id_version(attr_dict.get('gene_id', ''))
                start_codons[ensg].append(start if strand == '+' else end)

    return genes, transcripts, exons, dict(start_codons)


def assign_subexon_ids(genes, transcripts, exons_raw):
    """Assign AltAnalyze-style sub-exon IDs to exons."""
    # Group exons by gene, then by transcript
    gene_exons = defaultdict(list)
    gene_transcript_exons = defaultdict(lambda: defaultdict(list))

    for ensg, enst, chrom, start, end, strand, exon_number, exon_id in exons_raw:
        gene_exons[ensg].append((start, end, strand, enst, exon_number, exon_id))
        gene_transcript_exons[ensg][enst].append((start, end, exon_number, exon_id))

    exon_table_rows = []
    junction_rows = []

    for ensg, exon_list in gene_exons.items():
        if ensg not in genes:
            continue
        chrom, _, _, strand, gene_name = genes[ensg]

        # Collect all unique splice boundaries for this gene
        boundary_groups = defaultdict(list)
        all_boundaries = set()
        for start, end, s, enst, enum, eid in exon_list:
            all_boundaries.add(start)
            all_boundaries.add(end)

        # Sort boundaries
        sorted_boundaries = sorted(all_boundaries)

        # Assign exon group numbers based on boundaries
        # Exons sharing a boundary (donor or acceptor site) are in the same group
        exon_to_group = {}
        exon_group_counter = 1

        # Sort exons by position
        sorted_exons = sorted(exon_list, key=lambda x: (x[0], x[1]))

        # Group by (start, end) to find unique exon coordinates
        unique_coords = sorted(set((s, e) for s, e, _, _, _, _ in exon_list),
                               key=lambda x: (x[0], x[1]))

        # Map each unique coordinate to an exon group number
        # Group by shared splice sites
        coord_groups = {}
        current_group = 1
        for i, (start, end) in enumerate(unique_coords):
            # Check if this coordinate shares a boundary with the previous one
            if i > 0:
                prev_start, prev_end = unique_coords[i - 1]
                if start != prev_end:  # New exon group
                    current_group += 1
            coord_groups[(start, end)] = current_group

        # Assign sub-exon IDs within each group
        coord_subexon = {}
        group_subexon_counter = defaultdict(int)
        for start, end in unique_coords:
            group = coord_groups[(start, end)]
            group_subexon_counter[group] += 1
            subexon_idx = group_subexon_counter[group]
            coord_subexon[(start, end)] = f"E{group}.{subexon_idx}"

        # Collect Ensembl exon IDs per unique coordinate (pipe-separated)
        coord_ens_ids = {}
        coord_id_sets = defaultdict(set)
        for start, end, s, enst, enum, eid in exon_list:
            if eid:
                coord_id_sets[(start, end)].add(eid)
        for coord, id_set in coord_id_sets.items():
            sorted_ids = sorted(id_set)
            coord_ens_ids[coord] = "|".join(sorted_ids)

        # Now build the exon table
        for start, end, s, enst, enum, eid in sorted_exons:
            exon_id_str = coord_subexon.get((start, end), f"E0.1")
            constitutive = "yes" if len(set((s, e) for s, e, _, _, _, _ in exon_list)) == 1 else "no"
            exon_table_rows.append((ensg, exon_id_str, chrom, strand, start, end, constitutive, eid, "", "", "False"))

        # Generate junction records: pairs of adjacent exons connected by an intron
        # Format: gene  E1.1-E2.1  chrom  strand  start1|start2  stop1|stop2  constitutive  eid1|eid2
        for i in range(len(unique_coords) - 1):
            end_curr = unique_coords[i][1]
            start_next = unique_coords[i + 1][0]
            if end_curr < start_next:
                exon1_id = coord_subexon[unique_coords[i]]
                exon2_id = coord_subexon[unique_coords[i + 1]]
                junc_id = f"{exon1_id}-{exon2_id}"
                starts = f"{unique_coords[i][0]}|{unique_coords[i + 1][0]}"
                stops = f"{unique_coords[i][1]}|{unique_coords[i + 1][1]}"
                ens_ids_1 = coord_ens_ids.get(unique_coords[i], "")
                ens_ids_2 = coord_ens_ids.get(unique_coords[i + 1], "")
                if ens_ids_1 and ens_ids_2:
                    ens_ids = f"{ens_ids_1}|{ens_ids_2}"
                elif ens_ids_1:
                    ens_ids = ens_ids_1
                elif ens_ids_2:
                    ens_ids = ens_ids_2
                else:
                    ens_ids = ""
                junction_rows.append((ensg, junc_id, chrom, strand, starts, stops, "no", ens_ids, "", ""))

    # Build gene_transcript_exons with new exon IDs
    gene_transcript_map = {}
    for ensg in gene_transcript_exons:
        gene_transcript_map[ensg] = {}
        for enst, exon_list in gene_transcript_exons[ensg].items():
            exon_ids = []
            for start, end, enum, eid in sorted(exon_list, key=lambda x: (x[0], x[1])):
                if (start, end) in coord_subexon:
                    exon_ids.append(coord_subexon[(start, end)])
            gene_transcript_map[ensg][enst] = " ".join(exon_ids) if exon_ids else ""

    return exon_table_rows, gene_transcript_map, junction_rows


def build_exon_table(exon_table_rows):
    """Build the exon table with is_suffer column."""
    lines = ["gene\texon-id\tchromosome\tstrand\texon-region-start(s)\texon-region-stop(s)\tconstitutive_call\tens_exon_ids\tsplice_events\tsplice_junctions\tis_suffer"]
    for row in exon_table_rows:
        gene, exon_id, chrom, strand, start, end, constitutive, eid, se, sj, is_suffer = row
        lines.append(f"{gene}\t{exon_id}\t{chrom}\t{strand}\t{start}\t{end}\t{constitutive}\t{eid}\t{se}\t{sj}\t{is_suffer}")
    return "\n".join(lines) + "\n"


def build_transcript_db(gene_transcript_map, genes):
    """Build the transcript-to-exon mapping file."""
    lines = ["EnsGID\tEnsTID\tEnsPID\tExons"]
    for ensg in sorted(gene_transcript_map.keys()):
        for enst in sorted(gene_transcript_map[ensg].keys()):
            ensp = enst.replace("ENSMUST", "ENSMUSP")
            exons = gene_transcript_map[ensg][enst]
            lines.append(f"{ensg}\t{enst}\t{ensp}\t{exons}")
    return "\n".join(lines) + "\n"


def build_start_codon_table(start_codons, genes):
    """Build start codon position table."""
    lines = ["start_codon\tnon_redundant"]
    for ensg in sorted(start_codons.keys()):
        positions = sorted(start_codons[ensg])
        lines.append(f"{ensg}\t{positions}")
    return "\n".join(lines) + "\n"


def build_annotations(genes, transcripts, gene_transcript_map):
    """Build gene annotation files for AltAnalyze.

    Mm_Ensembl-annotations.txt: 4 columns (gene_id, description, gene_name, rna_processing)
    Mm_Ensembl-annotations_simple.txt: 3 columns (gene_id, description, gene_name)
    AltAnalyze reads _simple.txt via importGeneAnnotations() during expression building,
    then writes _annotations.txt via exportEnsemblAnnotations() with RNA_processing added.
    """
    # Mm_Ensembl-annotations.txt: gene_id, description, gene_name, rna_processing
    ann_lines = ["gene_id\tdescription\tgene_name\trna_processing"]
    for ensg in sorted(genes.keys()):
        chrom, start, end, strand, gene_name = genes[ensg]
        description = gene_name  # GENCODE GTF has no separate description field
        ann_lines.append(f"{ensg}\t{description}\t{gene_name}\t")

    # Mm_Ensembl-annotations_simple.txt: gene_id, description, gene_name
    ann_simple_lines = ["Ensembl Gene ID\tDescription\tGene name"]
    for ensg in sorted(genes.keys()):
        chrom, start, end, strand, gene_name = genes[ensg]
        description = gene_name
        ann_simple_lines.append(f"{ensg}\t{description}\t{gene_name}")

    # Mm_Ensembl_transcript-annotations.txt: gene_id, transcript_id, transcript_name, exon_ids
    trans_ann_lines = ["gene_id\ttranscript_id\ttranscript_name\texon_ids"]
    for ensg in sorted(gene_transcript_map.keys()):
        for enst, exon_ids in sorted(gene_transcript_map[ensg].items()):
            trans_ann_lines.append(f"{ensg}\t{enst}\t{enst}\t{exon_ids}")

    # Mm_Ensembl_transcript-biotypes.txt: gene_id, transcript_id, biotype
    bio_lines = ["gene_id\ttranscript_id\tbiotype"]
    for enst in sorted(transcripts.keys()):
        ensg, chrom, start, end, strand, biotype = transcripts[enst]
        bio_lines.append(f"{ensg}\t{enst}\t{biotype}")

    return "\n".join(ann_lines) + "\n", "\n".join(ann_simple_lines) + "\n", "\n".join(trans_ann_lines) + "\n", "\n".join(bio_lines) + "\n"


def build_altanalyze_exon_table(exon_table_rows):
    """Build the 10-column exon table (without is_suffer) for AltAnalyze."""
    lines = ["gene\texon-id\tchromosome\tstrand\texon-region-start(s)\texon-region-stop(s)\tconstitutive_call\tens_exon_ids\tsplice_events\tsplice_junctions"]
    for row in exon_table_rows:
        gene, exon_id, chrom, strand, start, end, constitutive, eid, se, sj, _ = row
        lines.append(f"{gene}\t{exon_id}\t{chrom}\t{strand}\t{start}\t{end}\t{constitutive}\t{eid}\t{se}\t{sj}")
    return "\n".join(lines) + "\n"


def build_junction_table(junction_rows):
    """Build exon-exon junction table for AltAnalyze.

    Each row represents a junction between two adjacent exons with format:
      gene  exon1_id-exon2_id  chrom  strand  start1|start2  stop1|stop2  constitutive  eid1|eid2  splice_events  splice_junctions
    """
    lines = ["gene\texon-id\tchromosome\tstrand\texon-region-start(s)\texon-region-stop(s)\tconstitutive_call\tens_exon_ids\tsplice_events\tsplice_junctions"]
    for row in junction_rows:
        gene, junc_id, chrom, strand, starts, stops, constitutive, eids, se, sj = row
        lines.append(f"{gene}\t{junc_id}\t{chrom}\t{strand}\t{starts}\t{stops}\t{constitutive}\t{eids}\t{se}\t{sj}")
    return "\n".join(lines) + "\n"


def build_bed_file(exon_table_rows):
    """Build BED file for AltAnalyze junction calling."""
    lines = []
    # Skip header
    for row in exon_table_rows:
        gene, exon_id, chrom, strand, start, end, constitutive, eid, se, sj, _ = row
        # BED is 0-based half-open
        bed_start = start - 1  # Convert from 1-based to 0-based
        bed_end = end
        name = f"{gene}:{exon_id}"
        lines.append(f"{chrom}\t{bed_start}\t{bed_end}\t{name}\t0\t{strand}")
    return "\n".join(lines) + "\n"


def index_fasta(fasta_path):
    """Build a byte-offset index for a FASTA file (like samtools faidx).

    Returns dict: chrom -> (file_offset, seq_length, line_bps, line_bytes)
    """
    idx = {}
    with open(fasta_path, 'rb') as f:
        offset = 0
        current_chrom = None
        current_offset = 0
        first_data_line_len = None
        line_bps = None
        line_bytes = None
        seq_length = 0
        for raw_line in f:
            line = raw_line.rstrip(b'\n\r')
            if line.startswith(b'>'):
                if current_chrom is not None:
                    idx[current_chrom] = (current_offset, seq_length, line_bps, line_bytes)
                current_chrom = line[1:].split()[0].decode()
                current_offset = f.tell()
                seq_length = 0
                first_data_line_len = None
                line_bps = None
                line_bytes = None
            else:
                if first_data_line_len is None and len(line) > 0:
                    first_data_line_len = len(raw_line)
                    line_bps = len(line)
                    line_bytes = first_data_line_len
                seq_length += len(line)
        if current_chrom is not None:
            idx[current_chrom] = (current_offset, seq_length, line_bps, line_bytes)
    return idx


COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def fetch_seq(fasta_fh, idx, chrom, start_0based, length):
    """Extract a subsequence from an indexed FASTA file.

    start_0based: 0-based start position
    length: number of bases to extract
    """
    if chrom not in idx:
        return None
    file_offset, seq_length, line_bps, line_bytes = idx[chrom]
    if start_0based + length > seq_length:
        length = seq_length - start_0based
    if length <= 0:
        return None

    # Calculate file position using line structure
    lines_before = start_0based // line_bps
    offset_in_line = start_0based % line_bps
    seek_pos = file_offset + lines_before * line_bytes + offset_in_line

    fasta_fh.seek(seek_pos)
    seq_parts = []
    remaining = length
    while remaining > 0:
        chunk = fasta_fh.read(min(remaining + 100, line_bytes + 1))
        if not chunk:
            break
        clean = chunk.replace(b'\n', b'').replace(b'\r', b'')
        if len(clean) > remaining:
            clean = clean[:remaining]
        seq_parts.append(clean.decode())
        remaining -= len(clean)

    return ''.join(seq_parts)


def build_gene_fasta(gtf_path, fasta_path, output_path, genes):
    """Build gene-seq-2000_flank.fa from genome FASTA and GTF coordinates.

    Each gene gets a FASTA entry with 2000bp flanking regions:
    >GENE_ID|chromosome|start-2000|end+2000
    """
    FLANK = 2000
    print(f"Indexing genome FASTA {fasta_path}...")
    idx = index_fasta(fasta_path)
    print(f"  Indexed {len(idx)} chromosomes")

    print(f"Building gene FASTA...")
    written = 0
    skipped = 0
    with open(fasta_path, 'rb') as fasta_fh, open(output_path, 'w') as out:
        for ensg in sorted(genes.keys()):
            chrom, start, end, strand, gene_name = genes[ensg]
            # Try both with and without chr prefix
            chrom_fa = chrom if chrom in idx else f"chr{chrom}"
            if chrom_fa not in idx:
                skipped += 1
                continue
            seq_length = idx[chrom_fa][1]
            flank_start = max(1, start - FLANK)
            flank_end = min(seq_length, end + FLANK)
            try:
                seq = fetch_seq(fasta_fh, idx, chrom_fa, flank_start - 1, flank_end - flank_start + 1)
                if seq is None:
                    skipped += 1
                    continue
                if strand == '-':
                    seq = seq.translate(COMPLEMENT)[::-1]
                out.write(f">{ensg}|{chrom}|{flank_start}|{flank_end}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i+80] + "\n")
                written += 1
            except Exception as e:
                skipped += 1
                if skipped <= 5:
                    print(f"  [warn] skipping {ensg}: {e}")

    print(f"Wrote {written} gene sequences, skipped {skipped}")


def main():
    parser = argparse.ArgumentParser(description="Build GRCm39/mm39 reference databases from GENCODE vM31")
    parser.add_argument("--output-dir", required=True, help="Output directory for reference files")
    parser.add_argument("--gtf", help="Pre-downloaded GENCODE GTF file (.gtf or .gtf.gz)")
    parser.add_argument("--fasta", help="Pre-downloaded genome FASTA file (.fa or .fa.gz)")
    parser.add_argument("--skip-fasta", action="store_true", help="Skip gene FASTA generation (very slow)")
    args = parser.parse_args()

    db_dir = os.path.join(args.output_dir, DB_DIR)
    alt_db_dir = os.path.join(args.output_dir, "AltDatabase", ENSMART_VERSION, "ensembl", "Mm")
    os.makedirs(db_dir, exist_ok=True)
    os.makedirs(alt_db_dir, exist_ok=True)

    # Download GTF
    gtf_path = args.gtf
    if not gtf_path:
        gtf_path = os.path.join(args.output_dir, "gencode.vM31.primary_assembly.annotation.gtf.gz")
        download_file(GTF_URL, gtf_path)

    # Download genome FASTA (skip entirely if --skip-fasta, saves ~800MB)
    fasta_path = args.fasta
    if not args.skip_fasta and not fasta_path:
        fasta_path = os.path.join(args.output_dir, "GRCm39.primary_assembly.genome.fa.gz")
        download_file(FASTA_URL, fasta_path)

    # Parse GTF
    print("Parsing GENCODE GTF...")
    genes, transcripts, exons_raw, start_codons = parse_gencode_gtf(gtf_path)
    print(f"  Found {len(genes)} genes, {len(transcripts)} transcripts, {len(exons_raw)} exons")

    # Assign subexon IDs
    print("Assigning subexon IDs...")
    exon_table_rows, gene_transcript_map, junction_rows = assign_subexon_ids(genes, transcripts, exons_raw)
    print(f"  Generated {len(exon_table_rows)} exon entries, {len(junction_rows)} junction entries")

    # Write mSNAF reference files
    print(f"Writing mSNAF reference files to {db_dir}...")

    # 1. Exon table with is_suffer
    exon_table = build_exon_table(exon_table_rows)
    with open(os.path.join(db_dir, "Mm_Ensembl_exon_add_col.txt"), 'w') as f:
        f.write(exon_table)
    print(f"  Wrote Mm_Ensembl_exon_add_col.txt ({len(exon_table_rows)} entries)")

    # 2. Transcript-to-exon mapping
    transcript_db = build_transcript_db(gene_transcript_map, genes)
    with open(os.path.join(db_dir, "Mm_mRNA-ExonIDs.txt"), 'w') as f:
        f.write(transcript_db)
    print(f"  Wrote Mm_mRNA-ExonIDs.txt")

    # 3. Start codon table
    start_codon_table = build_start_codon_table(start_codons, genes)
    with open(os.path.join(db_dir, "Mm_df_start_codon.txt"), 'w') as f:
        f.write(start_codon_table)
    print(f"  Wrote Mm_df_start_codon.txt ({len(start_codons)} genes)")

    # 4. Copy GTF
    gtf_dest = os.path.join(db_dir, f"Mus_musculus.{ASSEMBLY}.gtf")
    if gtf_path.endswith('.gz'):
        import shutil
        with gzip.open(gtf_path, 'rb') as f_in:
            with open(gtf_dest, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        import shutil
        shutil.copy2(gtf_path, gtf_dest)
    print(f"  Wrote {gtf_dest}")

    # 5. Gene FASTA (slow)
    if not args.skip_fasta:
        fasta_output = os.path.join(db_dir, "Mm_gene-seq-2000_flank.fa")
        # Decompress FASTA if needed
        decompressed_fasta = fasta_path
        if fasta_path.endswith('.gz'):
            decompressed_fasta = fasta_path.replace('.gz', '')
            if not os.path.exists(decompressed_fasta):
                print(f"  Decompressing genome FASTA...")
                import gzip as gz
                import shutil
                with gz.open(fasta_path, 'rb') as f_in:
                    with open(decompressed_fasta, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            fasta_path_for_build = decompressed_fasta
        else:
            fasta_path_for_build = fasta_path

        build_gene_fasta(gtf_path, fasta_path_for_build, fasta_output, genes)
    else:
        print("  Skipping gene FASTA (--skip-fasta)")

    # Write AltAnalyze database files
    print(f"Writing AltAnalyze database files to {alt_db_dir}...")

    # 1. Exon table (10-col, no is_suffer)
    alt_exon_table = build_altanalyze_exon_table(exon_table_rows)
    with open(os.path.join(alt_db_dir, "Mm_Ensembl_exon.txt"), 'w') as f:
        f.write(alt_exon_table)

    # 2. Junction table
    junction_table = build_junction_table(junction_rows)
    with open(os.path.join(alt_db_dir, "Mm_Ensembl_junction.txt"), 'w') as f:
        f.write(junction_table)

    # 3-6. Annotation files
    ann, ann_simple, trans_ann, bio = build_annotations(genes, transcripts, gene_transcript_map)
    with open(os.path.join(alt_db_dir, "Mm_Ensembl-annotations.txt"), 'w') as f:
        f.write(ann)
    with open(os.path.join(alt_db_dir, "Mm_Ensembl-annotations_simple.txt"), 'w') as f:
        f.write(ann_simple)
    with open(os.path.join(alt_db_dir, "Mm_Ensembl_transcript-annotations.txt"), 'w') as f:
        f.write(trans_ann)
    with open(os.path.join(alt_db_dir, "Mm_Ensembl_transcript-biotypes.txt"), 'w') as f:
        f.write(bio)

    # 6. BED file
    bed = build_bed_file(exon_table_rows)
    with open(os.path.join(alt_db_dir, "Mm.bed"), 'w') as f:
        f.write(bed)

    print(f"\nDone! Reference files written to {args.output_dir}")
    print(f"  mSNAF files: {db_dir}")
    print(f"  AltAnalyze files: {alt_db_dir}")


if __name__ == "__main__":
    main()