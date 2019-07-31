## Chunyu Zhao 20190319

import Bio.SeqIO, pysam, numpy as np
from collections import defaultdict
import argparse

def keep_read(aln):
    ## MIDAS default: discard reads
    global aln_stats
    aln_stats['aligned_reads'] += 1
    min_mapid = 94.0
    min_mapq = 20
    min_reafq = 20
    min_aln_cov = 0.75
    align_len = len(aln.query_alignment_sequence)
    query_len = aln.query_length
    if 100*(align_len-dict(aln.tags)['NM'])/float(align_len) < min_mapid:
        return False
    if np.mean(aln.query_qualities) < min_reafq:
        return False
    if aln.mapping_quality < min_mapq:
        return False
    if align_len/float(query_len) < min_aln_cov:
        return False
    aln_stats['mapped_reads'] += 1
    return True


def mpileup(genome_fp, bampath, pileup_file, stats_file):
    contigs_dict = {}
    for rec in Bio.SeqIO.parse(genome_fp, 'fasta'):
        contigs_dict[rec.id] = str(rec.seq)
    
    with open(pileup_file, "w") as out_file:
        with open(stats_file, "w") as fstat:
            header = ['ref_id', 'ref_pos', 'ref_allele', 'depth', 'count_a', 'count_c', 'count_g', 'count_t']
            out_file.write('\t'.join(header)+'\n')
            header2 = ['ref_id', 'genome_length', 'total_depth', 'covered_bases', 'aligned_reads', 'mapped_reads']
            fstat.write('\t'.join(header2)+'\n')

            with pysam.AlignmentFile(bampath, 'rb') as bamfile:
                for contig_id in sorted(list(contigs_dict.keys())):

                    global aln_stats
                    aln_stats = {'ref_id': contig_id, 'genome_length':0,'total_depth':0,'covered_bases':0,'aligned_reads':0,'mapped_reads':0}

                    contig_length = len(contigs_dict[contig_id])
                    counts = bamfile.count_coverage(
                                        contig_id,
                                        start=0,
                                        end=contig_length,
                                        quality_threshold=30, ## baseq
                                        read_callback=keep_read)
                    ## four array.arrays of the same length in order A C G T

                    for i in range(0, contig_length):
                        ref_pos = i + 1
                        ref_allele = contigs_dict[contig_id][i].upper()
                        depth = sum([counts[_][i] for _ in range(4)])
                        count_a = counts[0][i]
                        count_c = counts[1][i]
                        count_g = counts[2][i]
                        count_t = counts[3][i]
                        if depth > 0:
                            row = [contig_id, ref_pos, ref_allele, depth, count_a, count_c, count_g, count_t]
                            out_file.write('\t'.join([str(_) for _ in row])+'\n')
                            aln_stats['covered_bases'] += 1
                        aln_stats['genome_length'] += 1
                        aln_stats['total_depth'] += depth
                    ## write to file
                    fstat.write('\t'.join([str(_) for _ in aln_stats.values()])+'\n')

def main():
    p = argparse.ArgumentParser(prog="python midas_mpileup.py",
		 description='Do mpileup similar to MIDAS approache.')

    p.add_argument(
        "--reference-genome", required=True,
        type=argparse.FileType("r"),
        help="Reference genome used for BOWTIE2 alignment")
    p.add_argument(
        "--bam-file", required=True,
        type=argparse.FileType("r"),
        help="BAM alignment file")

    p.add_argument(
        "--pileup-file", required=True,
        help="Output pileup file")
    p.add_argument(
        "--stats-file", required=True,
        help="Output genome stats file")

    args = p.parse_args()
    mpileup(args.reference_genome, args.bam_file, args.pileup_file, args.stats_file)


main()
