
/*
 * Prepare sequences
 */
process prepareMultifasta {
    input:
    path fasta

    output:
    path "prepared.fa"

    script:
    """
    seqkit seq -g ${fasta} | seqtk seq -l 10000000 - > prepared.fa
    """
}


/*
 * Split fasta into pairwise combinations
 */
process splitSequences {
    input:
    path fasta

    output:
    path "split/*.fa"

    script:
    """
    python3 - <<EOF
import os
from itertools import combinations

os.makedirs("split", exist_ok=True)
with open("${fasta}") as f:
    lines = f.read().splitlines()

seqs = {lines[i][1:]: lines[i+1] for i in range(0, len(lines), 2)}
for (a, b) in combinations(seqs, 2):
    with open(f"split/{a}-{b}.fa", "w") as out:
        out.write(f">{a}\\n{seqs[a]}\\n>{b}\\n{seqs[b]}\\n")
EOF
    """
}

// process computeSimilarity {

//     publishDir params.outdir
    
//     input:
//     path aligned_file

//     output:
//     path "similarity_scores/*.tsv"

//     script:
//     """
//     mkdir -p similarity_scores

//     python3 - <<EOF > similarity_scores/\$(basename ${aligned_file} .fa).tsv
// from Bio import SeqIO
// import os

// input_file = "${aligned_file}"
// pair_name = os.path.basename(input_file).replace(".fa", "")
// dist = 0
// gaps = 0
// recs = []

// try:
//     with open(input_file, "r") as hndlr:
//         for record in SeqIO.parse(hndlr, "fasta"):
//             recs.append(record)

//     if len(recs) != 2:
//         raise ValueError("Expected 2 sequences, got %d" % len(recs))

//     seq1 = str(recs[0].seq)
//     seq2 = str(recs[1].seq)

//     if len(seq1) != len(seq2):
//         raise ValueError("Sequences are not aligned")

//     for i in range(len(seq1)):
//         if seq1[i] != "-" and seq2[i] != "-":
//             if seq1[i] != seq2[i]:
//                 dist += 1
//         else:
//             gaps += 1
//     pair1, pair2 = pair_name.split("-")
//     similarity = 1 - (float(dist) / (len(seq1) - gaps))
//     alnlength = len(seq1)
//     query_cover = float((alnlength - gaps)/alnlength)

//     print(f"{pair1},{pair2},{similarity:.6f},{query_cover:.6f},{alnlength}")

// except Exception as e:
//     print(f"{pair1},{pair2},{similarity:.6f},{query_cover:.6f},{alnlength}")
// EOF
//     """
// }

// process mergeSimilarityScores {

//     publishDir params.outdir

//     input:
//     path similarity_files

//     output:
//     path 'similarity_scores.csv'

//     script:
//     """
//     echo "seq1,seq2,identity,query_cover,alignment_length" > similarity_scores.csv
//     cat ${similarity_files.join(' ')} >> similarity_scores.csv
//     """
// }

process computeAndMergeSimilarity {
    publishDir params.outdir, mode: 'copy'

    input:
    path aligned_files

    output:
    path "similarity_scores.csv"

    script:
    """
    python3 <<EOF
import os
import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool, cpu_count

def compute_similarity(filepath):
    try:
        pair_name = os.path.basename(filepath).replace(".fa", "")
        recs = []
        with open(filepath, "r") as hndlr:
            for record in SeqIO.parse(hndlr, "fasta"):
                recs.append(record)
        if len(recs) != 2:
            raise ValueError(f"Expected 2 sequences, got {len(recs)} in {filepath}")
        seq1 = str(recs[0].seq)
        seq2 = str(recs[1].seq)
        if len(seq1) != len(seq2):
            raise ValueError(f"Sequences are not aligned in {filepath}")
        dist = 0
        gaps = 0
        for i in range(len(seq1)):
            if seq1[i] != "-" and seq2[i] != "-":
                if seq1[i] != seq2[i]:
                    dist += 1
            else:
                gaps += 1
        pair1, pair2 = pair_name.split("-")
        similarity = 1 - (float(dist) / (len(seq1) - gaps))
        alnlength = len(seq1)
        query_cover = float((alnlength - gaps)/alnlength)
        return (pair1, pair2, similarity, query_cover, alnlength)
    except Exception as e:
        return (pair_name, "ERROR", 0, 0, 0)

def main():
    files = [f for f in os.listdir('.') if f.endswith('.fa')]
    files = [os.path.abspath(f) for f in files]
    with Pool(cpu_count()) as pool:
        results = pool.map(compute_similarity, files)
    df = pd.DataFrame(results, columns=["seq1", "seq2", "identity", "query_cover", "alignment_length"])
    df.to_csv("similarity_scores.csv", index=False)

if __name__ == "__main__":
    main()
EOF
    """
}

process cleanupAlignments {
    input:
    val dummy

    script:
    """
    rm -rf alignments
    """
}

process cleanupIntermediaries {
    input:
    val dummy

    script:
    """
    rm -rf adjusted* rotated*
    """
}