rule get_reference_repeats:
  conda: "../env.yaml"
  input:
    mapping = "data/mapping.tsv"
  output: 
    temp("tmp/hs1.fa.out.gz"),
    temp("data/repeatsReferenceTE.bed"),
    temp("data/repeatsReferenceTE.bed.gz"),
    temp("data/repeatsReferenceTE.bed.gz.csi")
  shell:
    """
    wget -O {output[0]} https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz

    zcat {output[0]} | tail -n +4 | \
    grep -E "(S|L)INE|Retroposon|LTR|DNA" | \
    awk -v OFS="\\t" '{{print $5, $6, $7, $9, $10, $11}}' > tmp.repeats.bed

    awk 'BEGIN {{
      FS = OFS = "\\t";
      while ((getline < "{input.mapping}") > 0) {{
        map[$1] = $2
      }}
    }}
    {{
      if ($1 in map) {{
        $1 = map[$1]
        print
      }}
    }}' tmp.repeats.bed > {output[1]}

    bgzip -k {output[1]}
    tabix -p bed --csi {output[2]}
    """

rule prepare_te_nooverlap:
  conda: "../env.yaml"
  input:
    "data/repeatsReferenceTE.bed.gz"
  output:
    temp("tmp/repeatsReferenceTENoOverlap.bed")
  shell:
    """
    bedtools merge -c 4,5,6 -o collapse,collapse,collapse \
      -i {input[0]} > \
      {output[0]}
    """

rule get_reference_t2t_chm13:
  conda: "../env.yaml"
  output: 
    str(workflow.basedir) + "/data/GCA_009914755.4_T2T-CHM13v2.0_genomic.fa",
    "data/ref.fna.gz"
  shell: 
    """
    wget -O {output[1]} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
    gunzip {output[1]} > {output[0]}
    """