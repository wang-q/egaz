# Build alignments of several strains of *Saccharomyces cerevisiae*

As example for egaz and alignDB. Extract from Scer_wgs of
[`OPs-download.md`](https://github.com/wang-q/withncbi/blob/master/pop/OPs-download.md) and
[`OPs-align.md`](https://github.com/wang-q/withncbi/blob/master/pop/OPs-align.md).

## Download

1. Create `scer_example.tsv` manually.

    ```bash
    mkdir -p ~/data/alignment/example/scer      # operation directory
    mkdir -p ~/data/alignment/example/GENOMES   # sequence directory

    cd ~/data/alignment/example/GENOMES

    echo -e '#name\tprefix\torganism\tcontigs' > example_wgs.tsv
    echo -e "YJM789\tAAFW\tSaccharomyces cerevisiae YJM789\t258" >> example_wgs.tsv
    echo -e "Spar\tAABY\tSaccharomyces paradoxus NRRL Y-17217\t832" >> example_wgs.tsv
    ```

2. Create working directory and download WGS sequences.

    ```bash
    cd ~/data/alignment/example/GENOMES

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f example_wgs.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/example_wgs.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Download strains of *Saccharomyces cerevisiae* at good assembly status.

    ```bash
    mkdir -p ~/data/alignment/example/GENOMES/DOWNLOAD
    cd ~/data/alignment/example/GENOMES/DOWNLOAD

    # Download S288c and RM11_1a separately
    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000146045.2.assembly.txt \
        --nuclear -name S288c \
        > S288c.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000149365.1.assembly.txt \
        --genbank --scaffold -name RM11_1a \
        > RM11_1a.seq.csv

    echo "#strain_name,accession,strain_taxon_id,seq_name" > example.seq.csv
    cat S288c.seq.csv RM11_1a.seq.csv \
        | perl -nl -e '/^#/ and next; /^\s*$/ and next; print;' \
        >> example.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f example.seq.csv
    ```

## Align

In these steps, I don't use ensembl annotations. Just Assembly and RepeatMasker ones.

Stats on pairwise and multiple alignments are minimal.

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/example/scer
    cd ~/data/alignment/example/scer

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/example/GENOMES/WGS/example_wgs.data.yml \
        -o scer_test.yml \
        -d ~/data/alignment/example/GENOMES/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=scer \
        --opt base_dir='~/data/alignment/example' \
        --opt data_dir="~/data/alignment/example/scer" \
        --opt rm_species=Fungi \
        --dd ~/data/alignment/example/GENOMES/DOWNLOAD \
        --download "name=S288c;taxon=559292" \
        --download "name=RM11_1a;taxon=285006" \
        --plan 'name=Scer_n3_Spar;t=S288c;qs=RM11_1a,YJM789,Spar;o=Spar' \
        --plan 'name=Scer_n3_pop;t=S288c;qs=RM11_1a,YJM789' \
        --plan 'name=Scer_n2_Spar;t=S288c;qs=RM11_1a,Spar;o=Spar' \
        -y
    ```

2. Rest routing things.

    ```bash
    cd ~/data/alignment/example/scer

    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 8 -i scer_test.yml

    bash 01_file.sh
    bash 02_rm.sh
    bash 03_strain_info.sh

    # plan_ALL.sh
    bash plan_ALL.sh

    bash 1_real_chr.sh
    bash 3_pair_cmd.sh
    bash 4_rawphylo.sh
    bash 5_multi_cmd.sh

    # other plans
    bash plan_Scer_n3_Spar.sh
    bash 5_multi_cmd.sh

    bash plan_Scer_n3_pop.sh
    bash 5_multi_cmd.sh

    bash plan_Scer_n2_Spar.sh
    bash 5_multi_cmd.sh
    ```

## Raw alignments

```bash
mkdir -p ~/data/alignment/example/raw
cd ~/data/alignment/example/raw

# maf2fas
for pair in S288cvsRM11_1a S288cvsYJM789 S288cvsSpar ; do
    mkdir -p ~/data/alignment/example/raw/${pair};
    find ~/data/alignment/example/scer/Pairwise/${pair} -name "*.maf" -or -name "*.maf.gz" \
        | parallel --no-run-if-empty -j 4 \
            fasops maf2fas {} -o ~/data/alignment/example/raw/${pair}/{/}.fas
    fasops covers ~/data/alignment/example/raw/${pair}/*.fas -n S288c -l 1000 -t 10 -o ${pair}.yml
done

runlist compare --op intersect \
    S288cvsRM11_1a.yml \
    S288cvsYJM789.yml \
    S288cvsSpar.yml \
    -o intersect.raw.yml
runlist span --op excise -n 1000 intersect.raw.yml -o intersect.filter.yml
rm intersect.raw.yml

# slice
for pair in S288cvsRM11_1a S288cvsYJM789 S288cvsSpar ; do
    if [ -e ${pair}.slice.fas ];
    then
        rm ${pair}.slice.fas
    fi
    find ~/data/alignment/example/raw/${pair}/ -name "*.fas" -or -name "*.fas.gz" \
        | sort \
        | parallel --no-run-if-empty --keep-order -j 1 " \
            fasops slice {} intersect.filter.yml -n S288c -l 1000 -o stdout \
            >> ${pair}.slice.fas
        "
done

# join
fasops join \
    S288cvsRM11_1a.slice.fas \
    S288cvsYJM789.slice.fas \
    S288cvsSpar.slice.fas \
    -n S288c -o join1.fas

cat <<'EOF' > names.list
S288c
RM11_1a
YJM789
Spar
EOF

fasops subset join1.fas names.list --required -o join2.fas

# refine (can't use quick mode)
fasops refine join2.fas --msa mafft --parallel 8 -o join.fas
rm join1.fas join2.fas

# 12071326,9413812,0.7798
fasops covers -n S288c join.fas -l 1000 -o join.yml
runlist stat \
    --size ~/data/alignment/example/scer/Genomes/S288c/chr.sizes \
    --all -o stdout \
    join.yml

# 12071326,10356201,0.8579
fasops covers -n S288c ~/data/alignment/example/scer/plan_ALL_refined/*.gz -l 1000 -o mz.yml
runlist stat \
    --size ~/data/alignment/example/scer/Genomes/S288c/chr.sizes \
    --all -o stdout \
    mz.yml
```
