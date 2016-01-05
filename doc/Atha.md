# Self-alignment of Arabidopsis thaliana

## Detailed steps

### prepare seqs

```bash
mkdir -p ~/data/self_alignment
cp -R ~/data/alignment/Ensembl/Atha ~/data/self_alignment
```

## self alignment

```bash
cd ~/data/self_alignment

if [ -d Athavsselfalign ]
then
    rm -fr Athavsselfalign
fi

time perl ~/Scripts/egaz/bz.pl --parallel 8 \
    --is_self \
    -s set01 -C 0 --noaxt \
    -dt Atha \
    -dq Atha \
    -dl Athavsselfalign
# real    27m54.586s
# user    187m1.586s
# sys     0m11.192s

time perl ~/Scripts/egaz/lpcna.pl --parallel 8 \
    -dt Atha \
    -dq Atha \
    -dl Athavsselfalign
# real    1m31.218s
# user    1m49.880s
# sys     1m7.008s

```

###  blast and merge

```bash
cd ~/data/self_alignment

if [ ! -d Atha_proc ]
then
    mkdir -p Atha_proc
fi

if [ ! -d Atha_result ]
then
    mkdir -p Atha_result
fi

cd ~/data/self_alignment/Atha_proc

# genome
find ../Atha -type f -name "*.fa" \
    | sort | xargs cat \
    | perl -nl -e '/^>/ or $_ = uc; print' \
    > genome.fa
faops size genome.fa > chr.sizes

# Correct genomic positions
fasops axt2fas ../Athavsselfalign/axtNet/*.axt.gz -l 1000 -o stdout > axt.fas
fasops separate axt.fas --nodash -s .sep.fasta

time perl ~/Scripts/egas/sparsemem_exact.pl -f target.sep.fasta -g genome.fa \
    --length 500 -o replace.target.tsv
# real    1m18.015s
# user    2m15.346s
# sys     0m4.220s

fasops replace axt.fas replace.target.tsv -o axt.target.fas

time perl ~/Scripts/egas/sparsemem_exact.pl -f query.sep.fasta -g genome.fa \
    --length 500 -o replace.query.tsv
# real    1m16.346s
# user    2m11.253s
# sys     0m5.515s
fasops replace axt.target.fas replace.query.tsv -o axt.correct.fas

# coverage stats
fasops covers axt.correct.fas -o axt.correct.yml
runlist split axt.correct.yml -s .temp.yml
runlist compare --op union target.temp.yml query.temp.yml -o axt.union.yml
runlist stat --size chr.sizes axt.union.yml -o ../Atha_result/Atha.union.csv

# links by lastz-chain
fasops links axt.correct.fas -o stdout \
    | perl -nl -e 's/(target|query)\.//g; print;' \
    > links.lastz.tsv

# remove species names
fasops separate axt.correct.fas --nodash -o stdout \
    | perl -nl -e '/^>/ and s/^>(target|query)\./\>/; print;' \
    > axt.gl.fasta

# Get more paralogs
time perl ~/Scripts/egas/blastn_genome.pl -c 0.95 -f axt.gl.fasta -g genome.fa -o axt.bg.fasta
# real    8m46.682s
# user    19m46.518s
# sys     0m20.860s

if [ -e axt.bg.fasta ];
then
    cat axt.gl.fasta axt.bg.fasta > axt.all.fasta
else
    cat axt.gl.fasta > axt.all.fasta
fi

# link paralogs
time perl ~/Scripts/egas/blastn_paralog.pl -f axt.all.fasta -c 0.95 -o links.blast.tsv
# real    1m21.303s
# user    3m46.030s
# sys     0m1.738s

# merge
time perl ~/Scripts/egas/merge_node.pl    -v -f links.lastz.tsv -f links.blast.tsv -o Atha.merge.yml -c 0.95
# real    7m44.699s
# user    22m14.711s
# sys     0m0.865s
perl ~/Scripts/egas/paralog_graph.pl -v -f links.lastz.tsv -f links.blast.tsv -m Atha.merge.yml --nonself -o Atha.merge.graph.yml
perl ~/Scripts/egas/cc.pl               -f Atha.merge.graph.yml
time perl ~/Scripts/egas/proc_cc_chop.pl     -f Atha.cc.raw.yml --size chr.sizes --genome genome.fa --msa mafft
# real    25m31.983s
# user    13m42.967s
# sys     15m1.670s
perl ~/Scripts/egas/proc_cc_stat.pl     -f Atha.cc.yml --size chr.sizes

runlist stat --size chr.sizes Atha.cc.chr.runlist.yml
perl ~/Scripts/egas/cover_figure.pl --size chr.sizes -f Atha.cc.chr.runlist.yml
```

### result & clean

```bash
cd ~/Scripts/egas/data/Atha_proc

cp Atha.cc.yml ../Atha_result
mv Atha.cc.csv ../Atha_result
mv Atha.cc.chr.runlist.yml.csv ../Atha_result/Atha.chr.csv
mv Atha.cc.chr.runlist.png ../Atha_result/Atha.chr.png

# clean
find . -type f -name "*genome.fa*" | xargs rm
find . -type f -name "*all.fasta*" | xargs rm
find . -type f -name "*.sep.fasta" | xargs rm
find . -type f -name "axt.*" | xargs rm
find . -type f -name "replace.*.tsv" | xargs rm
find . -type f -name "*.temp.yml" | xargs rm
```
