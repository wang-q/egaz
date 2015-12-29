# Self-alignment of S288c

## Detailed steps

### prepare seqs

```bash
mkdir -p ~/Scripts/egas/data
cp -R ~/data/alignment/example/scer/Genomes/S288c ~/Scripts/egas/data/
```

## self alignment

```bash
cd ~/Scripts/egas/data

if [ -d S288cvsselfalign ]
then
    rm -fr S288cvsselfalign
fi

perl ~/Scripts/egaz/bz.pl --parallel 8 \
    --is_self \
    -s set01 -C 0 --noaxt \
    -dt S288c \
    -dq S288c \
    -dl S288cvsselfalign

perl ~/Scripts/egaz/lpcna.pl --parallel 8 \
    -dt S288c \
    -dq S288c \
    -dl S288cvsselfalign
```

###  blast and merge

```bash
cd ~/Scripts/egas/data

if [ ! -d S288c_proc ]
then
    mkdir S288c_proc
fi

if [ ! -d S288c_result ]
then
    mkdir S288c_result
fi

cd ~/Scripts/egas/data/S288c_proc

# genome
find ../S288c -type f -name "*.fa" \
    | sort | xargs cat \
    | perl -nl -e '/^>/ or $_ = uc; print' \
    > genome.fa
faops size genome.fa > chr.sizes

# Correct genomic positions
fasops axt2fas ../S288cvsselfalign/axtNet/*.axt.gz -l 1000 -o stdout > axt.fas
fasops separate axt.fas --nodash -s .sep.fasta

perl ~/Scripts/egas/sparsemem_exact.pl -f target.sep.fasta -g genome.fa \
    --length 500 -o replace.target.tsv
fasops replace axt.fas replace.target.tsv -o axt.target.fas

perl ~/Scripts/egas/sparsemem_exact.pl -f query.sep.fasta -g genome.fa \
    --length 500 -o replace.query.tsv
fasops replace axt.target.fas replace.query.tsv -o axt.correct.fas

# coverage stats
fasops covers axt.correct.fas -o axt.correct.yml
runlist split axt.correct.yml -s .temp.yml
runlist compare --op union target.temp.yml query.temp.yml -o axt.union.yml
runlist stat --size chr.sizes axt.union.yml -o ../S288c_result/S288c.union.csv

# links by lastz-chain
fasops links axt.correct.fas -o stdout \
    | perl -nl -e 's/(target|query)\.//g; print;' \
    > links.lastz.tsv

# remove species names
fasops separate axt.correct.fas --nodash -o stdout \
    | perl -nl -e '/^>/ and s/^>(target|query)\./\>/; print;' \
    > axt.gl.fasta

# Get more paralogs
perl ~/Scripts/egas/blastn_genome.pl -c 0.95 -f axt.gl.fasta -g genome.fa -o axt.bg.fasta

if [ -e axt.bg.fasta ];
then
    cat axt.gl.fasta axt.bg.fasta > axt.all.fasta
else
    cat axt.gl.fasta > axt.all.fasta
fi

# link paralogs
perl ~/Scripts/egas/blastn_paralog.pl -f axt.all.fasta -c 0.95 -o links.blast.tsv

# merge
perl ~/Scripts/egas/merge_node.pl    -v -f links.lastz.tsv -f links.blast.tsv -o S288c.merge.yml -c 0.9
perl ~/Scripts/egas/paralog_graph.pl -v -f links.lastz.tsv -f links.blast.tsv -m S288c.merge.yml --nonself -o S288c.merge.graph.yml
perl ~/Scripts/egas/cc.pl               -f S288c.merge.graph.yml
perl ~/Scripts/egas/proc_cc_chop.pl     -f S288c.cc.yml --size chr.sizes --genome genome.fa --msa mafft
perl ~/Scripts/egas/proc_cc_stat.pl     -f S288c.cc.yml --size chr.sizes

runlist stat --size chr.sizes S288c.cc.chr.runlist.yml
perl ~/Scripts/egas/cover_figure.pl --size chr.sizes -f S288c.cc.chr.runlist.yml
```

### result & clean

```bash
cd ~/Scripts/egas/data/S288c_proc

cp S288c.cc.yml ../S288c_result
mv S288c.cc.csv ../S288c_result
mv S288c.cc.chr.runlist.yml.csv ../S288c_result/S288c.chr.csv
mv S288c.cc.chr.runlist.png ../S288c_result/S288c.chr.png

# clean
find . -type f -name "*genome.fa*" | xargs rm
find . -type f -name "*all.fasta*" | xargs rm
find . -type f -name "*.sep.fasta" | xargs rm
find . -type f -name "axt.*" | xargs rm
find . -type f -name "replace.*.tsv" | xargs rm
find . -type f -name "*.temp.yml" | xargs rm
```
