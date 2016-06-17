# Self-alignment of Arabidopsis thaliana

## Detailed steps

On server with dual E5-2690 v3.

### Prepare sequences

```bash
mkdir -p ~/data/self_alignment
cp -R ~/data/alignment/Ensembl/Atha ~/data/self_alignment
```

### self alignment

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
# real    21m18.036s
# user    149m57.173s
# sys     0m33.244s

time perl ~/Scripts/egaz/lpcna.pl --parallel 8 \
    -dt Atha \
    -dq Atha \
    -dl Athavsselfalign
# real    1m23.318s
# user    1m50.571s
# sys     1m12.166s

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
fasops axt2fas ../Athavsselfalign/axtNet/*.axt.gz -l 1000 -s chr.sizes -o stdout > axt.fas
fasops separate axt.fas --nodash -s .sep.fasta

time perl ~/Scripts/egas/sparsemem_exact.pl -f target.sep.fasta -g genome.fa \
    --length 500 -o replace.target.tsv
# real    0m56.162s
# user    1m29.831s
# sys     0m5.504s

fasops replace axt.fas replace.target.tsv -o axt.target.fas

perl ~/Scripts/egas/sparsemem_exact.pl -f query.sep.fasta -g genome.fa \
    --length 500 -o replace.query.tsv
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
fasops separate axt.correct.fas --nodash --rc -o stdout \
    | perl -nl -e '/^>/ and s/^>(target|query)\./\>/; print;' \
    | faops filter -u stdin stdout \
    | faops filter -n 250 stdin stdout \
    > axt.gl.fasta

# Get more paralogs
time perl ~/Scripts/egas/fasta_blastn.pl  -f axt.gl.fasta -g genome.fa -o axt.bg.blast
# real    0m30.664s
# user    0m57.902s
# sys     0m3.393s
time perl ~/Scripts/egas/blastn_genome.pl -f axt.bg.blast -g genome.fa -o axt.bg.fasta -c 0.95
# real    0m55.320s
# user    0m33.419s
# sys     0m24.626s

cat axt.gl.fasta axt.bg.fasta \
    | faops filter -u stdin stdout \
    | faops filter -n 250 stdin stdout \
    > axt.all.fasta

# link paralogs
time perl ~/Scripts/egas/fasta_blastn.pl   -f axt.all.fasta -g axt.all.fasta -o axt.all.blast
# real    0m19.575s
# user    0m55.607s
# sys     0m1.503s
perl ~/Scripts/egas/blastn_paralog.pl -f axt.all.blast -c 0.95 -o links.blast.tsv

# merge
time perl ~/Scripts/egas/merge_node.pl    -v -f links.lastz.tsv -f links.blast.tsv -o Atha.merge.yml -c 0.95
# real    3m48.494s
# user    10m37.467s
# sys     0m1.047s
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
