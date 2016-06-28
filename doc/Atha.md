# Self-alignment of Arabidopsis thaliana

## Detailed steps

### Prepare sequences

```bash
mkdir -p ~/data/alignment/example/Atha_self
cp -R ~/data/alignment/Ensembl/Atha ~/data/self_alignment
```

### self alignment

```bash
cd ~/data/alignment/example/Atha_self

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

On Mac i7-6700k

```bash
cd ~/data/alignment/example/Atha_self

if [ ! -d Atha_proc ]
then
    mkdir -p Atha_proc
fi

if [ ! -d Atha_result ]
then
    mkdir -p Atha_result
fi

cd ~/data/alignment/example/Atha_self/Atha_proc

# genome
find ../Atha -type f -name "*.fa" \
    | sort | xargs cat \
    | perl -nl -e '/^>/ or $_ = uc; print' \
    > genome.fa
faops size genome.fa > chr.sizes

# Get exact copies in the genome
fasops axt2fas ../Athavsselfalign/axtNet/*.axt.gz -l 1000 -s chr.sizes -o stdout > axt.fas
fasops separate axt.fas --nodash -s .sep.fasta

time perl ~/Scripts/egaz/sparsemem_exact.pl -f target.sep.fasta -g genome.fa \
    --length 500 -o replace.target.tsv
# real    0m56.162s
# user    1m29.831s
# sys     0m5.504s

fasops replace axt.fas replace.target.tsv -o axt.target.fas

perl ~/Scripts/egaz/sparsemem_exact.pl -f query.sep.fasta -g genome.fa \
    --length 500 -o replace.query.tsv
fasops replace axt.target.fas replace.query.tsv -o axt.correct.fas

# coverage stats
fasops covers axt.correct.fas -o axt.correct.yml
runlist split axt.correct.yml -s .temp.yml
runlist compare --op union target.temp.yml query.temp.yml -o axt.union.yml
runlist stat --size chr.sizes axt.union.yml -o union.csv

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
time perl ~/Scripts/egaz/fasta_blastn.pl  -f axt.gl.fasta -g genome.fa -o axt.bg.blast
# real    0m30.664s
# user    0m57.902s
# sys     0m3.393s
time perl ~/Scripts/egaz/blastn_genome.pl -f axt.bg.blast -g genome.fa -o axt.bg.fasta -c 0.95
# real    0m55.320s
# user    0m33.419s
# sys     0m24.626s

cat axt.gl.fasta axt.bg.fasta \
    | faops filter -u stdin stdout \
    | faops filter -n 250 stdin stdout \
    > axt.all.fasta

# link paralogs
time perl ~/Scripts/egaz/fasta_blastn.pl   -f axt.all.fasta -g axt.all.fasta -o axt.all.blast
# real    0m19.575s
# user    0m55.607s
# sys     0m1.503s
perl ~/Scripts/egaz/blastn_paralog.pl -f axt.all.blast -c 0.95 -o links.blast.tsv
```

### merge

```bash
cd ~/data/alignment/example/Atha_self/Atha_proc

# merge
rangeops sort -o links.sort.tsv \
    links.lastz.tsv links.blast.tsv

time rangeops clean   links.sort.tsv         -o links.sort.clean.tsv
#real	3m27.115s
#user	3m23.631s
#sys	0m1.952s
time rangeops merge   links.sort.clean.tsv   -o links.merge.tsv       -c 0.95 -p 8
#real	2m3.180s
#user	6m22.328s
#sys	0m3.126s
time rangeops clean   links.sort.clean.tsv   -o links.clean.tsv       -r links.merge.tsv --bundle 500
#real	1m34.452s
#user	1m30.996s
#sys	0m0.885s
rangeops connect links.clean.tsv        -o links.connect.tsv     -r 0.9
rangeops filter  links.connect.tsv      -o links.filter.tsv      -r 0.8

# recreate links
time rangeops create links.filter.tsv    -o multi.temp.fas       -g genome.fa
#real	0m24.954s
#user	0m8.294s
#sys	0m18.919s
time fasops   refine multi.temp.fas      -o multi.refine.fas     --msa mafft -p 8 --chop 10
#real	4m17.647s
#user	16m39.523s
#sys	9m30.666s
fasops   links  multi.refine.fas    -o stdout \
    | rangeops sort stdin -o links.refine.tsv

fasops   links  multi.refine.fas    -o stdout     --best \
    | rangeops sort stdin -o links.best.tsv
time rangeops create links.best.tsv      -o pair.temp.fas    -g genome.fa
#real	0m25.920s
#user	0m8.568s
#sys	0m19.423s
time fasops   refine pair.temp.fas       -o pair.refine.fas  --msa mafft -p 8
#real	4m1.471s
#user	15m51.691s
#sys	9m58.307s

cat links.refine.tsv \
    | perl -nla -F"\t" -e 'print for @F' \
    | runlist cover stdin -o cover.yml

echo "key,count" > links.count.csv
for n in 2 3 4 5-50
do
    rangeops filter links.refine.tsv -n ${n} -o stdout \
        > links.copy${n}.tsv

    cat links.copy${n}.tsv \
        | perl -nla -F"\t" -e 'print for @F' \
        | runlist cover stdin -o copy${n}.yml

    wc -l links.copy${n}.tsv \
        | perl -nl -e '
            @fields = grep {/\S+/} split /\s+/;
            next unless @fields == 2;
            next unless $fields[1] =~ /links\.([\w-]+)\.tsv/;
            printf qq{%s,%s\n}, $1, $fields[0];
        ' \
        >> links.count.csv

    rm links.copy${n}.tsv
done

runlist merge copy2.yml copy3.yml copy4.yml copy5-50.yml -o copy.all.yml
runlist stat --size chr.sizes copy.all.yml --mk --all -o links.copy.csv

cat links.copy.csv links.count.csv \
    | perl ~/Scripts/withncbi/util/merge_csv.pl --concat -o links.csv

runlist stat --size chr.sizes cover.yml
perl ~/Scripts/egaz/cover_figure.pl --size chr.sizes -f cover.yml
```

### result & clean

```bash
cd ~/data/alignment/example/Atha_self/Atha_proc

cp cover.yml        ../Atha_result/Atha.cover.yml
cp links.refine.tsv ../Atha_result/Atha.links.tsv
mv links.csv        ../Atha_result/Atha.links.csv
mv cover.yml.csv    ../Atha_result/Atha.cover.csv
mv cover.png        ../Atha_result/Atha.cover.png
mv multi.refine.fas ../Atha_result/Atha.multi.fas
mv pair.refine.fas  ../Atha_result/Atha.pair.fas

# clean
find . -type f -name "*genome.fa*" | xargs rm
find . -type f -name "*all.fasta*" | xargs rm
find . -type f -name "*.sep.fasta" | xargs rm
find . -type f -name "axt.*" | xargs rm
find . -type f -name "replace.*.tsv" | xargs rm
find . -type f -name "*.temp.yml" | xargs rm
find . -type f -name "*.temp.fas" | xargs rm
find . -type f -name "copy*.yml" | xargs rm
```
