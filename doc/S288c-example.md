# Self-alignment of S288c

## Detailed steps

### Prepare sequences

```bash
mkdir -p ~/data/alignment/example/S288c_self
cp -R ~/data/alignment/example/scer/Genomes/S288c ~/data/alignment/example/S288c_self
```

### self alignment

```bash
cd ~/data/alignment/example/S288c_self

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
cd ~/data/alignment/example/S288c_self

if [ ! -d S288c_proc ]
then
    mkdir S288c_proc
fi

if [ ! -d S288c_result ]
then
    mkdir S288c_result
fi

cd ~/data/alignment/example/S288c_self/S288c_proc

# genome
find ../S288c -type f -name "*.fa" \
    | sort | xargs cat \
    | perl -nl -e '/^>/ or $_ = uc; print' \
    > genome.fa
faops size genome.fa > chr.sizes

# Get exact copies in the genome
fasops axt2fas ../S288cvsselfalign/axtNet/*.axt.gz -l 1000 -s chr.sizes -o stdout > axt.fas
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
perl ~/Scripts/egas/fasta_blastn.pl  -f axt.gl.fasta -g genome.fa -o axt.bg.blast
perl ~/Scripts/egas/blastn_genome.pl -f axt.bg.blast -g genome.fa -o axt.bg.fasta -c 0.95

cat axt.gl.fasta axt.bg.fasta \
    | faops filter -u stdin stdout \
    | faops filter -n 250 stdin stdout \
    > axt.all.fasta

# link paralogs
perl ~/Scripts/egas/fasta_blastn.pl   -f axt.all.fasta -g axt.all.fasta -o axt.all.blast
perl ~/Scripts/egas/blastn_paralog.pl -f axt.all.blast -c 0.95 -o links.blast.tsv

# merge
rangeops merge   links.lastz.tsv    links.blast.tsv -o links.merge.tsv -c 0.95 -p 8
rangeops sort    links.lastz.tsv    links.blast.tsv -o links.sort.tsv
rangeops clean   links.sort.tsv  -r links.merge.tsv -o links.clean.tsv
rangeops connect links.clean.tsv                    -o links.connect.tsv
rangeops filter  links.connect.tsv                  -o links.filter.tsv -r 0.8

cat links.filter.tsv \
    | perl -nla -F"\t" -e 'print for @F' \
    | runlist cover stdin -o cover.yml

# fasops create links.filter.tsv
# fasops refine
# fasops links --best
# fasops create links.pairwise.tsv

echo "key,count" > links.count.csv
for n in 2 3 4 5-50
do
    rangeops filter links.filter.tsv -n ${n} -o stdout \
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
perl ~/Scripts/egas/cover_figure.pl --size chr.sizes -f cover.yml
```

### result & clean

```bash
cd ~/data/alignment/example/S288c_selfS288c_proc

cp cover.yml        ../S288c_result/S288c.cover.yml
cp links.filter.tsv ../S288c_result/S288c.links.tsv
mv links.csv        ../S288c_result/S288c.links.csv
mv cover.yml.csv    ../S288c_result/S288c.cover.csv
mv cover.png        ../S288c_result/S288c.cover.png

# clean
find . -type f -name "*genome.fa*" | xargs rm
find . -type f -name "*all.fasta*" | xargs rm
find . -type f -name "*.sep.fasta" | xargs rm
find . -type f -name "axt.*" | xargs rm
find . -type f -name "replace.*.tsv" | xargs rm
find . -type f -name "*.temp.yml" | xargs rm
find . -type f -name "copy*.yml" | xargs rm
```

## Use `self_batch.pl`

```bash
mkdir -p ~/data/alignment/example/self_batch
cp -R ~/data/alignment/example/scer/Genomes/S288c ~/data/alignment/example/self_batch

cd ~/data/alignment/example/self_batch

perl ~/Scripts/withncbi/taxon/strain_info.pl \
    --file   yeast_taxon.csv \
    --id    559292 \
    --name 559292=S288c

perl ~/Scripts/egaz/self_batch.pl \
    -c ~/data/alignment/example/self_batch/yeast_taxon.csv \
    --working_dir ~/data/alignment/example/self_batch \
    --seq_dir ~/data/alignment/example/self_batch \
    --length 1000  \
    --norm \
    --name yeast \
    --parallel 8 \
    -t S288c

cd ~/data/alignment/example/self_batch
bash yeast/1_real_chr.sh
bash yeast/3_self_cmd.sh
bash yeast/4_proc_cmd.sh
# real	1m8.817s
# user	2m39.879s
# sys	1m18.005s
bash yeast/5_circos_cmd.sh
bash yeast/6_feature_cmd.sh

```
