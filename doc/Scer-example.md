# Several approaches to build genomic alignments of two *S. cerevisiae* strains, S288c and RM11_1a

## Detailed steps

### Prepare sequences

```bash
mkdir -p ~/data/alignment/egaz
cd ~/data/alignment/egaz

cp -R ~/data/alignment/example/scer/Genomes/S288c   .
cp -R ~/data/alignment/example/scer/Genomes/RM11_1a .

cat ./S288c/*.fa   > S288c.fa
cat ./RM11_1a/*.fa > RM11_1a.fa

faops size S288c.fa > S288c.chr.sizes
faops size RM11_1a.fa > RM11_1a.chr.sizes
```

### Vanilla lastz

First sequence in target file will be aligned against all sequences in query file.

```bash
cd ~/data/alignment/egaz
lastz S288c/I.fa RM11_1a.fa > I.lav
```

### With `bz.pl` 

Use lastz's chaining.

`bz.pl` has build-in lav2axt.

```bash
cd ~/data/alignment/egaz

perl ~/Scripts/egaz/bz.pl \
    -p 8 \
    -dt S288c \
    -dq RM11_1a \
    -dl S288cvsRM11_1a_bz
```

Check axt positions.

```bash
fasops axt2fas S288cvsRM11_1a_bz/*.axt -l 1000 -t S288c -q RM11_1a -s RM11_1a.chr.sizes \
    -o axt.bz.fas

fasops check axt.bz.fas S288c.fa -n S288c -o stdout | grep -v "OK"
fasops check axt.bz.fas RM11_1a.fa -n RM11_1a -o stdout | grep -v "OK"
```

### With `bz.pl` and `lpcna.pl` 

Use UCSC chaining.

```bash
cd ~/data/alignment/egaz

perl ~/Scripts/egaz/bz.pl \
    -p 8 -C 0 --noaxt \
    -dt S288c \
    -dq RM11_1a \
    -dl S288cvsRM11_1a_ucsc

perl ~/Scripts/egaz/lpcna.pl \
    -p 8 \
    -dt S288c \
    -dq RM11_1a \
    -dl S288cvsRM11_1a_ucsc
```

Check axt positions.

```bash
fasops axt2fas S288cvsRM11_1a_ucsc/axtNet/*.axt.gz -l 1000 -t S288c -q RM11_1a -s RM11_1a.chr.sizes \
    -o axt.ucsc.fas

fasops check axt.ucsc.fas S288c.fa -n S288c -o stdout | grep -v "OK"
fasops check axt.ucsc.fas RM11_1a.fa -n RM11_1a -o stdout | grep -v "OK"
```

### Self alignment

```bash
cd ~/data/alignment/egaz

perl ~/Scripts/egaz/bz.pl \
    --is_self \
    -p 8 -C 0 --noaxt \
    -dt S288c \
    -dq S288c \
    -dl S288cvsselfalign

perl ~/Scripts/egaz/lpcna.pl \
    -p 8 \
    -dt S288c \
    -dq S288c \
    -dl S288cvsselfalign
```

Check axt positions.

```bash
fasops axt2fas S288cvsselfalign/axtNet/*.axt.gz -l 1000 -t S288c -q S288c -s S288c.chr.sizes \
    -o axt.self.fas

fasops check axt.self.fas S288c.fa -n S288c -o stdout | grep -v "OK"
```

### Partitions

```bash
cd ~/data/alignment/egaz

rm -fr S288c_parted RM11_1a_parted
perl ~/Scripts/egaz/part_seq.pl -i S288c -o S288c_parted --chunk 500000 --overlap 10000
perl ~/Scripts/egaz/part_seq.pl -i RM11_1a  -o RM11_1a_parted  --chunk 500000 --overlap 0

rm -fr S288cvsRM11_1a_parted
perl ~/Scripts/egaz/bz.pl \
    -p 8 \
    -dt S288c_parted -tp \
    -dq RM11_1a_parted -qp \
    -dl S288cvsRM11_1a_parted
```

Check axt positions.

```bash
fasops axt2fas S288cvsRM11_1a_parted/*.axt -l 1000 -t S288c -q RM11_1a -s RM11_1a.chr.sizes \
    -o axt.parted.fas

fasops check axt.parted.fas S288c.fa -n S288c -o stdout | grep -v "OK"
fasops check axt.parted.fas RM11_1a.fa -n RM11_1a -o stdout | grep -v "OK"
```
