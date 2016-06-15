# Alignments of *S. cerevisiae* strains

## Detailed steps

### Prepare sequences

```bash
mkdir -p ~/Scripts/egaz/data
cp -R ~/data/alignment/example/scer/Genomes/S288c ~/Scripts/egaz/data/
cp -R ~/data/alignment/example/scer/Genomes/RM11_1a/ ~/Scripts/egaz/data/RM11

cat ~/Scripts/egaz/data/S288c/*.fa > ~/Scripts/egaz/data/S288c.fa
cat ~/Scripts/egaz/data/RM11/*.fa > ~/Scripts/egaz/data/RM11.fa

faops size S288c.fa > S288c.chr.sizes
faops size RM11.fa > RM11.chr.sizes
```

### Vanilla lastz

First sequence in target file will be aligned against all sequences in query file.

```bash
cd ~/Scripts/egaz/data
lastz S288c/I.fa RM11.fa > I.lav
```

### With `bz.pl` 

Use lastz's chaining.

`bz.pl` has build-in lav2axt.

```bash
cd ~/Scripts/egaz/data

perl ~/Scripts/egaz/bz.pl \
    -p 8 \
    -dt S288c \
    -dq RM11 \
    -dl S288cvsRM11_bz
```

Check axt positions.

```bash
fasops axt2fas S288cvsRM11_bz/*.axt -l 1000 -t S288c -q RM11 -s RM11.chr.sizes \
    -o stdout > axt.bz.fas

perl ~/Scripts/alignDB/util/check_header.pl \
    -i axt.bz.fas \
    -n S288c \
    -g S288c.fa \
    --detail
perl ~/Scripts/alignDB/util/check_header.pl \
    -i axt.bz.fas \
    -n RM11 \
    -g RM11.fa \
    --detail
```

### With `bz.pl` and `lpcna.pl` 

Use UCSC chaining.

```bash
cd ~/Scripts/egaz/data

perl ~/Scripts/egaz/bz.pl \
    -p 8 -C 0 --noaxt \
    -dt S288c \
    -dq RM11 \
    -dl S288cvsRM11_ucsc

perl ~/Scripts/egaz/lpcna.pl \
    -p 8 \
    -dt S288c \
    -dq RM11 \
    -dl S288cvsRM11_ucsc
```

Check axt positions.

```bash
fasops axt2fas S288cvsRM11_ucsc/axtNet/*.axt.gz -l 1000 -t S288c -q RM11 -s RM11.chr.sizes \
    -o stdout > axt.ucsc.fas

perl ~/Scripts/alignDB/util/check_header.pl \
    -i axt.ucsc.fas \
    -n S288c \
    -g S288c.fa \
    --detail
perl ~/Scripts/alignDB/util/check_header.pl \
    -i axt.ucsc.fas \
    -n RM11 \
    -g RM11.fa \
    --detail
```

### Self alignment

```bash
cd ~/Scripts/egaz/data

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
fasops axt2fas S288cvsselfalign/axtNet/*.axt.gz -l 1000 -t S288c -q RM11 -s S288c.chr.sizes \
    -o stdout > axt.self.fas

perl ~/Scripts/alignDB/util/check_header.pl \
    -i axt.self.fas \
    -n S288c \
    -g S288c.fa \
    --detail
```

### Partitions


```bash
cd ~/Scripts/egaz/data

rm -fr S288c_parted RM11_parted
perl ~/Scripts/egaz/part_seq.pl -i S288c -o S288c_parted --chunk 500000 --overlap 10000
perl ~/Scripts/egaz/part_seq.pl -i RM11  -o RM11_parted  --chunk 500000 --overlap 0

rm -fr S288cvsRM11_parted
perl ~/Scripts/egaz/bz.pl \
    -p 8 \
    -dt S288c_parted -tp \
    -dq RM11_parted -qp \
    -dl S288cvsRM11_parted

```

Check axt positions.

```bash
fasops axt2fas S288cvsRM11_parted/*.axt -l 1000 -t S288c -q RM11 -s RM11.chr.sizes \
    -o stdout > axt.parted.fas

perl ~/Scripts/alignDB/util/check_header.pl \
    -i axt.parted.fas \
    -n S288c \
    -g S288c.fa \
    --detail
perl ~/Scripts/alignDB/util/check_header.pl \
    -i axt.parted.fas \
    -n RM11 \
    -g RM11.fa \
    --detail
```
