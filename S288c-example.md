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

### coverage stat

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

fasops axt2fas ../S288cvsselfalign/axtNet/*.axt.gz -l 1000 -o stdout > S288cvsselfalign.axt.fas
fasops covers S288cvsselfalign.axt.fas -n target -o S288cvsselfalign.target.yml
runlist stat --size ../S288c/chr.sizes S288cvsselfalign.target.yml -o ../S288c_result/S288c.1000.csv

```

###  blast and merge

```bash
cd ~/Scripts/egas/data/S288c_proc

# genome
find ../S288c -type f -name "*.fa" \
    | sort | xargs cat \
    | perl -nl -e '/^>/ or uc; print' \
    > S288c.genome.fa

~/share/blast/bin/formatdb -p F -o T -i S288c.genome.fa

~/share/blast/bin/blastall -p blastn -F "m D" -m 0 -b 10 -v 10 -e 1e-3 -a 8 -i S288cvsselfalign.axt.fas -d S288c.genome.fa -o S288cvsselfalign.axt.blast

# Coreect queries positions in axt files
perl ~/Scripts/egas/blastn_genome_location.pl -f S288cvsselfalign.axt.blast -m 0 -i 90 -c 0.95

# paralog
~/share/blast/bin/formatdb -p F -o T -i S288cvsselfalign.gl.fasta

~/share/blast/bin/blastall -p blastn -F "m D" -m 0 -b 10 -v 10 -e 1e-3 -a 8 -i S288cvsselfalign.gl.fasta -d S288cvsselfalign.gl.fasta -o S288cvsselfalign.gl.blast

# merge
perl ~/Scripts/egas/blastn_paralog.pl -f S288cvsselfalign.gl.blast -m 0 -i 90 -c 0.9

perl ~/Scripts/egas/merge_node.pl    -v -f S288cvsselfalign.blast.tsv -o S288cvsselfalign.merge.yml -c 0.9
perl ~/Scripts/egas/paralog_graph.pl -v -f S288cvsselfalign.blast.tsv -m S288cvsselfalign.merge.yml --nonself -o S288cvsselfalign.merge.graph.yml
perl ~/Scripts/egas/cc.pl               -f S288cvsselfalign.merge.graph.yml
perl ~/Scripts/egas/proc_cc_chop.pl     -f S288cvsselfalign.cc.yml --size ../S288c/chr.sizes --msa mafft
perl ~/Scripts/egas/proc_cc_stat.pl     -f S288cvsselfalign.cc.yml --size ../S288c/chr.sizes

runlist stat --size ../S288c/chr.sizes S288cvsselfalign.cc.chr.runlist.yml;
```

### result & clean

```bash
cd ~/Scripts/egas/data/S288c_proc

cp S288cvsselfalign.cc.yml ../S288c_result
mv S288cvsselfalign.cc.csv ../S288c_result
cp S288cvsselfalign.cc.chr.runlist.yml.csv ../S288c_result/S288cvsselfalign.chr.csv

# clean
find . -type f -name "*genome.fa*" | xargs rm
find . -type f -name "*gl.fasta*" | xargs rm
find . -type f -name "*.blast" | xargs rm
```
