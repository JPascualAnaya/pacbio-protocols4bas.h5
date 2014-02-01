FILES=$(shell cat input.fofn)
FILENAMES=$(notdir $(FILES))
SMRT_CELL_IDS=$(subst .bax.h5,,$(FILENAMES))
FASTA_FILES=$(addsuffix .fasta,$(SMRT_CELL_IDS))
MASKED_FILES=$(addsuffix .masked.bax.h5,$(SMRT_CELL_IDS))
UPDATED_REGIONS=$(addsuffix .regions.log,$(SMRT_CELL_IDS))

all : $(UPDATED_REGIONS) quiver.fasta

.SECONDARY :

.PHONY : source clean

source :
	rsync -arv --bwlimit=10000 $(FILES) .

quiver.fasta : data/polished_assembly.fasta.gz
	gunzip -c $< | sed 's/|quiver//' > $@

data/polished_assembly.fasta.gz : input.masked.xml
	smrtpipe.py --distribute --params=$(MAKEDIR)/HGAP_Assembly_Advanced.1.xml xml:$<

input.masked.xml : input.masked.fofn
	fofnToSmrtpipeInput.py $< > $@

input.masked.fofn : $(UPDATED_REGIONS)
	find ./ -name "*.bax.h5" -exec readlink -f {} \; > $@

%.regions.log : %.non_vector.bed
	python $(MAKEDIR)/update_regions.py $* $< 2> $@

%.non_vector.bed : %.screen.bed %.read_lengths.screen.tab
	bedtools complement -i $*.screen.bed -g $*.read_lengths.screen.tab > $@
	cut -f 1 $< | sort | uniq > $*.tmp_all_holes.txt
	cut -f 1 $@ | sort | uniq > $*.tmp_holes_with_non_vector.txt
	comm -23 $*.tmp_all_holes.txt $*.tmp_holes_with_non_vector.txt \
		| awk 'OFS="\t" { print $$1,0,0 }' >> $@
	rm -f $*.tmp_*.txt

%.read_lengths.screen.tab : %.reads.txt %.read_lengths.tab
	join -j 1 $^ | sed 's/\s\+/\t/g' > $@

%.reads.txt : %.screen.bed
	cut -f 1 $^ | sort -k 1,1 | uniq > $@

%.screen.bed : vectors.txt %.screen.out
	grep -f $^ \
		| awk 'OFS="\t" { if ($$7-$$6 >= 100) { print $$5,$$6,$$7 } }' \
		| sort -k 1,1 -k 2,2n \
		| mergeBed -i stdin -d 500 > $@

vectors.txt :
	grep "^>" $(VECTOR) | sed 's/>//' > $@

%.screen.out : %.fasta
	cross_match $^ $(VECTOR) -minmatch 10 -minscore 20 -screen > $@

%.read_lengths.tab : %.fasta
	fastalength $^ | awk 'OFS="\t" { print $$2,$$1 }' | sort -k 1,1 > $@

%.fasta : %.masked.bax.h5
	pls2fasta $^ $@ -noSplitSubreads

%.masked.bax.h5 :
	mv $*.bax.h5 $@

clean :
	rm -f filtered_longreads.fasta index* metadata.rdf noop quiver.fasta toc.xml variants.*
	rm -rf data/ log/ movie_metadata/ reference/ results/ workflow/
