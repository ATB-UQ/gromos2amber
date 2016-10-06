
TEST_CASES := \
    1177_H2O_vacuum 1177_H2O_liquid \
    3798_H2O_vacuum 3798_H2O_liquid \
    32084_aa_H2O_vacuum \
    34231_aa_H2O_vacuum \
    32785_aa_H2O_vacuum \
    31084_aa_H2O_vacuum \
    24070_aa_H2O_vacuum \
    4482_aa_H2O_vacuum \
    34239_aa_H2O_vacuum \
    32133_aa_H2O_vacuum \
    33786_aa_H2O_vacuum \
    820_aa_H2O_vacuum \
    34011_aa_H2O_vacuum \
    32084_ua_H2O_vacuum \
    34231_ua_H2O_vacuum \
    32785_ua_H2O_vacuum \
    31084_ua_H2O_vacuum \
    24070_ua_H2O_vacuum \
    4482_ua_H2O_vacuum \
    34239_ua_H2O_vacuum \
    32133_ua_H2O_vacuum \
    33786_ua_H2O_vacuum \
    820_ua_H2O_vacuum \
    34011_ua_H2O_vacuum \
    DLPC_H2O_512_bilayer \
    DLPC_noH2O_512_bilayer
 
PRMTOPS := $(foreach X,$(TEST_CASES),temp/$X.amber.prmtop)
AENERGY := $(foreach X,$(TEST_CASES),temp/$X.amber.energy)
TRES := $(foreach X,$(TEST_CASES),temp/$X.gromos.tre)
IMDS := $(foreach X,$(TEST_CASES),temp/$X.gromos.imd)
YAML := $(foreach X,$(TEST_CASES),temp/$X.gromos.yml) \
    $(foreach X,$(TEST_CASES),temp/$X.amber.yml)
COMPARE := $(foreach X,$(TEST_CASES),$X_compare)
GROMOS2AMBER := ../gromos2amber

.PHONY : test
test : $(PRMTOPS) $(AENERGY) $(TRES) $(YAML) $(COMPARE) $(IMDS)

.PHONY : temp
temp :
	mkdir -p temp

temp/% : % | temp
	cp $< $@ 

temp/%_vacuum.amber.in : vacuum.amber.in
	cp $< $@

temp/%_liquid.amber.in : liquid.amber.in
	cp $< $@

temp/%_bilayer.amber.in : bilayer.amber.in
	cp $< $@

temp/%_vacuum.gromos.imd : vacuum.gromos.imd %_vacuum.gromos.cnf
	./make_imd $(word 2,$^) 3 < $< > $@

temp/%_liquid.gromos.imd : liquid.gromos.imd %_liquid.gromos.cnf
	./make_imd $(word 2,$^) 3 < $< > $@

temp/%_bilayer.gromos.imd : bilayer.gromos.imd %_bilayer.gromos.cnf
	./make_imd $(word 2,$^) 3 < $< > $@

temp/%.amber.prmtop : $(GROMOS2AMBER) temp/%.gromos.top temp/%.gromos.cnf
	./$< temp/$*.gromos.top temp/$*.gromos.cnf \
		temp/$*.amber.prmtop temp/$*.amber.inpcrd

temp/%.energy: temp/%.in temp/%.prmtop #temp/%.inpcrd
	cd temp && \
	sander -O -i $*.in \
	    -o $*.out -p $*.prmtop -c $*.inpcrd \
	    -r $*.crd -inf $*.mdinfo -e $*.energy


temp/%.tre : temp/%.imd temp/%.top temp/%.cnf
	md \
	    \@input temp/$*.imd \
	    \@topo temp/$*.top \
	    \@conf temp/$*.cnf \
	    \@fin  temp/$*.final.cnf \
	    \@trc  temp/$*.trc \
	    \@tre  temp/$*.tre > temp/$*.log

gromos_format.py : ../lib/gromos_format.py
	cp $< $@

temp/%.gromos.yml : temp/%.gromos.tre gromos_format.py
	./gromos_energy < $< > $@

temp/%.amber.yml : temp/%.amber.energy
	./amber_energy < $< > $@

#.PHONY : $(COMPARE)
%_compare : temp/%.amber.yml temp/%.gromos.yml
	echo $* && \
	    ./compare_energies $< $(word 2,$^)


.PHONY : clean
clean :
	rm -r temp gromos_format.py
