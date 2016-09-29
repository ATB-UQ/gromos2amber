
TEST_CASES := 1177_H2O_vacuum 1177_H2O_liquid
PRMTOPS := $(foreach X,$(TEST_CASES),temp/$X.amber.prmtop)
MDINFOS := $(foreach X,$(TEST_CASES),temp/$X.amber.mdinfo)
GLOGS := $(foreach X,$(TEST_CASES),temp/$X.gromos.log)
GROMOS2AMBER := ../gromos2amber

.PHONY : test
test : $(PRMTOPS) $(MDINFOS) $(GLOGS)

.PHONY : temp
temp :
	mkdir -p temp

temp/% : % | temp
	cp $< $@ 

temp/%.amber.prmtop : $(GROMOS2AMBER) temp/%.gromos.top temp/%.gromos.cnf
	./$< temp/$*.gromos.top temp/$*.gromos.cnf \
		temp/$*.amber.prmtop temp/$*.amber.inpcrd

temp/%.mdinfo: temp/%.in temp/%.prmtop #temp/%.inpcrd
	cd temp && \
	sander -O -i $*.in \
	    -o $*.out -p $*.prmtop -c $*.inpcrd \
	    -r $*.crd -inf $*.mdinfo


temp/%.log : temp/%.imd temp/%.top temp/%.cnf
	md \
	    \@input temp/$*.imd \
	    \@topo temp/$*.top \
	    \@conf temp/$*.cnf \
	    \@fin  temp/$*.final.cnf \
	    \@trc  temp/$*.trc \
	    \@tre  temp/$*.tre > temp/$*.log

.PHONY : clean
clean :
	rm -r temp
