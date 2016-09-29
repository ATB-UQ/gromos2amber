
TEST_CASES := 1177_H2O 1177_H2O_melt
PRMTOPS := $(foreach X,$(TEST_CASES),temp/$X.parm7)
MDINFOS := $(foreach X,$(TEST_CASES),temp/$X.mdinfo)
GLOGS := $(foreach X,$(TEST_CASES),temp/$X.log)

.PHONY : test
test : $(PRMTOPS) $(MDINFOS) $(GLOGS)

.PHONY : temp
temp :
	mkdir -p temp

temp/% : test/% | temp
	cp $< $@ 

temp/%.parm7 : gromos2amber temp/%.top temp/%.g96
	./$< temp/$*.top temp/$*.g96 \
		temp/$*.parm7 temp/$*.rst7

temp/%.mdinfo: temp/%.in temp/%.parm7 #temp/%.rst7
	cd temp && \
	sander -O -i $*.in \
	    -o $*_amber.out -p $*.parm7 -c $*.rst7 \
	    -r $*_amber.crd -inf $*_amber.mdinfo


temp/%.log : temp/%.imd temp/%.top temp/%.g96
	md \
	    \@input temp/$*.imd \
	    \@topo temp/$*.top \
	    \@conf temp/$*.g96 \
	    \@fin  temp/$*.g96 \
	    \@trc  temp/$*.trc \
	    \@tre  temp/$*.tre > temp/$*.log
	


.PHONY : clean
clean :
	rm -r temp
