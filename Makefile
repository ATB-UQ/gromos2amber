
TEST_CASES := 1177_H2O 1177_H2O_melt
PRMTOPS := $(foreach X,$(TEST_CASES),temp/$X.parm7)
MDINFOS := $(foreach X,$(TEST_CASES),temp/$X.mdinfo)

.PHONY : test
test : $(PRMTOPS) $(MDINFOS)

.PHONY : temp
temp :
	mkdir -p temp

temp/% : test/% | temp
	cp $< $@ 

temp/%.parm7 : gromos2amber temp/%.top temp/%.g96
	./$< temp/$*.top temp/$*.g96 \
		temp/$*.parm7 temp/$*.rst7

temp/%.mdinfo: temp/measure_energy.in temp/%.parm7 #temp/%.rst7
	cd temp && \
	sander -O -i measure_energy.in \
	    -o $*_amber.out -p $*.parm7 -c $*.rst7 \
	    -r $*_amber.crd -inf $*_amber.mdinfo

temp/%.blah : temp/%.top temp/%.g96
	cd temp && \
	    md \
	    \@input measure_energy.run \
	    \@topo $*.top \
	    \@conf $*.g96 \
	    \@out  $*.out \
	    \@fin  $*.fin
	    \@trc  $*.trc
	


.PHONY : clean
clean :
	rm -r temp
