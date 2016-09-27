
TEST_CASES := 1177_H2O
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
	    -o $*.out -p $*.parm7 -c $*.rst7 -r $*.crd -inf $*.mdinfo

.PHONY : clean
clean :
	rm -r temp
