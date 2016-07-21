###################################################
#                                                 #
#   LICHEM: Layered Interacting CHEmical Models   #
#                                                 #
#        Symbiotic Computational Chemistry        #
#                                                 #
###################################################

### Standard compiler settings ###

CXX=g++
CXXFLAGS=-static -O3 -fopenmp

### Libarary settings ###

#The local copy of Eigen is located in ./Eigen3/
LDFLAGS=-I./Eigen3/

### Python settings ###

PYPATH=/usr/bin/python

### Sed commands ###

#In-place flag (GNU: -i, OSX: -i "")
SEDI=-i

### LaTeX settings ###

#For OSX these can be replaced with dummy calls to cat
TEX=pdflatex
BIB=bibtex

### Advanced compiler settings for developers ###

DEVFLAGS=-g -Wall -std=c++14
GPUFLAGS=-fopenacc

#####################################################

### Compile rules for users and devs ###

install:	title binary testexe manual compdone

Dev:	title devbin devtest manual stats compdone

GPUDev:	title gpubin devtest manual stats compdone

clean:	title delbin compdone

#####################################################

### Combine settings variables ###

# NB: Do not modify this section

FLAGSBIN=$(CXXFLAGS) $(LDFLAGS) -I./src/ -I./include/
FLAGSDEV=$(CXXFLAGS) $(DEVFLAGS) $(LDFLAGS) -I./src/ -I./include/
FLAGSGPU=$(CXXFLAGS) $(DEVFLAGS) $(GPUFLAGS) $(LDFLAGS) -I./src/ -I./include/

#####################################################

### Rules for building various parts of the code ###

binary:	
	@echo ""; \
	echo "### Compiling the LICHEM binary ###"; \
	mkdir -p bin
	$(CXX) ./src/LICHEM.cpp -o ./bin/lichem $(FLAGSBIN)
	@strip ./bin/lichem

devbin:	
	@echo ""; \
	echo "### Compiling the LICHEM development binary ###"; \
	mkdir -p bin
	$(CXX) ./src/LICHEM.cpp -o ./bin/lichem $(FLAGSDEV)

gpubin:	
	@echo ""; \
	echo "### Compiling the LICHEM GPU binary ###"; \
	mkdir -p bin
	$(CXX) ./src/LICHEM.cpp -o ./bin/lichem $(FLAGSGPU)

testexe:	
	@echo ""; \
	echo "### Creating test suite executable ###"
	@echo 'echo "#!$(PYPATH)" > ./tests/runtests'; \
	echo "!!$(PYPATH)" > ./tests/runtests
	cat ./src/runtests.py >> ./tests/runtests
	@sed $(SEDI) 's/\#.*//g' ./tests/runtests; \
	sed $(SEDI) 's/\s*$$//g' ./tests/runtests; \
	sed $(SEDI) '/^$$/d' ./tests/runtests; \
	sed $(SEDI) 's/\!\!/\#\!/g' ./tests/runtests; \
	chmod a+x ./tests/runtests

devtest:	
	@echo ""; \
	echo "### Creating development test suite executable ###"
	@echo 'echo "#!$(PYPATH)" > ./tests/runtests'; \
	echo "!!$(PYPATH)" > ./tests/runtests
	cat ./src/runtests.py >> ./tests/runtests
	@sed $(SEDI) 's/\#.*//g' ./tests/runtests; \
	sed $(SEDI) 's/\s*$$//g' ./tests/runtests; \
	sed $(SEDI) '/^$$/d' ./tests/runtests; \
	sed $(SEDI) 's/\!\!/\#\!/g' ./tests/runtests; \
	sed $(SEDI) 's/updateResults = 0/updateResults = 1/g' ./tests/runtests; \
	sed $(SEDI) 's/forceAll = 0/forceAll = 1/g' ./tests/runtests; \
	chmod a+x ./tests/runtests

checksyntax:	title
	@echo ""; \
	echo "### Checking for warnings and syntax errors ###"
	$(CXX) -fsyntax-only ./src/LICHEM.cpp -o lichem $(FLAGSDEV)
	@echo ""; \
	echo "### Source code statistics ###"; \
	echo "Number of LICHEM source code files:"; \
	ls -al include/* src/* | wc -l; \
	echo "Total length of LICHEM (lines):"; \
	cat include/* src/* | wc -l; \

manual:	
	@echo ""; \
	echo "### Compiling the documentation ###"; \
	cd src/; \
	echo "$(TEX) manual"; \
	$(TEX) manual > doclog.txt; \
	$(BIB) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	echo "$(BIB) manual"; \
	$(BIB) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	mv manual.pdf ../doc/LICHEM_manual.pdf; \
	rm -f manual.aux manual.bbl manual.blg; \
	rm -f manual.log manual.out manual.toc; \
	rm -f doclog.txt

title:	
	@echo ""; \
	echo "###################################################"; \
	echo "#                                                 #"; \
	echo "#   LICHEM: Layered Interacting CHEmical Models   #"; \
	echo "#                                                 #"; \
	echo "#        Symbiotic Computational Chemistry        #"; \
	echo "#                                                 #"; \
	echo "###################################################"

stats:	
	@echo ""; \
        echo "### Source code statistics ###"; \
	echo "Number of LICHEM source code files:"; \
	ls -al include/* src/* | wc -l; \
	echo "Total length of LICHEM (lines):"; \
	cat include/* src/* | wc -l

compdone:	
	@echo ""; \
	echo "Done."; \
	echo ""

delbin:	
	@echo ""; \
	if grep -q "JOKES = 1" include/LICHEM_options.h; then \
	echo '     ___'; \
	echo '    |_  |'; \
	echo '      \ \'; \
	echo '      |\ \'; \
	echo '      | \ \'; \
	echo '      \  \ \'; \
	echo '       \  \ \'; \
	echo '        \  \ \       <wrrr vroooom wrrr> '; \
	echo '         \__\ \________'; \
	echo '             |_________\'; \
	echo '             |__________|  ..,  ,.,. .,.,, ,..'; \
	echo ""; \
 	fi; \
        echo ""; \
	echo "Removing binary and manual..."; \
	rm -rf lichem ./doc/LICHEM_manual.pdf ./tests/runtests ./bin

