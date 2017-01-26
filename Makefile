
all: 
	cd source; $(MAKE) all

tester:
	cd source; $(MAKE) tester

getplots:
	cd source; $(MAKE) getplots

install:
	mv dmbayes_chord tester bin

clean:
	cd source;  $(MAKE) cleanall
