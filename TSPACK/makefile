# Makefile for TSPACK

TSSRC = ENDSLP.f SIGS.f SNHCSH.f STORE.f YPCOEF.f YPC1.f YPC1P.f YPC2.f YPC2P.f TSPSI.f \
 INTRVL.f HVAL.f HPVAL.f HPPVAL.f TSINTL.f TSVAL1.f

all : TSPACK.o

TSPACK.o : $(TSSRC) makefile
	cat $(TSSRC) > TSPACK.f
	$(FF) $(FOPT) -c -o TSPACK.o TSPACK.f
	ar rv libTSPACK.a *.o
	ranlib libTSPACK.a
	rm -f TSPACK.f *.o

clean:
	rm -f TSPACK.f *.o *.a
