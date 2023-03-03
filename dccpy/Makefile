#
#        DCCPY module makefile
#
#----------------------------------------------------------------------------
# Project specific path defintions.
#----------------------------------------------------------------------------

M_BIN_DIR   = ../bin
M_LIB_DIR   = ../lib

PROJ_DIR    = .

SRC_DIR     = $(PROJ_DIR)

#----------------------------------------------------------------------------
# Include the appropriate compiler/platform definitions ...
#----------------------------------------------------------------------------
include ../etc/Makefile.platform

PYTHON_FILES   =   cifparse.py \
		   cns_inp.py \
		   config.py \
		   ligand.py \
		   matt.py \
		   ncs.py \
		   parse.py \
		   prog.py \
		   tls.py \
		   util.py \
		   edstat.pl \
		   write_cif.py \
                   sdsc_ligmap.py

PYTHON_PROG = dcc_calc.py

all: install

install: $(M_BIN_DIR)/dcc

$(M_BIN_DIR)/dcc:
	cd $(M_BIN_DIR) ; ln -s ../dccpy/dcc_calc.py dcc

clean: clean_build

clean_build:
	rm -f *.pyc
	cd $(M_BIN_DIR) ; rm -f dcc


export:
	mkdir -p $(EXPORT_DIR)
	@cp Makefile $(EXPORT_DIR)
	@$(INSTALL) -m 644 $(PYTHON_FILES) $(EXPORT_DIR)
	@$(INSTALL) -m 755 $(PYTHON_PROG) $(EXPORT_DIR)
