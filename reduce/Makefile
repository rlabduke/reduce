#
# reduce top level Makefile
#

SUBDIRS1 = toolclasses libpdb reduce_src
SUBDIRS = toolclasses libpdb reduce_src nosym_reduce_src

all:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	echo "Start Compling:"
	for i in $(SUBDIRS1) ;\
	do \
		(cd $$i ; echo "making" all "in $(CURRENT_DIR)/$$i..."; \
			$(MAKE) $(MFLAGS) 'CXXDEBUGFLAGS=$(CXXDEBUGFLAGS)' all); \
	done

nosym:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	echo "Start Compling:"
	for i in $(SUBDIRS1) ;\
	do \
		(cd $$i ; echo "making" all "in $(CURRENT_DIR)/$$i..."; \
			$(MAKE) $(MFLAGS) USE_SYM=false 'CXXDEBUGFLAGS=$(CXXDEBUGFLAGS)' all); \
	done

install:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) ;\
	do \
		(cd $$i ; echo "installing" "$(CURRENT_DIR)/$$i..."; \
			$(MAKE) $(MFLAGS) 'CXXDEBUGFLAGS=$(CXXDEBUGFLAGS)' install); \
	done

clean:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) ;\
	do \
		(cd $$i ; echo "cleaning" "in $(CURRENT_DIR)/$$i..."; \
			$(MAKE) $(MFLAGS)  clean); \
	done
