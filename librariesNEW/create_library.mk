# -*- makefile -*-

# Function to create a static library - as a serialised operation
# using a locking mechanism and a random wait time at the beginning
# (as we don't have another approach yet to catch two processes that
#  execute these instructions exactly concurrently - which both set
#  the lock file simultaneously.)
CREATE_LIB=sleep `echo $${PPID} | sed -e 's/^.*\(.\)$$/\1/'`; \
	if test -f $(LOCKFILE); then \
	    loop=0; \
	    while test -f $(LOCKFILE) -a $${loop} -lt $(RETRIES); do \
		echo; \
		echo "\#"; \
		echo "\# Detected a concurrent attempt to create"; \
		echo "\# <$@>."; \
		echo "\# Waiting for the other process to finish."; \
		echo "\#"; \
		echo "\# If you are sure that this is a false alarm, simply remove the file"; \
		echo "\# <"`pwd`"/$(LOCKFILE)>."; \
		echo "\# Then, compilation will automatically proceed."; \
		echo "\#"; \
		echo; \
		sleep $(WAITTIME); \
		loop=`expr $${loop} + 1`; \
	    done; \
	fi; \
	errorcode=0; \
	if test -f $(LOCKFILE); then \
	    echo; \
	    echo "\#"; \
	    echo "\# Lock for creating $@"; \
	    echo "\# was not released within "`expr $(WAITTIME) \* $(RETRIES)`" seconds."; \
	    echo "\# Giving up."; \
	    echo "\#"; \
	    echo "\# In case of a false alarm please manually remove the file"; \
	    echo "\# <$(LOCKFILE)>."; \
	    echo "\#"; \
	    echo; \
	    exit 1; \
	else \
	    trap "rm -f $(LOCKFILE)" 2 3 9; \
	    touch $(LOCKFILE); \
	    echo $(ARCH) $@ $^; \
	    ( $(ARCH) $@ $^ && \
	      echo $(RANLIB) $@ && \
	      $(RANLIB) $@ ) || errorcode=1; \
	    rm -f $(LOCKFILE); \
	fi; \
	trap - 2 3 9; \
	exit $${errorcode};
