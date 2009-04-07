# -*- makefile -*-

# Function to create a static library - as a serialised operation
# using a locking mechanism
CREATE_LIB=trap "rm -f $(LOCKFILE)" 2 3 9; \
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
	    touch $(LOCKFILE); \
	    echo $(AR) $@ $^; \
	    $(ARCH) $@ $^ && \
	    echo $(RANLIB) $@; \
	    $(RANLIB) $@ && \
	    rm -f $(LOCKFILE); \
	fi; \
	trap - 2 3 9;
