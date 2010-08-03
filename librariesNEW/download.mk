# -*- makefile -*-

# Function to download a tarball - as a serialised operation
# using a locking mechanism
DOWNLOAD=if test -f $(LOCKFILE); then \
	    loop=0; \
	    while test -f $(LOCKFILE) -a $${loop} -lt $(RETRIES); do \
		echo; \
		echo "\#"; \
		echo "\# Detected a concurrent attempt to download $(TARBALL)."; \
		echo "\# Waiting for the other process to finish."; \
		echo "\#"; \
		echo "\# If you are sure that this is a false alarm, simply remove the file"; \
		echo "\# <"`pwd`"/$(LOCKFILE)>."; \
		echo "\# Then, compilation will proceed automatically."; \
		echo "\#"; \
		echo; \
		sleep $(WAITTIME); \
		loop=`expr $${loop} + 1`; \
	    done; \
	    if test -f $(LOCKFILE); then \
		echo; \
		echo "\#"; \
		echo "\# Lock for downloading $(TARBALL) has not been released"; \
		echo "\# within "`expr $(WAITTIME) \* $(RETRIES)`" seconds. Giving up."; \
		echo "\#"; \
		echo "\# In case of a false alarm please manually remove the file"; \
		echo "\# <$(LOCKFILE)>."; \
		echo "\#"; \
		echo; \
		exit 1; \
	    fi; \
	else \
	    trap "rm -f $(LOCKFILE)" 2 3 9; \
	    touch $(LOCKFILE); \
	    wget -O $(TARBALL) $(URL); \
	    rm -f $(LOCKFILE); \
	fi; \
	trap - 2 3 9;
