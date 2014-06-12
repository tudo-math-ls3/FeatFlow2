# -*- makefile -*-

# Function to unpack a tarball - as a serialised operation
# using a locking mechanism
#
# Logic:
# If a lock file exists, there is already a concurrent process unpacking the
# tarball. In this case, we merely need to inform the user why this process
# won't extract the tarball - to catch case where the lock file wrongfully
# lingers. But by design it is not necessary to afterwards start extracting
# the tarball, too.
# If no lock file exists, extract the tarball and patch the sources, if
# necessary.
UNPACK=if test -f $(LOCKFILE); then \
	    loop=0; \
	    while test -f $(LOCKFILE) -a $${loop} -lt $(RETRIES); do \
		echo; \
		echo "\#"; \
		echo "\# Detected a concurrent attempt to unpack $(TARBALL)."; \
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
		echo "\# Lock for unpacking $(TARBALL) has not been released"; \
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
	    echo "\# Extracting $(NAME)..."; \
	    gzip -dc $(TARBALL) | ( cd ..; tar -xvf - ); \
	    if test "echoX" != '$(PATCHCMD)X'; then \
		test -n '$(PATCHTEXT1)' && echo $(PATCHTEXT1); \
		test -n '$(PATCHTEXT2)' && echo $(PATCHTEXT2); \
		test -n '$(PATCHTEXT3)' && echo $(PATCHTEXT3); \
		echo "$(PATCHCMD)"; \
		$(PATCHCMD); \
		echo "\# Sources patched."; \
	    fi; \
	    rm -f $(LOCKFILE); \
	fi; \
	trap - 2 3 9;
