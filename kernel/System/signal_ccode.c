#include <signal.h>
#include <unistd.h>

#define NSIGS 32

typedef void (*sighandler_t)(int);
static struct sigaction actions[NSIGS];

/**
 * Helper routine. This routine is registered as signal handler
 * for all registered handlers. If it is raised, it determines
 * the actual signal handler and calls it, whereby the signal
 * number is provided "by-reference".
 */
void signal_raised(int sig) {
    (void) actions[sig].sa_handler((int) &sig);
}

/**
 * Register a signal handler to a given signal. If another signal handler is
 * already registered to the signal, then the old one is deregistered and
 * the new one is used instead.
 */
static sighandler_t
register_handler(int sig, sighandler_t func, sighandler_t catch) {
    struct sigaction new_sig, old_sig;

    new_sig.sa_handler = catch;
    sigemptyset (&new_sig.sa_mask);
    new_sig.sa_flags = SA_RESTART;
    if (sigaction (sig, &new_sig, &old_sig) < 0)
    return SIG_ERR;

    if (sig > NSIGS)
    return SIG_ERR;
    actions[sig].sa_handler = func;
    return old_sig.sa_handler;
}

/**
 * Deregister a signal handler for a given signal
 */
static sighandler_t
deregister_handler (int sig) {
    struct sigaction del_sig;

    if (sigaction (sig, NULL, &del_sig) < 0)
    return SIG_ERR;
    del_sig.sa_handler = SIG_DFL;
    sigdelset (&del_sig.sa_mask, sig);
    del_sig.sa_flags = SA_RESTART;
    if (sigaction (sig, &del_sig, NULL) < 0)
    return SIG_ERR;

    actions[sig].sa_handler = SIG_DFL;
    return del_sig.sa_handler;
}

/**
 * Registers a signal handler with a given signal
 */
void signal_register(int *sig, void func(int)) {
    register_handler(*sig, func, signal_raised);
}
    
/**
 * Registers a signal handler with a given signal
 */
void signal_register_(int *sig, void func(int)) {
    register_handler(*sig, func, signal_raised);
}

/**
 * Registers a signal handler with a given signal
 */
void signal_register__(int *sig, void func(int)) {
    register_handler(*sig, func, signal_raised);
}

/**
 * Deregisters a signal handler for a given signal
 */
void signal_deregister(int *sig) {
    deregister_handler(*sig);
}

/**
 * Deregisters a signal handler for a given signal
 */
void signal_deregister_(int *sig) {
    deregister_handler(*sig);
}

/**
 * Deregisters a signal handler for a given signal
 */
void signal_deregister__(int *sig) {
    deregister_handler(*sig);
}
