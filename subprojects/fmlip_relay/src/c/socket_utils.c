/*
 * socket_utils.c
 * --------------
 * Thin C wrappers around POSIX sockets for use from Fortran via ISO_C_BINDING.
 * Compile with:  gcc -O2 -c socket_utils.c -o socket_utils.o
 */

#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <stdio.h>

/* ── connect ──────────────────────────────────────────────────────────────── */

/**
 * sock_connect_c(host_c, port) -> file descriptor, or -1 on error
 * host_c must be a null-terminated ASCII string (e.g. "127.0.0.1")
 */
int sock_connect_c(const char *host_c, int port)
{
    int fd = socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) { perror("socket"); return -1; }

    /* disable Nagle for low-latency request/response */
    int one = 1;
    setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one));

    struct sockaddr_in addr;
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port   = htons((uint16_t)port);
    if (inet_pton(AF_INET, host_c, &addr.sin_addr) <= 0) {
        fprintf(stderr, "inet_pton failed for '%s'\n", host_c);
        close(fd);
        return -1;
    }

    if (connect(fd, (struct sockaddr *)&addr, sizeof(addr)) < 0) {
        perror("connect");
        close(fd);
        return -1;
    }
    return fd;
}

/* ── send / recv (blocking, exactly n bytes) ─────────────────────────────── */

/**
 * Returns 0 on success, -1 on error/disconnect.
 */
int sock_send_c(int fd, const void *buf, int nbytes)
{
    const char *p = (const char *)buf;
    int remaining = nbytes;
    while (remaining > 0) {
        ssize_t n = send(fd, p, (size_t)remaining, MSG_NOSIGNAL);
        if (n <= 0) {
            if (n < 0) perror("send");
            return -1;
        }
        p         += n;
        remaining -= (int)n;
    }
    return 0;
}

int sock_recv_c(int fd, void *buf, int nbytes)
{
    char *p = (char *)buf;
    int remaining = nbytes;
    while (remaining > 0) {
        ssize_t n = recv(fd, p, (size_t)remaining, MSG_WAITALL);
        if (n <= 0) {
            if (n < 0) perror("recv");
            return -1;
        }
        p         += n;
        remaining -= (int)n;
    }
    return 0;
}

/* ── close ───────────────────────────────────────────────────────────────── */

void sock_close_c(int fd)
{
    if (fd >= 0) close(fd);
}

/* ── process management ──────────────────────────────────────────────────── */

/**
 * Spawn a command in a new process group.
 * Returns the child PID, or -1 on error.
 * stdout_fd_out is set to the read end of a pipe connected to child stdout,
 * so the caller can wait for the "READY" line.
 */
#include <fcntl.h>
#include <signal.h>
#include <time.h>

int spawn_process_c(const char *cmd, int *stdout_fd_out)
{
    int pipefd[2];
    if (pipe(pipefd) < 0) { perror("pipe"); return -1; }

    pid_t pid = fork();
    if (pid < 0) { perror("fork"); return -1; }

    if (pid == 0) {
        /* child */
        close(pipefd[0]);                    /* close read end */
        dup2(pipefd[1], STDOUT_FILENO);      /* stdout -> pipe write end */
        close(pipefd[1]);
        /* new session so we can kill the whole group */
        setsid();
        execl("/bin/sh", "sh", "-c", cmd, (char *)NULL);
        _exit(127);
    }

    /* parent */
    close(pipefd[1]);          /* close write end */
    *stdout_fd_out = pipefd[0];
    return (int)pid;
}

/**
 * Read from fd until a line matching prefix is found (or error/EOF).
 * Returns 0 if found, -1 otherwise.
 */
int wait_for_line_c(int fd, const char *prefix, int timeout_sec)
{
    /* set non-blocking read with simple timeout loop */
    char buf[512];
    int  pos = 0;
    time_t deadline = time(NULL) + timeout_sec;
    int   plen = (int)strlen(prefix);

    while (time(NULL) < deadline) {
        ssize_t n = read(fd, buf + pos, 1);
        if (n < 0) {
            if (errno == EAGAIN || errno == EWOULDBLOCK) {
                usleep(10000); /* 10 ms */
                continue;
            }
            return -1;
        }
        if (n == 0) return -1;   /* EOF */
        if (buf[pos] == '\n') {
            buf[pos] = '\0';
            if (strncmp(buf, prefix, (size_t)plen) == 0) return 0;
            pos = 0;             /* next line */
        } else {
            pos++;
            if (pos >= 511) pos = 511; /* guard */
        }
    }
    return -1; /* timeout */
}

/**
 * Send SIGTERM to pid.
 */
void kill_process_c(int pid)
{
    if (pid > 0) kill((pid_t)pid, SIGTERM);
}


/**
 * Thin shim exposing getpid() with a plain C int return type.
 * Fortran's iso_c_binding can map this directly via c_int.
 * Cast is explicit to silence any implicit-narrowing warnings
 * (pid_t is typically a 32-bit int on modern platforms anyway). 
 */
int get_pid_shim(void) {
    return (int)getpid();
}
