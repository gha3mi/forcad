#ifndef _WIN32
#define _POSIX_C_SOURCE 200809L
#endif

#include <signal.h>
#include <stdio.h>

static volatile sig_atomic_t interrupted;

static void interrupt_stdio(int number)
{
    interrupted = number;
}

int studio_ipc_interrupted(void)
{
    return (int)interrupted;
}

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <conio.h>
#include <fcntl.h>
#include <io.h>

static void (*previous_interrupt)(int);
static void (*previous_terminate)(int);
static int previous_output_mode;

void studio_ipc_begin(void)
{
    interrupted = 0;
    previous_interrupt = signal(SIGINT, interrupt_stdio);
    previous_terminate = signal(SIGTERM, interrupt_stdio);
    previous_output_mode = _setmode(_fileno(stdout), _O_BINARY);
}

void studio_ipc_end(void)
{
    if (previous_interrupt != SIG_ERR)
        signal(SIGINT, previous_interrupt);
    if (previous_terminate != SIG_ERR)
        signal(SIGTERM, previous_terminate);
    if (previous_output_mode != -1)
        _setmode(_fileno(stdout), previous_output_mode);
}

int studio_ipc_read(char *buffer, int capacity, int milliseconds)
{
    HANDLE input = GetStdHandle(STD_INPUT_HANDLE);
    DWORD input_type, available, received;
    ULONGLONG deadline;
    int character;

    if (input == NULL || input == INVALID_HANDLE_VALUE)
        return -1;
    input_type = GetFileType(input);
    if (input_type == FILE_TYPE_UNKNOWN)
        return -1;
    deadline = GetTickCount64() + (DWORD)milliseconds;
    for (;;) {
        if (interrupted)
            return 0;
        if (input_type == FILE_TYPE_PIPE) {
            if (!PeekNamedPipe(input, NULL, 0, NULL, &available, NULL))
                return -1;
            if (available > 0) {
                if (available > (DWORD)capacity)
                    available = (DWORD)capacity;
                if (!ReadFile(input, buffer, available, &received, NULL) && GetLastError() != ERROR_MORE_DATA)
                    return -1;
                return received == 0 ? -1 : (int)received;
            }
        } else if (input_type == FILE_TYPE_DISK) {
            if (!ReadFile(input, buffer, (DWORD)capacity, &received, NULL))
                return -1;
            return received == 0 ? -1 : (int)received;
        } else {
            DWORD console_mode;
            if (!GetConsoleMode(input, &console_mode))
                return -1;
            if (_kbhit()) {
                character = _getch();
                if (character == 26)
                    return -1;
                buffer[0] = (char)(character == '\r' ? '\n' : character);
                return 1;
            }
        }
        ULONGLONG now = GetTickCount64();
        if (now >= deadline)
            return 0;
        Sleep(deadline - now < 10 ? (DWORD)(deadline - now) : 10);
    }
}
#else
#include <errno.h>
#include <poll.h>
#include <unistd.h>

static struct sigaction previous_interrupt;
static struct sigaction previous_terminate;
static struct sigaction previous_pipe;
static int restore_interrupt;
static int restore_terminate;
static int restore_pipe;

void studio_ipc_begin(void)
{
    struct sigaction action = {0};

    interrupted = 0;
    action.sa_handler = interrupt_stdio;
    sigemptyset(&action.sa_mask);
    restore_interrupt = sigaction(SIGINT, &action, &previous_interrupt) == 0;
    restore_terminate = sigaction(SIGTERM, &action, &previous_terminate) == 0;
    action.sa_handler = SIG_IGN;
    restore_pipe = sigaction(SIGPIPE, &action, &previous_pipe) == 0;
}

void studio_ipc_end(void)
{
    if (restore_interrupt)
        sigaction(SIGINT, &previous_interrupt, NULL);
    if (restore_terminate)
        sigaction(SIGTERM, &previous_terminate, NULL);
    if (restore_pipe)
        sigaction(SIGPIPE, &previous_pipe, NULL);
}

int studio_ipc_read(char *buffer, int capacity, int milliseconds)
{
    struct pollfd input = {STDIN_FILENO, POLLIN, 0};
    int ready;
    ssize_t received;

    if (interrupted)
        return 0;
    ready = poll(&input, 1, milliseconds);
    if (ready < 0)
        return errno == EINTR ? 0 : -1;
    if (ready == 0)
        return 0;
    if (input.revents & POLLNVAL)
        return -1;
    received = read(input.fd, buffer, (size_t)capacity);
    if (received < 0 && (errno == EINTR || errno == EAGAIN || errno == EWOULDBLOCK))
        return 0;
    return received <= 0 ? -1 : (int)received;
}
#endif
