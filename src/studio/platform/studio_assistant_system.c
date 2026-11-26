#ifndef _WIN32
#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE
#define _DEFAULT_SOURCE
#endif

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <io.h>
#include <wchar.h>

struct agent_worker {
    HANDLE process;
    HANDLE job;
};

static wchar_t *native_path(const char *utf8)
{
    int count = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, utf8, -1, NULL, 0);
    if (!count) return NULL;
    wchar_t *wide = malloc((size_t)count * sizeof(*wide));
    if (!wide) return NULL;
    if (!MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, utf8, -1, wide, count)) {
        free(wide);
        return NULL;
    }
    return wide;
}

void *studio_agent_start(const char *application, const char *directory)
{
    static const wchar_t arguments[] = L" --assistant-worker ";
    wchar_t *program = native_path(application);
    wchar_t *folder = native_path(directory);
    wchar_t *command = NULL;
    struct agent_worker *worker = NULL;
    JOBOBJECT_EXTENDED_LIMIT_INFORMATION job_limits = {0};
    STARTUPINFOW startup = { .cb = sizeof(startup) };
    PROCESS_INFORMATION child = {0};
    HANDLE job = NULL;

    if (!program || !folder) goto finish;
    size_t capacity = wcslen(program) + wcslen(folder) + sizeof(arguments) / sizeof(*arguments) + 4;
    command = malloc(capacity * sizeof(*command));
    if (!command) goto finish;
    if (swprintf(command, capacity, L"\"%ls\"%ls\"%ls\"", program, arguments, folder) < 0) goto finish;

    worker = malloc(sizeof(*worker));
    if (!worker) goto finish;
    job = CreateJobObjectW(NULL, NULL);
    if (!job) goto failed;
    job_limits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
    if (!SetInformationJobObject(job, JobObjectExtendedLimitInformation, &job_limits, sizeof(job_limits))) goto failed;
    if (!CreateProcessW(program, command, NULL, NULL, FALSE, CREATE_SUSPENDED | CREATE_NO_WINDOW,
                        NULL, NULL, &startup, &child)) goto failed;
    if (!AssignProcessToJobObject(job, child.hProcess)) goto failed;
    if (ResumeThread(child.hThread) == (DWORD)-1) goto failed;

    worker->process = child.hProcess;
    worker->job = job;
    CloseHandle(child.hThread);
    goto finish;

failed:
    if (child.hProcess) {
        TerminateProcess(child.hProcess, 1);
        WaitForSingleObject(child.hProcess, INFINITE);
        CloseHandle(child.hThread);
        CloseHandle(child.hProcess);
    }
    if (job) CloseHandle(job);
    free(worker);
    worker = NULL;

finish:
    free(command);
    free(folder);
    free(program);
    return worker;
}

int studio_agent_poll(void *handle, int *successful)
{
    struct agent_worker *worker = handle;
    DWORD code = 1;
    DWORD state = WaitForSingleObject(worker->process, 0);
    *successful = 0;
    if (state == WAIT_TIMEOUT) return 0;
    if (state != WAIT_OBJECT_0 || !GetExitCodeProcess(worker->process, &code)) return -1;
    *successful = (code == 0);
    return 1;
}

void studio_agent_close(void *handle)
{
    struct agent_worker *worker = handle;
    CloseHandle(worker->job);
    WaitForSingleObject(worker->process, INFINITE);
    CloseHandle(worker->process);
    free(worker);
}

int studio_append_event(const char *path, const char *frame, size_t length)
{
    if (length > UINT32_MAX) return -1;
    wchar_t *name = native_path(path);
    if (!name) return -1;
    HANDLE file = CreateFileW(name, GENERIC_READ | FILE_APPEND_DATA, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                              NULL, OPEN_EXISTING, FILE_FLAG_OPEN_REPARSE_POINT, NULL);
    free(name);
    if (file == INVALID_HANDLE_VALUE) return -1;

    BY_HANDLE_FILE_INFORMATION info;
    OVERLAPPED lock = {0};
    DWORD written = 0;
    int result = -1;
    lock.Offset = 0xfffffffe;
    lock.OffsetHigh = 0x7fffffff;
    if (GetFileInformationByHandle(file, &info) &&
        !(info.dwFileAttributes & (FILE_ATTRIBUTE_REPARSE_POINT | FILE_ATTRIBUTE_DIRECTORY)) &&
        LockFileEx(file, LOCKFILE_EXCLUSIVE_LOCK, 0, 1, 0, &lock)) {
        if (WriteFile(file, frame, (DWORD)length, &written, NULL) && written == length) result = 0;
        UnlockFileEx(file, 0, 1, 0, &lock);
    }
    if (!CloseHandle(file)) result = -1;
    return result;
}

int studio_agent_seek(int file, int64_t offset)
{
    return _lseeki64(file, offset, SEEK_SET) < 0 ? -1 : 0;
}

int64_t studio_epoch_microseconds(void)
{
    FILETIME clock;
    GetSystemTimePreciseAsFileTime(&clock);
    uint64_t ticks = ((uint64_t)clock.dwHighDateTime << 32) | clock.dwLowDateTime;
    return (int64_t)(ticks / 10) - INT64_C(11644473600000000);
}

int studio_set_environment(const char *name, const char *value)
{
    return _putenv_s(name, value);
}

int studio_unset_environment(const char *name)
{
    return _putenv_s(name, "");
}

#else

#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <spawn.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#ifdef __APPLE__
#include <crt_externs.h>
#else
extern char **environ;
#endif

struct agent_worker {
    pid_t pid;
    int reaped;
};

void *studio_agent_start(const char *program, const char *directory)
{
    struct agent_worker *worker = calloc(1, sizeof(*worker));
    if (!worker) return NULL;

    char *const arguments[] = {(char *)program, "--assistant-worker", (char *)directory, NULL};
#ifdef __APPLE__
    char **environment = *_NSGetEnviron();
#else
    char **environment = environ;
#endif
    posix_spawnattr_t attributes;
    if (posix_spawnattr_init(&attributes) != 0) {
        free(worker);
        return NULL;
    }
    int status = posix_spawnattr_setflags(&attributes, POSIX_SPAWN_SETPGROUP);
    if (status == 0) status = posix_spawnattr_setpgroup(&attributes, 0);
    if (status == 0) status = posix_spawn(&worker->pid, program, NULL, &attributes, arguments, environment);
    posix_spawnattr_destroy(&attributes);
    if (status != 0) {
        free(worker);
        return NULL;
    }
    return worker;
}

int studio_agent_poll(void *handle, int *successful)
{
    struct agent_worker *worker = handle;
    int status = 0;
    pid_t result;
    *successful = 0;
    do {
        result = waitpid(worker->pid, &status, WNOHANG);
    } while (result == -1 && errno == EINTR);
    if (result == 0) return 0;
    if (result < 0 && errno != ECHILD) return -1;

    worker->reaped = 1;
    *successful = (result > 0 && WIFEXITED(status) && WEXITSTATUS(status) == 0);
    return 1;
}

void studio_agent_close(void *handle)
{
    struct agent_worker *worker = handle;
    if (!worker->reaped) {
        int status;
        struct timespec grace = {0, 300000000L};
        if (kill(-worker->pid, SIGTERM) == -1) kill(worker->pid, SIGTERM);
        while (nanosleep(&grace, &grace) == -1 && errno == EINTR) {}
        if (kill(-worker->pid, SIGKILL) == -1) kill(worker->pid, SIGKILL);
        while (waitpid(worker->pid, &status, 0) == -1 && errno == EINTR) {}
    }
    free(worker);
}

int studio_append_event(const char *path, const char *frame, size_t length)
{
    int file = open(path, O_WRONLY | O_APPEND | O_NOFOLLOW | O_CLOEXEC | O_NONBLOCK);
    if (file < 0) return -1;

    int status = -1;
    struct stat info;
    if (fstat(file, &info) == 0 && S_ISREG(info.st_mode)) {
        do {
            status = flock(file, LOCK_EX);
        } while (status < 0 && errno == EINTR);
        if (status == 0) {
            size_t offset = 0;
            while (offset < length) {
                ssize_t written = write(file, frame + offset, length - offset);
                if (written < 0 && errno == EINTR) continue;
                if (written <= 0) break;
                offset += (size_t)written;
            }
            if (offset != length) status = -1;
        }
    }
    if (close(file) != 0) status = -1;
    return status;
}

int studio_agent_seek(int file, int64_t offset)
{
    return lseek(file, offset, SEEK_SET) < 0 ? -1 : 0;
}

int64_t studio_epoch_microseconds(void)
{
    struct timeval clock;
    if (gettimeofday(&clock, NULL) != 0) return 0;
    return (int64_t)clock.tv_sec * INT64_C(1000000) + clock.tv_usec;
}

int studio_set_environment(const char *name, const char *value)
{
    return setenv(name, value, 1);
}

int studio_unset_environment(const char *name)
{
    return unsetenv(name);
}

#endif

int studio_assistant_platform(void)
{
#if defined(_WIN32)
    return 2;
#elif defined(__APPLE__)
    return 1;
#else
    return 0;
#endif
}
