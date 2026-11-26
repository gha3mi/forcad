#include <stdint.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#elif defined(__linux__)
#include <inttypes.h>
#include <stdio.h>
#include <string.h>

static int64_t memory_quantity(const char *group, const char *field)
{
    char path[8192];
    int length = snprintf(path, sizeof path, "/sys/fs/cgroup%s/%s", group, field);
    if (length < 0 || (size_t)length >= sizeof path) return -1;
    FILE *input = fopen(path, "r");
    if (!input) return -1;
    int64_t amount = -1;
    if (fscanf(input, "%" SCNd64, &amount) != 1) amount = -1;
    fclose(input);
    return amount;
}

static int64_t linux_available_memory(void)
{
    char line[4096];
    int64_t bytes = -1;
    FILE *input = fopen("/proc/meminfo", "r");
    if (!input) return -1;
    while (fgets(line, sizeof line, input)) {
        int64_t kib;
        if (sscanf(line, "MemAvailable: %" SCNd64, &kib) == 1 && kib >= 0 && kib <= INT64_MAX / 1024) {
            bytes = kib * 1024;
            break;
        }
    }
    fclose(input);

    input = fopen("/proc/self/cgroup", "r");
    if (!input) return bytes;
    while (fgets(line, sizeof line, input)) {
        if (strncmp(line, "0::", 3) != 0) continue;
        char *group = line + 3;
        group[strcspn(group, "\r\n")] = '\0';
        for (;;) {
            int64_t limit = memory_quantity(group, "memory.max");
            if (limit >= 0) {
                int64_t used = memory_quantity(group, "memory.current");
                if (used >= 0) {
                    int64_t remaining = limit > used ? limit - used : 0;
                    if (bytes < 0) bytes = INT64_C(1024) * 1024 * 1024;
                    if (remaining < bytes) bytes = remaining;
                }
            }
            char *slash = strrchr(group, '/');
            if (!slash) break;
            *slash = '\0';
        }
        break;
    }
    fclose(input);
    return bytes;
}
#endif

int64_t studio_available_memory(void)
{
#ifdef _WIN32
    MEMORYSTATUSEX status = {0};
    status.dwLength = sizeof status;
    if (!GlobalMemoryStatusEx(&status)) return -1;
    return (int64_t)(status.ullAvailPhys < status.ullAvailPageFile ? status.ullAvailPhys : status.ullAvailPageFile);
#elif defined(__APPLE__)
    vm_statistics64_data_t stats;
    vm_size_t page_size;
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
    mach_port_t host = mach_host_self();
    kern_return_t result = host_page_size(host, &page_size);
    if (result == KERN_SUCCESS)
        result = host_statistics64(host, HOST_VM_INFO64, (host_info64_t)&stats, &count);
    mach_port_deallocate(mach_task_self(), host);
    if (result != KERN_SUCCESS) return -1;
    return ((int64_t)stats.free_count + stats.inactive_count) * page_size;
#elif defined(__linux__)
    return linux_available_memory();
#else
    return -1;
#endif
}
