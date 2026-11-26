#ifndef _WIN32
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define _DARWIN_C_SOURCE
#endif

#include <errno.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { current_directory = -100, regular_file = 1, directory_file = 2 };

struct studio_metadata {
    int kind, owned, is_private, links;
    int64_t volume, inode, bytes;
};

#ifdef _WIN32

#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0a00
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <winternl.h>
#include <aclapi.h>
#include <bcrypt.h>
#include <fcntl.h>
#include <io.h>
#include <wchar.h>

static const DWORD share_all = FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE;

static wchar_t *wide_path(const char *path)
{
    int count = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, NULL, 0);
    if (!count) return NULL;
    wchar_t *wide = malloc((size_t)count * sizeof *wide);
    if (!wide) return NULL;
    MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, wide, count);
    for (wchar_t *c = wide; *c; ++c)
        if (*c == L'/') *c = L'\\';
    return wide;
}

static wchar_t *absolute_path(const char *path)
{
    wchar_t *wide = wide_path(path);
    if (!wide) return NULL;
    DWORD capacity = GetFullPathNameW(wide, 0, NULL, NULL);
    wchar_t *full = capacity ? malloc((capacity + 8) * sizeof *full) : NULL;
    if (!full) {
        free(wide);
        return NULL;
    }
    DWORD length = GetFullPathNameW(wide, capacity, full + 8, NULL);
    free(wide);
    if (!length || length >= capacity) {
        free(full);
        return NULL;
    }
    wchar_t *name = full + 8;
    if (wcsncmp(name, L"\\\\?\\", 4) == 0) {
        memmove(full, name, (length + 1) * sizeof *full);
    } else if (wcsncmp(name, L"\\\\.\\", 4) == 0) {
        free(full);
        return NULL;
    } else if (wcsncmp(name, L"\\\\", 2) == 0) {
        memmove(full + 8, name + 2, (length - 1) * sizeof *full);
        memcpy(full, L"\\\\?\\UNC\\", 8 * sizeof *full);
    } else {
        memmove(full + 4, name, (length + 1) * sizeof *full);
        memcpy(full, L"\\\\?\\", 4 * sizeof *full);
    }
    return full;
}

static char *utf8_path(wchar_t *name)
{
    if (wcsncmp(name, L"\\\\?\\UNC\\", 8) == 0) {
        name += 6;
        name[0] = L'\\';
        name[1] = L'\\';
    } else if (wcsncmp(name, L"\\\\?\\", 4) == 0) {
        name += 4;
    }
    for (wchar_t *c = name; *c; ++c)
        if (*c == L'\\') *c = L'/';

    int count = WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS, name, -1, NULL, 0, NULL, NULL);
    char *result = count ? malloc((size_t)count) : NULL;
    if (result) WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS, name, -1, result, count, NULL, NULL);
    return result;
}

static HANDLE handle_from_fd(int fd)
{
    return fd < 0 ? INVALID_HANDLE_VALUE : (HANDLE)_get_osfhandle(fd);
}

static HANDLE open_absolute(const char *path, DWORD access, DWORD flags, DWORD disposition, SECURITY_ATTRIBUTES *security)
{
    wchar_t *wide = absolute_path(path);
    if (!wide) return INVALID_HANDLE_VALUE;
    HANDLE handle = CreateFileW(wide, access, share_all, security, disposition, flags, NULL);
    free(wide);
    return handle;
}

typedef NTSTATUS (NTAPI *open_native_file)(PHANDLE, ACCESS_MASK, POBJECT_ATTRIBUTES, PIO_STATUS_BLOCK,
    PLARGE_INTEGER, ULONG, ULONG, ULONG, ULONG, PVOID, ULONG);

static HANDLE open_relative(int fd, const char *leaf, ACCESS_MASK access, ULONG disposition,
                            ULONG options, PSECURITY_DESCRIPTOR security)
{
    HANDLE parent = handle_from_fd(fd);
    size_t length = strlen(leaf);
    if (parent == INVALID_HANDLE_VALUE || length == 0 || length > 1024 ||
        strpbrk(leaf, "/\\:") || leaf[length - 1] == '.' || leaf[length - 1] == ' ') return INVALID_HANDLE_VALUE;

    wchar_t *wide = wide_path(leaf);
    if (!wide) return INVALID_HANDLE_VALUE;
    UNICODE_STRING name;
    name.Length = (USHORT)(wcslen(wide) * sizeof *wide);
    name.MaximumLength = (USHORT)(name.Length + sizeof *wide);
    name.Buffer = wide;

    OBJECT_ATTRIBUTES attributes = {0};
    attributes.Length = sizeof attributes;
    attributes.RootDirectory = parent;
    attributes.ObjectName = &name;
    attributes.Attributes = OBJ_CASE_INSENSITIVE | 0x1000;
    attributes.SecurityDescriptor = security;
    IO_STATUS_BLOCK io;
    HANDLE result = INVALID_HANDLE_VALUE;
    open_native_file open_file = (open_native_file)(void *)GetProcAddress(GetModuleHandleW(L"ntdll.dll"), "NtCreateFile");
    if (open_file && open_file(&result, access | SYNCHRONIZE, &attributes, &io, NULL, FILE_ATTRIBUTE_NORMAL,
                              share_all, disposition, options | FILE_SYNCHRONOUS_IO_NONALERT, NULL, 0) < 0)
        result = INVALID_HANDLE_VALUE;
    free(wide);
    return result;
}

static HANDLE open_entry(int fd, const char *path, DWORD access)
{
    if (fd == current_directory)
        return open_absolute(path, access, FILE_FLAG_OPEN_REPARSE_POINT | FILE_FLAG_BACKUP_SEMANTICS, OPEN_EXISTING, NULL);
    return open_relative(fd, path, access, FILE_OPEN, FILE_OPEN_REPARSE_POINT, NULL);
}

struct private_access {
    TOKEN_USER *user;
    PACL acl;
    SECURITY_DESCRIPTOR descriptor;
    SECURITY_ATTRIBUTES attributes;
};

static TOKEN_USER *current_user(void)
{
    HANDLE token;
    if (!OpenProcessToken(GetCurrentProcess(), TOKEN_QUERY, &token)) return NULL;
    DWORD size = 0;
    GetTokenInformation(token, TokenUser, NULL, 0, &size);
    TOKEN_USER *result = size ? malloc(size) : NULL;
    if (result && !GetTokenInformation(token, TokenUser, result, size, &size)) {
        free(result);
        result = NULL;
    }
    CloseHandle(token);
    return result;
}

static void release_access(struct private_access *access)
{
    free(access->user);
    if (access->acl) LocalFree(access->acl);
}

static int owner_access(struct private_access *access, int directory)
{
    memset(access, 0, sizeof *access);
    access->user = current_user();
    if (!access->user) return -1;
    EXPLICIT_ACCESSW entry = {0};
    entry.grfAccessPermissions = FILE_ALL_ACCESS;
    entry.grfAccessMode = SET_ACCESS;
    entry.grfInheritance = directory ? SUB_CONTAINERS_AND_OBJECTS_INHERIT : NO_INHERITANCE;
    BuildTrusteeWithSidW(&entry.Trustee, access->user->User.Sid);
    if (SetEntriesInAclW(1, &entry, NULL, &access->acl) != ERROR_SUCCESS ||
        !InitializeSecurityDescriptor(&access->descriptor, SECURITY_DESCRIPTOR_REVISION) ||
        !SetSecurityDescriptorOwner(&access->descriptor, access->user->User.Sid, FALSE) ||
        !SetSecurityDescriptorDacl(&access->descriptor, TRUE, access->acl, FALSE) ||
        !SetSecurityDescriptorControl(&access->descriptor, SE_DACL_PROTECTED, SE_DACL_PROTECTED)) return -1;

    access->attributes = (SECURITY_ATTRIBUTES){sizeof access->attributes, &access->descriptor, FALSE};
    return 0;
}

static void inspect_owner(HANDLE handle, struct studio_metadata *metadata)
{
    PSID owner = NULL;
    PACL acl = NULL;
    PSECURITY_DESCRIPTOR security = NULL;
    TOKEN_USER *user = current_user();
    if (!user) return;
    if (GetSecurityInfo(handle, SE_FILE_OBJECT, OWNER_SECURITY_INFORMATION | DACL_SECURITY_INFORMATION,
                        &owner, NULL, &acl, NULL, &security) == ERROR_SUCCESS) {
        metadata->owned = owner && EqualSid(owner, user->User.Sid);
        metadata->is_private = metadata->owned && acl && acl->AceCount > 0;
        for (DWORD i = 0; metadata->is_private && i < acl->AceCount; ++i) {
            ACCESS_ALLOWED_ACE *ace = NULL;
            if (!GetAce(acl, i, (void **)&ace) ||
                ace->Header.AceType != ACCESS_ALLOWED_ACE_TYPE ||
                !EqualSid(&ace->SidStart, user->User.Sid)) metadata->is_private = 0;
        }
    }
    if (security) LocalFree(security);
    free(user);
}

static int inspect_handle(HANDLE handle, struct studio_metadata *metadata)
{
    BY_HANDLE_FILE_INFORMATION information;
    memset(metadata, 0, sizeof *metadata);
    if (GetFileType(handle) != FILE_TYPE_DISK || !GetFileInformationByHandle(handle, &information)) return -1;
    if (!(information.dwFileAttributes & FILE_ATTRIBUTE_REPARSE_POINT)) {
        metadata->kind = information.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY ? directory_file : regular_file;
    }
    metadata->volume = information.dwVolumeSerialNumber;
    metadata->inode = (int64_t)(((uint64_t)information.nFileIndexHigh << 32) | information.nFileIndexLow);
    metadata->bytes = (int64_t)(((uint64_t)information.nFileSizeHigh << 32) | information.nFileSizeLow);
    metadata->links = (int)information.nNumberOfLinks;
    inspect_owner(handle, metadata);
    return 0;
}

static int file_descriptor(HANDLE handle, int flags)
{
    if (handle == INVALID_HANDLE_VALUE) return -1;
    int fd = _open_osfhandle((intptr_t)handle, flags | _O_BINARY | _O_NOINHERIT);
    if (fd < 0) CloseHandle(handle);
    return fd;
}

int studio_open_directory(int fd, const char *path)
{
    HANDLE handle = open_entry(fd, path, FILE_GENERIC_READ);
    BY_HANDLE_FILE_INFORMATION information;
    if (handle == INVALID_HANDLE_VALUE) return -1;
    if (!GetFileInformationByHandle(handle, &information) ||
        (information.dwFileAttributes & (FILE_ATTRIBUTE_DIRECTORY | FILE_ATTRIBUTE_REPARSE_POINT)) != FILE_ATTRIBUTE_DIRECTORY) {
        CloseHandle(handle);
        return -1;
    }
    return file_descriptor(handle, _O_RDONLY);
}

int studio_open_regular(int fd, const char *path)
{
    HANDLE handle = open_entry(fd, path, FILE_GENERIC_READ);
    BY_HANDLE_FILE_INFORMATION information;
    if (handle == INVALID_HANDLE_VALUE) return -1;
    if (GetFileType(handle) != FILE_TYPE_DISK || !GetFileInformationByHandle(handle, &information) ||
        (information.dwFileAttributes & (FILE_ATTRIBUTE_DIRECTORY | FILE_ATTRIBUTE_REPARSE_POINT))) {
        CloseHandle(handle);
        return -1;
    }
    return file_descriptor(handle, _O_RDONLY);
}

int studio_create_file(int fd, const char *path)
{
    struct private_access access;
    HANDLE handle = INVALID_HANDLE_VALUE;
    if (owner_access(&access, 0) == 0) {
        if (fd == current_directory)
            handle = open_absolute(path, FILE_GENERIC_WRITE, FILE_ATTRIBUTE_NORMAL, CREATE_NEW, &access.attributes);
        else
            handle = open_relative(fd, path, FILE_GENERIC_WRITE, FILE_CREATE, FILE_NON_DIRECTORY_FILE, &access.descriptor);
    }
    release_access(&access);
    return file_descriptor(handle, _O_WRONLY);
}

int studio_mkdir_private(int fd, const char *path)
{
    struct private_access access;
    int status = -1;
    if (owner_access(&access, 1) == 0) {
        if (fd == current_directory) {
            wchar_t *wide = absolute_path(path);
            if (wide && CreateDirectoryW(wide, &access.attributes)) status = 0;
            free(wide);
        } else {
            HANDLE handle = open_relative(fd, path, FILE_GENERIC_READ, FILE_CREATE, FILE_DIRECTORY_FILE, &access.descriptor);
            if (handle != INVALID_HANDLE_VALUE) {
                CloseHandle(handle);
                status = 0;
            }
        }
    }
    release_access(&access);
    return status;
}

int studio_inspect_fd(int fd, struct studio_metadata *metadata)
{
    return inspect_handle(handle_from_fd(fd), metadata);
}

int studio_inspect_entry(int fd, const char *path, struct studio_metadata *metadata)
{
    memset(metadata, 0, sizeof *metadata);
    HANDLE handle = open_entry(fd, path, FILE_READ_ATTRIBUTES | READ_CONTROL);
    if (handle == INVALID_HANDLE_VALUE) return -1;
    int status = inspect_handle(handle, metadata);
    CloseHandle(handle);
    return status;
}

int studio_close_fd(int fd) { return _close(fd); }

int64_t studio_read_fd(int fd, void *buffer, size_t length)
{
    return _read(fd, buffer, (unsigned int)(length > INT32_MAX ? INT32_MAX : length));
}

int64_t studio_write_fd(int fd, const void *buffer, size_t length)
{
    return _write(fd, buffer, (unsigned int)(length > INT32_MAX ? INT32_MAX : length));
}

int studio_rename_child(int source, const char *old_name, int target, const char *new_name, int replace)
{
    HANDLE parent = handle_from_fd(target);
    if (parent == INVALID_HANDLE_VALUE || strpbrk(new_name, "/\\:")) return -1;
    wchar_t *wide = wide_path(new_name);
    if (!wide) return -1;
    size_t length = wcslen(wide) * sizeof *wide;
    size_t capacity = sizeof(FILE_RENAME_INFO) + length;
    FILE_RENAME_INFO *info = calloc(1, capacity);
    HANDLE handle = open_entry(source, old_name, DELETE | FILE_READ_ATTRIBUTES);
    int status = -1;
    if (info && handle != INVALID_HANDLE_VALUE && GetFileType(handle) == FILE_TYPE_DISK) {
        info->ReplaceIfExists = replace != 0;
        info->RootDirectory = parent;
        info->FileNameLength = (DWORD)length;
        memcpy(info->FileName, wide, length);
        if (SetFileInformationByHandle(handle, FileRenameInfo, info, (DWORD)capacity)) status = 0;
    }
    if (handle != INVALID_HANDLE_VALUE) CloseHandle(handle);
    free(info);
    free(wide);
    return status;
}

int studio_remove_child(int fd, const char *name, int directory)
{
    HANDLE handle = open_entry(fd, name, FILE_READ_ATTRIBUTES | DELETE);
    if (handle == INVALID_HANDLE_VALUE) return -1;
    BY_HANDLE_FILE_INFORMATION metadata;
    FILE_DISPOSITION_INFO disposition = {TRUE};
    int status = -1;
    if (GetFileInformationByHandle(handle, &metadata) &&
        !!(metadata.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) == !!directory &&
        SetFileInformationByHandle(handle, FileDispositionInfo, &disposition, sizeof disposition)) status = 0;
    CloseHandle(handle);
    return status;
}

int studio_replace_path(const char *source, const char *target)
{
    wchar_t *old = absolute_path(source);
    wchar_t *new = absolute_path(target);
    int status = old && new && MoveFileExW(old, new, MOVEFILE_REPLACE_EXISTING) ? 0 : -1;
    free(old);
    free(new);
    return status;
}

char *studio_realpath(const char *path)
{
    HANDLE handle = open_absolute(path, FILE_READ_ATTRIBUTES, FILE_FLAG_BACKUP_SEMANTICS, OPEN_EXISTING, NULL);
    if (handle == INVALID_HANDLE_VALUE) return NULL;
    DWORD capacity = GetFinalPathNameByHandleW(handle, NULL, 0, FILE_NAME_NORMALIZED);
    wchar_t *name = capacity ? malloc((capacity + 1) * sizeof *name) : NULL;
    DWORD length = name ? GetFinalPathNameByHandleW(handle, name, capacity + 1, FILE_NAME_NORMALIZED) : 0;
    char *result = length && length <= capacity ? utf8_path(name) : NULL;
    free(name);
    CloseHandle(handle);
    return result;
}

int studio_restrict_file(const char *path)
{
    HANDLE handle = open_entry(current_directory, path, FILE_READ_ATTRIBUTES | READ_CONTROL | WRITE_DAC);
    if (handle == INVALID_HANDLE_VALUE) return -1;
    struct studio_metadata metadata;
    struct private_access access;
    int status = -1;
    if (inspect_handle(handle, &metadata) == 0 && metadata.kind == regular_file &&
        metadata.owned && metadata.links == 1) {
        if (owner_access(&access, 0) == 0 &&
            SetSecurityInfo(handle, SE_FILE_OBJECT, DACL_SECURITY_INFORMATION | PROTECTED_DACL_SECURITY_INFORMATION,
                            NULL, NULL, access.acl, NULL) == ERROR_SUCCESS) status = 0;
        release_access(&access);
    }
    CloseHandle(handle);
    return status;
}

int studio_random_fill(void *buffer, size_t length)
{
    if (length > UINT32_MAX) return -1;
    return BCryptGenRandom(NULL, buffer, (ULONG)length, BCRYPT_USE_SYSTEM_PREFERRED_RNG) == 0 ? 0 : -1;
}

#else

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef __linux__
#include <sys/random.h>
#endif
#ifdef __APPLE__
#include <sys/acl.h>

static int clear_extended_acl(int fd)
{
    acl_t empty = acl_init(0);
    if (!empty) return -1;
    int status = acl_set_fd_np(fd, empty, ACL_TYPE_EXTENDED);
    acl_free(empty);
    return status;
}
#endif

static int native_dir(int fd) { return fd == current_directory ? AT_FDCWD : fd; }

static void inspect_stat(const struct stat *info, struct studio_metadata *metadata)
{
    metadata->kind = S_ISREG(info->st_mode) ? regular_file : (S_ISDIR(info->st_mode) ? directory_file : 0);
    metadata->owned = info->st_uid == geteuid();
    metadata->is_private = metadata->owned && (info->st_mode & 077) == 0;
    metadata->links = (int)info->st_nlink;
    metadata->volume = (int64_t)info->st_dev;
    metadata->inode = (int64_t)info->st_ino;
    metadata->bytes = (int64_t)info->st_size;
}

int studio_inspect_fd(int fd, struct studio_metadata *metadata)
{
    struct stat information;
    memset(metadata, 0, sizeof *metadata);
    if (fstat(fd, &information) < 0) return -1;
    inspect_stat(&information, metadata);
#ifdef __APPLE__
    if (metadata->is_private) {
        acl_t acl = acl_get_fd_np(fd, ACL_TYPE_EXTENDED);
        if (acl) {
            acl_entry_t entry;
            if (acl_get_entry(acl, ACL_FIRST_ENTRY, &entry) == 0) metadata->is_private = 0;
            acl_free(acl);
        }
    }
#endif
    return 0;
}

int studio_inspect_entry(int fd, const char *path, struct studio_metadata *metadata)
{
    struct stat information;
    memset(metadata, 0, sizeof *metadata);
    if (fstatat(native_dir(fd), path, &information, AT_SYMLINK_NOFOLLOW) < 0) return -1;
    inspect_stat(&information, metadata);
    return 0;
}

int studio_open_directory(int fd, const char *path)
{
    return openat(native_dir(fd), path, O_RDONLY | O_DIRECTORY | O_NOFOLLOW | O_CLOEXEC);
}

int studio_open_regular(int fd, const char *path)
{
    int file = openat(native_dir(fd), path, O_RDONLY | O_NONBLOCK | O_NOFOLLOW | O_CLOEXEC);
    if (file < 0) return -1;
    struct stat information;
    if (fstat(file, &information) < 0 || !S_ISREG(information.st_mode)) {
        close(file);
        return -1;
    }
    return file;
}

int studio_create_file(int fd, const char *path)
{
    int file = openat(native_dir(fd), path, O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW | O_CLOEXEC, 0600);
#ifdef __APPLE__
    if (file >= 0 && clear_extended_acl(file) < 0) {
        close(file);
        unlinkat(native_dir(fd), path, 0);
        return -1;
    }
#endif
    return file;
}

int studio_mkdir_private(int fd, const char *path)
{
    if (mkdirat(native_dir(fd), path, 0700) < 0) return -1;
#ifdef __APPLE__
    int folder = studio_open_directory(fd, path);
    if (folder < 0 || clear_extended_acl(folder) < 0) {
        if (folder >= 0) close(folder);
        unlinkat(native_dir(fd), path, AT_REMOVEDIR);
        return -1;
    }
    close(folder);
#endif
    return 0;
}

int studio_close_fd(int fd) { return close(fd); }

int64_t studio_read_fd(int fd, void *buffer, size_t length)
{
    ssize_t count;
    do { count = read(fd, buffer, length); } while (count < 0 && errno == EINTR);
    return count;
}

int64_t studio_write_fd(int fd, const void *buffer, size_t length)
{
    ssize_t count;
    do { count = write(fd, buffer, length); } while (count < 0 && errno == EINTR);
    return count;
}

int studio_rename_child(int source, const char *old_name, int target, const char *new_name, int replace)
{
    if (replace) return renameat(native_dir(source), old_name, native_dir(target), new_name);
#ifdef __APPLE__
    return renameatx_np(native_dir(source), old_name, native_dir(target), new_name, RENAME_EXCL);
#elif defined(__linux__)
    return renameat2(native_dir(source), old_name, native_dir(target), new_name, RENAME_NOREPLACE);
#endif
}

int studio_remove_child(int fd, const char *name, int directory)
{
    return unlinkat(native_dir(fd), name, directory ? AT_REMOVEDIR : 0);
}

int studio_replace_path(const char *source, const char *target) { return rename(source, target); }

char *studio_realpath(const char *path) { return realpath(path, NULL); }

int studio_restrict_file(const char *path)
{
    int fd = studio_open_regular(current_directory, path);
    if (fd < 0) return -1;
    struct studio_metadata metadata;
    int status = -1;
    if (studio_inspect_fd(fd, &metadata) == 0 && metadata.owned && metadata.links == 1) {
        status = fchmod(fd, 0600);
#ifdef __APPLE__
        if (status == 0 && clear_extended_acl(fd) < 0) status = -1;
#endif
    }
    close(fd);
    return status;
}

int studio_random_fill(void *buffer, size_t length)
{
#ifdef __APPLE__
    arc4random_buf(buffer, length);
#else
    size_t offset = 0;
    while (offset < length) {
        ssize_t count = getrandom((char *)buffer + offset, length - offset, 0);
        if (count < 0 && errno == EINTR) continue;
        if (count <= 0) return -1;
        offset += (size_t)count;
    }
#endif
    return 0;
}

#endif

int studio_is_private_directory(const char *path)
{
    int fd = studio_open_directory(current_directory, path);
    if (fd < 0) return 0;
    struct studio_metadata metadata;
    int result = studio_inspect_fd(fd, &metadata) == 0 && metadata.is_private;
    studio_close_fd(fd);
    return result;
}

int studio_fill_suffix(char *name, size_t length)
{
    const char alphabet[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    unsigned char bytes[32];
    size_t written = 0;
    while (written < length) {
        if (studio_random_fill(bytes, sizeof bytes) < 0) return -1;
        for (size_t i = 0; i < sizeof bytes && written < length; ++i)
            if (bytes[i] < 248) name[written++] = alphabet[bytes[i] % (sizeof alphabet - 1)];
    }
    return 0;
}

static int copy_file(int input, int output)
{
    char buffer[65536];
    for (;;) {
        int64_t read = studio_read_fd(input, buffer, sizeof buffer);
        if (read <= 0) return read == 0 ? 0 : -1;
        int64_t offset = 0;
        while (offset < read) {
            int64_t written = studio_write_fd(output, buffer + offset, (size_t)(read - offset));
            if (written <= 0) return -1;
            offset += written;
        }
    }
}

int studio_copy_from_fd(int input, const char *path)
{
    int output = studio_create_file(current_directory, path);
    if (output < 0) return -1;
    int status = copy_file(input, output);
    if (studio_close_fd(output) < 0) status = -1;
    if (status < 0) studio_remove_child(current_directory, path, 0);
    return status;
}

int studio_publish_file(const char *source, int parent, const char *target, int replace)
{
    int input = studio_open_regular(current_directory, source);
    if (input < 0) return -1;
    char staged[] = ".forcad-publish.0000000000000000";
    int output = -1;
    if (studio_fill_suffix(staged + sizeof ".forcad-publish." - 1, 16) == 0)
        output = studio_create_file(parent, staged);
    int status = -1;
    if (output >= 0) {
        status = copy_file(input, output);
        if (studio_close_fd(output) < 0) status = -1;
        if (status == 0) status = studio_rename_child(parent, staged, parent, target, replace);
        if (status < 0) studio_remove_child(parent, staged, 0);
    }
    studio_close_fd(input);
    return status;
}

size_t studio_path_separator(const char *path, size_t length)
{
    while (length > 0) {
        if (path[length - 1] == '/') return length;
#ifdef _WIN32
        if (path[length - 1] == '\\') return length;
#endif
        --length;
    }
    return 0;
}
