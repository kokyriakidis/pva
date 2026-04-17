// pva_utils.h — internal helpers shared across pva subcommand implementations.
// Not part of the public API; do not install or expose this header.

#pragma once

#include <string>
#include <vector>
#include <cstring>
#include <cerrno>
#include <iostream>

#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <linux/limits.h>

namespace pva::utils {

// ── Binary discovery ──────────────────────────────────────────────────────────
// Returns the path to a binary. Preference order:
//   1. hint (caller-supplied, if non-empty)
//   2. <dir_of_running_pva_binary>/<relative_to_self>
//   3. path_fallback (looked up via PATH by execvp)
inline std::string find_bin(const std::string& hint,
                             const std::string& relative_to_self,
                             const std::string& path_fallback)
{
    if (!hint.empty()) return hint;

    char self[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", self, sizeof(self) - 1);
    if (len > 0) {
        self[len] = '\0';
        std::string dir(self);
        auto slash = dir.rfind('/');
        if (slash != std::string::npos) {
            std::string candidate = dir.substr(0, slash + 1) + relative_to_self;
            if (access(candidate.c_str(), X_OK) == 0)
                return candidate;
        }
    }
    return path_fallback;
}

// ── Fork + exec, redirect stdout to a file ────────────────────────────────────
// Returns the child's exit code, or -1 on system error.
inline int run_to_file(const std::vector<std::string>& args,
                       const std::string&               out_path)
{
    std::vector<const char*> argv;
    for (auto& s : args) argv.push_back(s.c_str());
    argv.push_back(nullptr);

    int fd = open(out_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) {
        std::cerr << "pva: cannot open " << out_path << ": " << strerror(errno) << "\n";
        return -1;
    }

    pid_t pid = fork();
    if (pid < 0) {
        std::cerr << "pva: fork failed: " << strerror(errno) << "\n";
        close(fd);
        return -1;
    }
    if (pid == 0) {
        if (dup2(fd, STDOUT_FILENO) < 0) _exit(1);
        close(fd);
        execvp(argv[0], const_cast<char* const*>(argv.data()));
        std::cerr << "pva: exec failed (" << args[0] << "): " << strerror(errno) << "\n";
        _exit(127);
    }
    close(fd);

    int status = 0;
    waitpid(pid, &status, 0);
    if (!WIFEXITED(status)) return -1;
    return WEXITSTATUS(status);
}

// ── Fork + exec, capture stdout into a string ─────────────────────────────────
// Returns the child's exit code, or -1 on system error. Output appended to out.
inline int run_capture(const std::vector<std::string>& args, std::string& out)
{
    std::vector<const char*> argv;
    for (auto& s : args) argv.push_back(s.c_str());
    argv.push_back(nullptr);

    int pipefd[2];
    if (pipe(pipefd) < 0) {
        std::cerr << "pva: pipe failed: " << strerror(errno) << "\n";
        return -1;
    }

    pid_t pid = fork();
    if (pid < 0) {
        std::cerr << "pva: fork failed: " << strerror(errno) << "\n";
        close(pipefd[0]); close(pipefd[1]);
        return -1;
    }
    if (pid == 0) {
        close(pipefd[0]);
        if (dup2(pipefd[1], STDOUT_FILENO) < 0) _exit(1);
        close(pipefd[1]);
        execvp(argv[0], const_cast<char* const*>(argv.data()));
        std::cerr << "pva: exec failed (" << args[0] << "): " << strerror(errno) << "\n";
        _exit(127);
    }
    close(pipefd[1]);

    char buf[4096];
    ssize_t n;
    while ((n = read(pipefd[0], buf, sizeof(buf))) > 0)
        out.append(buf, n);
    close(pipefd[0]);

    int status = 0;
    waitpid(pid, &status, 0);
    if (!WIFEXITED(status)) return -1;
    return WEXITSTATUS(status);
}

// ── Fork + exec, pipe a string as stdin, redirect stdout to a file ────────────
// Used to re-encode filtered NDJSON snarls back to binary via `vg view -JR -`.
// Returns the child's exit code, or -1 on system error.
inline int run_pipe_to_file(const std::vector<std::string>& args,
                             const std::string&               input,
                             const std::string&               out_path)
{
    std::vector<const char*> argv;
    for (auto& s : args) argv.push_back(s.c_str());
    argv.push_back(nullptr);

    int in_pipe[2];
    if (pipe(in_pipe) < 0) {
        std::cerr << "pva: pipe failed: " << strerror(errno) << "\n";
        return -1;
    }

    int fd = open(out_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) {
        std::cerr << "pva: cannot open " << out_path << ": " << strerror(errno) << "\n";
        close(in_pipe[0]); close(in_pipe[1]);
        return -1;
    }

    pid_t pid = fork();
    if (pid < 0) {
        std::cerr << "pva: fork failed: " << strerror(errno) << "\n";
        close(in_pipe[0]); close(in_pipe[1]); close(fd);
        return -1;
    }
    if (pid == 0) {
        close(in_pipe[1]);
        if (dup2(in_pipe[0], STDIN_FILENO) < 0) _exit(1);
        close(in_pipe[0]);
        if (dup2(fd, STDOUT_FILENO) < 0) _exit(1);
        close(fd);
        execvp(argv[0], const_cast<char* const*>(argv.data()));
        std::cerr << "pva: exec failed (" << args[0] << "): " << strerror(errno) << "\n";
        _exit(127);
    }
    close(in_pipe[0]);
    close(fd);

    // Write input to child's stdin in chunks
    const char* p = input.data();
    size_t remaining = input.size();
    while (remaining > 0) {
        ssize_t n = write(in_pipe[1], p, std::min(remaining, (size_t)65536));
        if (n <= 0) break;
        p += n;
        remaining -= n;
    }
    close(in_pipe[1]); // signal EOF to child

    int status = 0;
    waitpid(pid, &status, 0);
    if (!WIFEXITED(status)) return -1;
    return WEXITSTATUS(status);
}

} // namespace pva::utils
