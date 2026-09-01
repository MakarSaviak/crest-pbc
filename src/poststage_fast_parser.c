/*
 * Fast, fail-closed parser for canonical process-MTD XYZ trajectories.
 *
 * The legacy path performs formatted Fortran I/O and an internal list-directed
 * READ for every coordinate record.  Under many OpenMP parser workers the
 * Fortran runtime becomes the dominant post-MTD cost.  This implementation
 * reads each completed worker trajectory into one private byte buffer, parses
 * disjoint files without shared parser state, and writes directly into the
 * caller-owned final coordinate slice.
 */
#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

static const char *const element_symbols[119] = {
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};

enum {
    CREST_PARSE_OK = 0,
    CREST_PARSE_ARGS = 1,
    CREST_PARSE_IO = 2,
    CREST_PARSE_INPUT = 3,
    CREST_PARSE_MEMORY = 4
};

static double monotonic_seconds(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
        return 0.0;
    }
    return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

static void set_message(char *dst, int dst_len, const char *text) {
    if (dst == NULL || dst_len <= 0) {
        return;
    }
    if (text == NULL) {
        dst[0] = '\0';
        return;
    }
    (void)snprintf(dst, (size_t)dst_len, "%s", text);
    dst[dst_len - 1] = '\0';
}

static void set_context_message(char *dst, int dst_len, const char *what,
                                int frame, int atom) {
    if (dst == NULL || dst_len <= 0) {
        return;
    }
    if (atom > 0) {
        (void)snprintf(dst, (size_t)dst_len, "%s at frame %d, atom %d",
                       what, frame, atom);
    } else if (frame > 0) {
        (void)snprintf(dst, (size_t)dst_len, "%s at frame %d", what, frame);
    } else {
        (void)snprintf(dst, (size_t)dst_len, "%s", what);
    }
    dst[dst_len - 1] = '\0';
}

static const char *skip_horizontal_space(const char *p, const char *end) {
    while (p < end && (*p == ' ' || *p == '\t')) {
        ++p;
    }
    return p;
}

static const char *skip_all_space(const char *p, const char *end) {
    while (p < end && isspace((unsigned char)*p)) {
        ++p;
    }
    return p;
}

static int next_line(const char **cursor, const char *end,
                     const char **line_begin, const char **line_end) {
    const char *p;
    const char *newline;

    if (cursor == NULL || *cursor == NULL || *cursor >= end) {
        return 0;
    }
    p = *cursor;
    newline = memchr(p, '\n', (size_t)(end - p));
    if (newline != NULL) {
        *line_begin = p;
        *line_end = newline;
        *cursor = newline + 1;
    } else {
        *line_begin = p;
        *line_end = end;
        *cursor = end;
    }
    if (*line_end > *line_begin && (*line_end)[-1] == '\r') {
        --(*line_end);
    }
    return 1;
}

static int parse_positive_integer_line(const char *begin, const char *end,
                                       int *value) {
    char *after = NULL;
    long parsed;
    const char *p = skip_horizontal_space(begin, end);

    if (p >= end) {
        return 0;
    }
    errno = 0;
    parsed = strtol(p, &after, 10);
    if (errno != 0 || after == p || after > end || parsed < 1 ||
        parsed > INT32_MAX) {
        return 0;
    }
    p = skip_horizontal_space(after, end);
    if (p != end) {
        return 0;
    }
    *value = (int)parsed;
    return 1;
}

static int parse_finite_double(const char **cursor, const char *line_end,
                               double *value) {
    char *after = NULL;
    const char *p = skip_horizontal_space(*cursor, line_end);
    double parsed;

    if (p >= line_end) {
        return 0;
    }
    errno = 0;
    parsed = strtod(p, &after);
    if (after == p || after > line_end || errno == ERANGE || !isfinite(parsed)) {
        return 0;
    }
    *cursor = after;
    *value = parsed;
    return 1;
}

static int parse_comment_line(const char *begin, const char *end,
                              double *energy) {
    const char *p = skip_horizontal_space(begin, end);

    if ((size_t)(end - p) < 4U || memcmp(p, "Epot", 4U) != 0) {
        return 0;
    }
    p += 4;
    if (p < end && !isspace((unsigned char)*p) && *p != '=') {
        return 0;
    }
    p = skip_horizontal_space(p, end);
    if (p >= end || *p != '=') {
        return 0;
    }
    ++p;
    if (!parse_finite_double(&p, end, energy)) {
        return 0;
    }
    /* The legacy comment validation reads only Epot, '=', and the first
       finite number; additional comment fields are therefore accepted. */
    return 1;
}

static int parse_symbol(const char **cursor, const char *line_end,
                        int atomic_number) {
    const char *p = skip_horizontal_space(*cursor, line_end);
    const char *start = p;
    const char *expected;
    size_t token_len;
    size_t expected_len;

    if (atomic_number < 1 || atomic_number > 118) {
        return 0;
    }
    while (p < line_end && !isspace((unsigned char)*p)) {
        ++p;
    }
    token_len = (size_t)(p - start);
    expected = element_symbols[atomic_number];
    expected_len = strlen(expected);
    if (token_len != expected_len || memcmp(start, expected, expected_len) != 0) {
        return 0;
    }
    *cursor = p;
    return 1;
}

static void store_comment(char *comments, int comment_len, int frame_index,
                          const char *begin, const char *end) {
    char *dst = comments + (size_t)frame_index * (size_t)comment_len;
    size_t line_len;
    size_t copy_len;

    while (end > begin && end[-1] == ' ') {
        --end;
    }
    line_len = (size_t)(end - begin);
    copy_len = line_len < (size_t)comment_len ? line_len : (size_t)comment_len;
    memset(dst, ' ', (size_t)comment_len);
    if (copy_len > 0U) {
        memcpy(dst, begin, copy_len);
    }
}

int crest_read_canonical_xyz_fast(const char *path, const int *expected_at,
                                  int nat, int expected_frames, double *xyz,
                                  double *energies, char *comments,
                                  int comment_len, int64_t *file_bytes,
                                  double *read_seconds, double *decode_seconds,
                                  int *error_frame, int *error_atom,
                                  char *message, int message_len) {
    int fd = -1;
    struct stat st;
    char *buffer = NULL;
    size_t file_size;
    size_t received = 0U;
    double read_start;
    double read_finish;
    double parse_start;
    double parse_finish;
    const char *cursor;
    const char *end;
    const char *line_begin;
    const char *line_end;
    int frame;
    int atom;
    int natoms;
    int rc = CREST_PARSE_OK;

    if (file_bytes != NULL) {
        *file_bytes = 0;
    }
    if (read_seconds != NULL) {
        *read_seconds = 0.0;
    }
    if (decode_seconds != NULL) {
        *decode_seconds = 0.0;
    }
    if (error_frame != NULL) {
        *error_frame = 0;
    }
    if (error_atom != NULL) {
        *error_atom = 0;
    }
    set_message(message, message_len, "");

    if (path == NULL || path[0] == '\0' || expected_at == NULL || nat < 1 ||
        expected_frames < 1 || xyz == NULL || energies == NULL ||
        comments == NULL || comment_len < 1 || file_bytes == NULL) {
        set_message(message, message_len, "invalid fast trajectory parser arguments");
        return CREST_PARSE_ARGS;
    }
    for (atom = 0; atom < nat; ++atom) {
        if (expected_at[atom] < 1 || expected_at[atom] > 118) {
            set_message(message, message_len,
                        "expected atom list is outside atomic numbers 1:118");
            return CREST_PARSE_ARGS;
        }
    }

    read_start = monotonic_seconds();
    fd = open(path, O_RDONLY | O_CLOEXEC);
    if (fd < 0) {
        set_message(message, message_len, "cannot open canonical trajectory");
        return CREST_PARSE_IO;
    }
    if (fstat(fd, &st) != 0 || st.st_size < 1) {
        set_message(message, message_len,
                    "cannot inspect or empty canonical trajectory");
        rc = CREST_PARSE_IO;
        goto cleanup;
    }
    if ((uintmax_t)st.st_size > (uintmax_t)(SIZE_MAX - 1U)) {
        set_message(message, message_len,
                    "canonical trajectory exceeds addressable buffer size");
        rc = CREST_PARSE_MEMORY;
        goto cleanup;
    }
    file_size = (size_t)st.st_size;
    buffer = (char *)malloc(file_size + 1U);
    if (buffer == NULL) {
        set_message(message, message_len,
                    "cannot allocate canonical trajectory byte buffer");
        rc = CREST_PARSE_MEMORY;
        goto cleanup;
    }
    while (received < file_size) {
        ssize_t got = read(fd, buffer + received, file_size - received);
        if (got < 0) {
            if (errno == EINTR) {
                continue;
            }
            set_message(message, message_len,
                        "cannot read canonical trajectory byte buffer");
            rc = CREST_PARSE_IO;
            goto cleanup;
        }
        if (got == 0) {
            set_message(message, message_len,
                        "canonical trajectory ended during byte-buffer read");
            rc = CREST_PARSE_IO;
            goto cleanup;
        }
        received += (size_t)got;
    }
    buffer[file_size] = '\0';
    *file_bytes = (int64_t)file_size;
    read_finish = monotonic_seconds();
    if (read_seconds != NULL) {
        *read_seconds = read_finish - read_start;
    }

    parse_start = monotonic_seconds();
    cursor = buffer;
    end = buffer + file_size;
    for (frame = 0; frame < expected_frames; ++frame) {
        double energy;

        if (error_frame != NULL) {
            *error_frame = frame + 1;
        }
        if (error_atom != NULL) {
            *error_atom = 0;
        }
        if (!next_line(&cursor, end, &line_begin, &line_end) ||
            !parse_positive_integer_line(line_begin, line_end, &natoms)) {
            set_context_message(message, message_len,
                                "malformed or missing atom count", frame + 1, 0);
            rc = CREST_PARSE_INPUT;
            goto cleanup;
        }
        if (natoms != nat) {
            set_context_message(message, message_len,
                                "canonical frame has unexpected atom count",
                                frame + 1, 0);
            rc = CREST_PARSE_INPUT;
            goto cleanup;
        }
        if (!next_line(&cursor, end, &line_begin, &line_end) ||
            !parse_comment_line(line_begin, line_end, &energy)) {
            set_context_message(message, message_len,
                                "noncanonical Epot comment", frame + 1, 0);
            rc = CREST_PARSE_INPUT;
            goto cleanup;
        }
        energies[frame] = energy;
        store_comment(comments, comment_len, frame, line_begin, line_end);

        for (atom = 0; atom < nat; ++atom) {
            const char *p;
            double x;
            double y;
            double z;
            size_t base;

            if (error_atom != NULL) {
                *error_atom = atom + 1;
            }
            if (!next_line(&cursor, end, &line_begin, &line_end)) {
                set_context_message(message, message_len,
                                    "missing canonical coordinate",
                                    frame + 1, atom + 1);
                rc = CREST_PARSE_INPUT;
                goto cleanup;
            }
            p = line_begin;
            if (!parse_symbol(&p, line_end, expected_at[atom]) ||
                !parse_finite_double(&p, line_end, &x) ||
                !parse_finite_double(&p, line_end, &y) ||
                !parse_finite_double(&p, line_end, &z)) {
                set_context_message(message, message_len,
                                    "noncanonical coordinate",
                                    frame + 1, atom + 1);
                rc = CREST_PARSE_INPUT;
                goto cleanup;
            }
            p = skip_horizontal_space(p, line_end);
            if (p != line_end) {
                set_context_message(message, message_len,
                                    "trailing token in canonical coordinate",
                                    frame + 1, atom + 1);
                rc = CREST_PARSE_INPUT;
                goto cleanup;
            }
            base = 3U * ((size_t)atom + (size_t)nat * (size_t)frame);
            xyz[base] = x;
            xyz[base + 1U] = y;
            xyz[base + 2U] = z;
        }
    }

    cursor = skip_all_space(cursor, end);
    if (cursor != end) {
        set_message(message, message_len,
                    "canonical trajectory contains extra trailing content");
        rc = CREST_PARSE_INPUT;
        goto cleanup;
    }
    parse_finish = monotonic_seconds();
    if (decode_seconds != NULL) {
        *decode_seconds = parse_finish - parse_start;
    }
    if (error_frame != NULL) {
        *error_frame = 0;
    }
    if (error_atom != NULL) {
        *error_atom = 0;
    }

cleanup:
    if (buffer != NULL) {
        free(buffer);
    }
    if (fd >= 0) {
        (void)close(fd);
    }
    return rc;
}

int crest_write_canonical_xyz_fast(const char *path,
                                   const int *atomic_numbers,
                                   int nat, int nframes,
                                   const double *xyz,
                                   const char *comments,
                                   int comment_len,
                                   int64_t *file_bytes,
                                   double *write_seconds,
                                   int *error_frame,
                                   int *error_atom,
                                   char *message, int message_len) {
    FILE *stream = NULL;
    char *io_buffer = NULL;
    const size_t io_buffer_size = 16U * 1024U * 1024U;
    double start;
    double finish;
    int frame;
    int atom;
    int rc = CREST_PARSE_OK;

    if (file_bytes != NULL) {
        *file_bytes = 0;
    }
    if (write_seconds != NULL) {
        *write_seconds = 0.0;
    }
    if (error_frame != NULL) {
        *error_frame = 0;
    }
    if (error_atom != NULL) {
        *error_atom = 0;
    }
    set_message(message, message_len, "");

    if (path == NULL || path[0] == '\0' || atomic_numbers == NULL ||
        xyz == NULL || comments == NULL || nat < 1 || nframes < 1 ||
        comment_len < 1) {
        set_message(message, message_len,
                    "invalid fast canonical writer arguments");
        return CREST_PARSE_ARGS;
    }
    for (atom = 0; atom < nat; ++atom) {
        if (atomic_numbers[atom] < 1 || atomic_numbers[atom] > 118) {
            if (error_atom != NULL) {
                *error_atom = atom + 1;
            }
            set_context_message(message, message_len,
                                "invalid atomic number", 0, atom + 1);
            return CREST_PARSE_ARGS;
        }
    }

    stream = fopen(path, "wb");
    if (stream == NULL) {
        set_message(message, message_len, strerror(errno));
        return CREST_PARSE_IO;
    }
    io_buffer = malloc(io_buffer_size);
    if (io_buffer != NULL) {
        (void)setvbuf(stream, io_buffer, _IOFBF, io_buffer_size);
    }

    start = monotonic_seconds();
    for (frame = 0; frame < nframes; ++frame) {
        const char *comment = comments + (size_t)frame * (size_t)comment_len;
        int comment_size = comment_len;
        if (fprintf(stream, "  %d\n", nat) < 0) {
            rc = CREST_PARSE_IO;
            if (error_frame != NULL) {
                *error_frame = frame + 1;
            }
            set_context_message(message, message_len,
                                "cannot write atom-count record", frame + 1, 0);
            break;
        }
        while (comment_size > 0 && comment[comment_size - 1] == ' ') {
            --comment_size;
        }
        if ((comment_size > 0 &&
             fwrite(comment, 1U, (size_t)comment_size, stream) !=
                 (size_t)comment_size) ||
            fputc('\n', stream) == EOF) {
            rc = CREST_PARSE_IO;
            if (error_frame != NULL) {
                *error_frame = frame + 1;
            }
            set_context_message(message, message_len,
                                "cannot write comment record", frame + 1, 0);
            break;
        }
        for (atom = 0; atom < nat; ++atom) {
            const size_t offset =
                ((size_t)frame * (size_t)nat + (size_t)atom) * 3U;
            const double x = xyz[offset];
            const double y = xyz[offset + 1U];
            const double z = xyz[offset + 2U];
            if (!isfinite(x) || !isfinite(y) || !isfinite(z)) {
                rc = CREST_PARSE_INPUT;
                if (error_frame != NULL) {
                    *error_frame = frame + 1;
                }
                if (error_atom != NULL) {
                    *error_atom = atom + 1;
                }
                set_context_message(message, message_len,
                                    "non-finite coordinate", frame + 1,
                                    atom + 1);
                break;
            }
            if (fprintf(stream, " %-2s %20.10f%20.10f%20.10f\n",
                        element_symbols[atomic_numbers[atom]], x, y, z) < 0) {
                rc = CREST_PARSE_IO;
                if (error_frame != NULL) {
                    *error_frame = frame + 1;
                }
                if (error_atom != NULL) {
                    *error_atom = atom + 1;
                }
                set_context_message(message, message_len,
                                    "cannot write coordinate record",
                                    frame + 1, atom + 1);
                break;
            }
        }
        if (rc != CREST_PARSE_OK) {
            break;
        }
    }

    if (rc == CREST_PARSE_OK && fflush(stream) != 0) {
        rc = CREST_PARSE_IO;
        set_message(message, message_len, strerror(errno));
    }
    if (rc == CREST_PARSE_OK) {
        const long position = ftell(stream);
        if (position < 0L) {
            rc = CREST_PARSE_IO;
            set_message(message, message_len, strerror(errno));
        } else if (file_bytes != NULL) {
            *file_bytes = (int64_t)position;
        }
    }
    if (fclose(stream) != 0 && rc == CREST_PARSE_OK) {
        rc = CREST_PARSE_IO;
        set_message(message, message_len, strerror(errno));
    }
    stream = NULL;
    finish = monotonic_seconds();
    if (write_seconds != NULL) {
        *write_seconds = finish - start;
    }
    free(io_buffer);
    if (rc != CREST_PARSE_OK) {
        (void)unlink(path);
    }
    return rc;
}
