/* Exact cyclic Jacobi spectra through compiled libpari, with no GP interpreter.
 *
 * The mathematics is the same classical transfer/Chebyshev reduction as in
 * proof.md and pari.gp.  Both sequential and binary-doubling recurrences use
 * PARI's specialized integer-polynomial kernels, including ZX_mul and ZX_sqr.
 *
 * Usage: pari_native linear|doubling
 * Standard input: L iterations, then L diagonal integers, then L bond integers.
 * Standard output: one JSON record with ascending full polynomial coefficients.
 * Parsing, pari_init and output are outside the process CPU timer. Validation,
 * least-period discovery, reduction and coefficient calculation occur afresh
 * on every measured iteration. PARI stack storage is reclaimed between solves
 * and within a solve when needed; no previous polynomial answer is reused.
 */
#define _POSIX_C_SOURCE 200809L
#include "pari.h"
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <time.h>

static void fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

static char *read_word(void)
{
    size_t used = 0, capacity = 64;
    char *word = malloc(capacity);
    int c;
    if (!word) fail("input allocation failed");
    do { c = getchar(); } while (c != EOF && isspace((unsigned char)c));
    if (c == EOF) fail("unexpected end of input");
    do {
        if (used + 1 == capacity) {
            char *grown;
            if (capacity > SIZE_MAX / 2) fail("integer token too long");
            capacity *= 2;
            grown = realloc(word, capacity);
            if (!grown) fail("input allocation failed");
            word = grown;
        }
        word[used++] = (char)c;
        c = getchar();
    } while (c != EOF && !isspace((unsigned char)c));
    word[used] = '\0';
    return word;
}

static unsigned long read_count(void)
{
    char *word = read_word(), *end;
    unsigned long result;
    size_t i;
    for (i = 0; word[i]; ++i)
        if (!isdigit((unsigned char)word[i])) fail("bad unsigned header");
    errno = 0;
    result = strtoul(word, &end, 10);
    if (errno || *end) fail("header exceeds machine range");
    free(word);
    return result;
}

static GEN read_integer(void)
{
    char *word = read_word();
    size_t offset = (word[0] == '-' || word[0] == '+') ? 1 : 0, i;
    GEN value;
    if (!word[offset]) fail("bad integer token");
    for (i = offset; word[i]; ++i)
        if (!isdigit((unsigned char)word[i])) fail("bad integer token");
    /* strtoi reads unsigned decimal digits; apply the sign explicitly. */
    value = strtoi(word + offset);
    if (word[0] == '-') value = negi(value);
    free(word);
    return value;
}

static long least_period(GEN a, GEN b, long length)
{
    long i, k, period;
    long *prefix = calloc((size_t)length, sizeof(long));
    if (!prefix) fail("period allocation failed");
    for (i = 1; i < length; ++i) {
        k = prefix[i - 1];
        while (k && (!equalii(gel(a, i + 1), gel(a, k + 1))
                     || !equalii(gel(b, i + 1), gel(b, k + 1))))
            k = prefix[k - 1];
        if (equalii(gel(a, i + 1), gel(a, k + 1))
            && equalii(gel(b, i + 1), gel(b, k + 1)))
            ++k;
        prefix[i] = k;
    }
    period = length - prefix[length - 1];
    free(prefix);
    return length % period ? length : period;
}

static GEN cell_discriminant(GEN a, GEN b, long period, GEN *bond_product)
{
    pari_sp checkpoint = avma;
    GEN p = pol_1(0), q = pol_0(0), r = pol_0(0), s = pol_1(0);
    GEN product = gen_1;
    long i;
    for (i = 1; i <= period; ++i) {
        GEN previous_bond = gel(b, i == 1 ? period : i - 1);
        GEN factor = ZX_Z_sub(pol_x(0), gel(a, i));
        GEN next_p = ZX_sub(ZX_mul(factor, p), ZX_Z_mul(r, previous_bond));
        GEN next_q = ZX_sub(ZX_mul(factor, q), ZX_Z_mul(s, previous_bond));
        r = ZX_Z_mul(p, gel(b, i));
        s = ZX_Z_mul(q, gel(b, i));
        p = next_p;
        q = next_q;
        product = mulii(product, gel(b, i));
        if (gc_needed(checkpoint, 2))
            gerepileall(checkpoint, 5, &p, &q, &r, &s, &product);
    }
    *bond_product = product;
    return ZX_add(p, s);
}

static GEN trace_linear(GEN discriminant, GEN product, unsigned long repeats)
{
    GEN square = sqri(product);
    pari_sp checkpoint = avma;
    GEN previous = scalar_ZX_shallow(gen_2, 0), current = discriminant;
    unsigned long i;
    for (i = 2; i <= repeats; ++i) {
        GEN next = ZX_sub(ZX_mul(discriminant, current), ZX_Z_mul(previous, square));
        previous = current;
        current = next;
        /* discriminant, product and square precede checkpoint and stay valid. */
        if (gc_needed(checkpoint, 2))
            gerepileall(checkpoint, 2, &previous, &current);
    }
    return current;
}

/* Return [R_n,R_(n+1)]. det(cell)=B^2, so
 * R_(2k)=R_k^2-2B^(2k),
 * R_(2k+1)=R_k R_(k+1)-D B^(2k).
 * Collect the result at each recursion boundary to reclaim temporary powers.
 */
static GEN trace_pair(GEN discriminant, GEN product, unsigned long n)
{
    pari_sp checkpoint = avma;
    GEN pair, power, twice_power, even, odd, next;
    unsigned long k;
    if (!n)
        return mkvec2(scalar_ZX_shallow(gen_2, 0), discriminant);
    k = n / 2;
    pair = trace_pair(discriminant, product, k);
    power = powiu(product, 2 * k);
    twice_power = shifti(power, 1);
    odd = ZX_sub(ZX_mul(gel(pair, 1), gel(pair, 2)),
                 ZX_Z_mul(discriminant, power));
    if (n % 2) {
        next = ZX_Z_sub(ZX_sqr(gel(pair, 2)), mulii(twice_power, sqri(product)));
        return gerepilecopy(checkpoint, mkvec2(odd, next));
    }
    even = ZX_Z_sub(ZX_sqr(gel(pair, 1)), twice_power);
    return gerepilecopy(checkpoint, mkvec2(even, odd));
}

/* The final caller needs only R_n, not R_(n+1). In particular n=1 must
 * return D immediately for aperiodic inputs rather than square a degree-L D.
 */
static GEN trace_power(GEN discriminant, GEN product, unsigned long n)
{
    pari_sp checkpoint = avma;
    GEN value, pair, answer;
    if (!n) return scalar_ZX_shallow(gen_2, 0);
    if (n == 1) return discriminant;
    if (!(n % 2)) {
        value = trace_power(discriminant, product, n / 2);
        answer = ZX_Z_sub(ZX_sqr(value), shifti(powiu(product, n), 1));
    } else {
        pair = trace_pair(discriminant, product, n / 2);
        answer = ZX_sub(ZX_mul(gel(pair, 1), gel(pair, 2)),
                        ZX_Z_mul(discriminant, powiu(product, n - 1)));
    }
    return gerepilecopy(checkpoint, answer);
}

static GEN solve(GEN a, GEN b, int doubling, long *period)
{
    long length = lg(a) - 1, i;
    unsigned long repeats;
    GEN discriminant, product, trace;
    if (length < 3 || lg(b) != lg(a)) fail("bad cycle length");
    for (i = 1; i <= length; ++i)
        if (typ(gel(a, i)) != t_INT || typ(gel(b, i)) != t_INT || !signe(gel(b, i)))
            fail("expected integer diagonals and nonzero integer bonds");
    *period = least_period(a, b, length);
    repeats = (unsigned long)(length / *period);
    discriminant = cell_discriminant(a, b, *period, &product);
    trace = doubling ? trace_power(discriminant, product, repeats)
                     : trace_linear(discriminant, product, repeats);
    return ZX_Z_sub(trace, shifti(powiu(product, repeats), 1));
}

static double cpu_time(void)
{
    struct timespec value;
    if (clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &value)) fail("process CPU clock failed");
    return (double)value.tv_sec + (double)value.tv_nsec * 1e-9;
}

int main(int argc, char **argv)
{
    unsigned long input_length, iterations, iteration;
    long length, i, period = 0;
    int doubling;
    double started, elapsed;
    GEN a, b, result = NULL;
    pari_sp input_end;
    if (argc != 2 || (strcmp(argv[1], "linear") && strcmp(argv[1], "doubling")))
        fail("usage: pari_native linear|doubling");
    doubling = !strcmp(argv[1], "doubling");
    input_length = read_count();
    iterations = read_count();
    if (input_length < 3 || input_length >= (unsigned long)LONG_MAX || !iterations)
        fail("bad cycle length or iteration count");
    length = (long)input_length;
    /* Same initial stack allocation as the GP runner; prime tables are unused. */
    pari_init(256000000, 2);
    a = cgetg(length + 1, t_VEC);
    for (i = 1; i <= length; ++i) gel(a, i) = read_integer();
    b = cgetg(length + 1, t_VEC);
    for (i = 1; i <= length; ++i) gel(b, i) = read_integer();
    input_end = avma;
    started = cpu_time();
    for (iteration = 0; iteration < iterations; ++iteration) {
        set_avma(input_end);
        result = solve(a, b, doubling, &period);
    }
    elapsed = cpu_time() - started;
    if (typ(result) != t_POL || degpol(result) != length
        || !equali1(gel(result, length + 2)))
        fail("incomplete characteristic polynomial");
    printf("{\"period\":%ld,\"iterations\":%lu,\"cpu_seconds\":%.12g,\"coefficients\":[",
           period, iterations, elapsed);
    for (i = 0; i <= length; ++i) {
        if (i) putchar(',');
        pari_printf("\"%Ps\"", gel(result, i + 2));
    }
    puts("]}");
    pari_close();
    return EXIT_SUCCESS;
}
