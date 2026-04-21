#pragma once

// ============================================================================
// ABC/ABCD/ABC2 File Format Parser
//
// Supports srsieve2, gcwsieve, mfsieve output formats for batch processing.
// Provides a CandidateSource abstraction for lazy or materialized iteration
// over candidate expressions parsed from these formats.
//
// ABC   — sieved output with template:  ABC $a*2^$b+1
// ABCD  — delta-encoded:               ABCD $a*2^$b+1 [3 1000]
// ABC2  — iterative range definitions:  ABC2 $a*2^$b+1
// ============================================================================

#include <string>
#include <vector>
#include <memory>
#include <cstdint>

class Logging;

// ============================================================================
// Constants
// ============================================================================

// Maximum sieve range for ABC2 "primes from X to Y" (100 million)
static const int64_t MAX_SIEVE_RANGE = 100000000LL;

// Maximum total ABC2 Cartesian product candidates
static const size_t MAX_ABC2_CANDIDATES = 100000000ULL;

// Maximum values in an ABC2 "in { ... }" explicit list
static const size_t MAX_ABC2_IN_LIST_VALUES = 1000;

// ============================================================================
// File format detection
// ============================================================================

enum FileFormat { FORMAT_UNKNOWN, FORMAT_ABC, FORMAT_ABCD, FORMAT_ABC2 };

FileFormat detect_format(const std::string& first_line);

// ============================================================================
// Candidate: a single number expression with optional k-value for tracking
// ============================================================================

struct Candidate
{
    std::string expression;
    std::string k_value;   // empty if not applicable (non-ABC formats)
};

// ============================================================================
// CandidateSource: abstract interface for batch candidate iteration
//
// Supports both materialized (vector-backed) and lazy (on-demand) sources.
// All implementations are thread-safe for const get() calls.
// ============================================================================

class CandidateSource
{
public:
    virtual ~CandidateSource() = default;

    // Total number of candidates
    virtual size_t size() const = 0;

    // Get candidate at given index (0-based). Returns false on invalid index.
    virtual bool get(size_t index, Candidate& out) const = 0;

    // Whether this source uses ABC-family format (has k-values)
    virtual bool is_abc() const = 0;
};

// ============================================================================
// Factory: parse a batch file and return appropriate CandidateSource
//
// Reads the file, auto-detects format, and returns a CandidateSource.
// For ABC/ABCD: materializes all candidates into memory (typical for sieve output).
// For ABC2: uses lazy generation to avoid materializing the Cartesian product.
// For raw (unknown format): wraps lines as-is with no k-values.
//
// Returns nullptr on failure (empty file, parse error, etc.).
// ============================================================================

std::unique_ptr<CandidateSource> parse_batch_file(
    const std::string& filename,
    Logging& logging);
