# MT_Forward_Modeling C++ Unified Coding Standard

## 1. Purpose and Scope

This standard applies to all C++ source code in this repository, including:

- `src/*.h`
- `src/*.cpp`
- future unit tests and tools

Goals:

- Easy to read
- Easy to understand
- Easy to encapsulate
- Easy to extend
- Easy to maintain

## 2. Baseline Requirements

1. Language standard: **C++17**
2. Formatting tool: **clang-format**
3. Static analysis tool: **clang-tidy**
4. CI/PR gate: formatting and lint checks must pass before merge

## 3. Naming Conventions

### 3.1 Type Names

- Use **PascalCase**
- Applies to: class, struct, enum, type alias

Examples:

- `FiniteElementSolver`
- `SolverParameters`
- `PointLocation`

### 3.2 Variables and Functions

- Use **snake_case**
- Applies to: local variables, member variables, function names, parameters

Examples:

- `polynomial_order`
- `source_identifier`
- `assemble_source_term()`

### 3.3 Constants

- Use **kCamelCase**
- Applies to: `constexpr`, `const` globals, file-scope constants, static class constants

Examples:

- `constexpr double kPi = 3.141592653589793;`
- `static constexpr int kMpiTagStations = 101;`

### 3.4 Semantic Accuracy

- Names must match real meaning.
- Avoid vague names like `tmp`, `data1`, `val`.
- Prefer domain-explicit names:
  - `h_iter` -> `max_refinement_iterations`
  - `beta` -> `refinement_threshold_ratio`
  - `p_order` -> `polynomial_order`

## 4. Header File Rules

1. **Forbidden in headers**: `using namespace ...;`
2. In headers, always use fully qualified names (`std::`, `mfem::`).
3. Keep includes minimal; prefer forward declarations when practical.
4. Every header must be self-contained and compilable when included alone.

## 5. Include Order

Use the following order in each `.cpp`:

1. Corresponding header (if any)
2. C++ standard library headers
3. Third-party headers
4. Project headers

Within each group, sort alphabetically.

## 6. Formatting Rules (clang-format)

1. No manual alignment by spaces.
2. Braces use a single consistent style across the project.
3. Max line length: **100**.
4. One declaration per line when it improves readability.
5. Keep functions visually compact, but do not sacrifice clarity.

All formatting decisions are delegated to `clang-format`; do not hand-format against it.

## 7. clang-tidy Rules

Enable checks that enforce:

1. Modernize (`modernize-*`)
2. Readability (`readability-*`)
3. Bug-prone patterns (`bugprone-*`)
4. Performance (`performance-*`)
5. Core C++ guidelines (`cppcoreguidelines-*`, with project-specific exceptions)

At minimum, disallow:

- raw `new`/`delete` in ordinary business logic
- C-style casts
- overly complex functions without decomposition

## 8. Class and Function Design

1. Member fields default to `private`.
2. Public API should be minimal and intention-revealing.
3. Prefer RAII and smart pointers (`std::unique_ptr`) over manual lifetime management.
4. Single-responsibility functions:
   - If a function grows too large, split by stage/intent.
5. Use `const` correctness aggressively:
   - `const` parameters/reference where mutation is not intended
   - `const` member functions when state does not change

## 9. Error Handling and Logging

1. Runtime errors: use consistent exception strategy.
2. `assert` is only for debug-time invariant checks, not user/input error handling.
3. Use unified logging format:
   - `[rank][module][iteration] message`
4. Do not scatter ad-hoc `std::cout` prints for production flow control.

## 10. Comments and Documentation

1. Comments should explain **why**, not restate **what** obvious code does.
2. Public APIs should have concise doc comments for inputs/outputs/units.
3. Keep scientific formulas and references, but attach them to relevant logic blocks.
4. Remove stale comments when code changes.

## 11. Magic Numbers and Configuration

1. No unexplained magic numbers.
2. Promote repeated literals to named constants (`kCamelCase`).
3. MPI tags, tolerances, and thresholds must be named constants.

## 12. Suggested Tooling Files

Repository should include:

- `.clang-format`
- `.clang-tidy`
- optional `CMakePresets.json` or build scripts that run format/lint targets

## 13. Local Developer Workflow

Before commit:

1. Build with C++17
2. Run `clang-format` on changed files
3. Run `clang-tidy` on changed files
4. Run tests / smoke run

Pull request checklist:

1. Naming follows PascalCase/snake_case/kCamelCase
2. No `using namespace` in headers
3. Format/lint checks pass
4. No newly introduced raw pointer ownership

## 14. Definition of Done (Style)

A change is style-complete only when:

1. It conforms to this document
2. It passes automated formatting/lint checks
3. Names accurately represent real semantics
4. It does not increase coupling or reduce encapsulation
