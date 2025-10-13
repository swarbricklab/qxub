# Command Handling Edge Cases - Testing Matrix

## Overview

This document provides a comprehensive taxonomy of edge cases for testing qxub's command handling across all execution modes (conda, modules, singularity, default). All modes use the same base64 encoding strategy, so edge cases apply universally.

## Current Architecture Summary

**Encoding Flow:**
1. Command parsed by Click: `command: Tuple[str, ...]`
2. Join with spaces: `cmd_str = " ".join(command)`
3. Base64 encode: `cmd_b64 = base64.b64encode(cmd_str.encode("utf-8")).decode("ascii")`
4. Submit to PBS: `qsub -v cmd_b64="{encoded}" template.pbs`
5. Decode in job: `cmd=$(echo "$cmd_b64" | base64 -d)`
6. Execute: `bash -c "$cmd"`

## Edge Case Taxonomy

### 1. Quoting and Escaping Edge Cases

#### 1.1 Nested Quoting Scenarios
**Current Coverage:** ✅ Basic nested quotes tested
**Gaps:** Complex nesting, mixed quote types

```bash
# Test cases needed:
qxub --env base -- echo "outer 'single' quotes"
qxub --env base -- echo 'outer "double" quotes'
qxub --env base -- python -c "print('single inside double: \"nested double\"')"
qxub --env base -- bash -c 'echo "level1 '\''level2'\'' back to level1"'

# Pathological case:
qxub --env base -- python -c "exec(\"print('deeply \\\"nested\\\" quotes')\")"
```

#### 1.2 Backslash Escaping Chain
**Current Coverage:** ✅ Some escaping tested
**Gaps:** Multi-level escaping through the toolchain

```bash
# Each level of toolchain may process backslashes:
# User → Shell → Click → Python join → Base64 → PBS → Bash → Final command

# Test cases:
qxub --env base -- echo "\\n"                    # Literal \n
qxub --env base -- echo "\\\\n"                  # Literal \\n
qxub --env base -- python -c "print(\"\\n\")"   # Python newline
qxub --env base -- python -c "print(\"\\\\n\")" # Python literal \n
qxub --env base -- sed 's/\t/\\t/g'              # Regex escaping
```

#### 1.3 Quote Character Edge Cases
**Current Coverage:** ❌ Not tested
**Gaps:** Special quote characters, Unicode quotes

```bash
# Test cases:
qxub --env base -- echo "backtick: \`date\`"
qxub --env base -- echo 'curly quotes: "text"'  # Unicode quotes
qxub --env base -- echo "apostrophe: don't"
qxub --env base -- python -c "print(\"triple quotes: '''text'''\")"
```

### 2. Shell Metacharacter Edge Cases

#### 2.1 Pipeline and Redirection Chains
**Current Coverage:** ✅ Basic pipes tested
**Gaps:** Complex pipelines, advanced redirection

```bash
# Test cases needed:
qxub --env base -- bash -c "echo test | tee /tmp/out | wc -l | tee -a /tmp/log"
qxub --env base -- bash -c "cmd 2>&1 | grep error >&2"
qxub --env base -- bash -c "exec 3>&1 4>&2; echo output >&3; echo error >&4"

# Here documents:
qxub --env base -- bash -c "cat <<EOF > output.txt
Line 1
Line 2 with \$VAR
EOF"

# Process substitution:
qxub --env base -- bash -c "diff <(sort file1) <(sort file2)"
```

#### 2.2 Command Substitution Edge Cases
**Current Coverage:** ❌ Not tested systematically
**Gaps:** Nested substitution, mixed syntax

```bash
# Test cases:
qxub --env base -- bash -c "echo \"Today is \$(date)\""
qxub --env base -- bash -c "echo \"Today is \`date\`\""
qxub --env base -- bash -c "echo \"Files: \$(ls | wc -l)\""

# Nested substitution:
qxub --env base -- bash -c "echo \"\$(basename \$(dirname \$(pwd)))\""
qxub --env base -- bash -c "echo \"\$(echo \"\$(date)\" | cut -d' ' -f1)\""
```

#### 2.3 Globbing and Pattern Matching
**Current Coverage:** ❌ Not tested
**Gaps:** Glob patterns, bracket expressions

```bash
# Test cases:
qxub --env base -- bash -c "ls *.py"
qxub --env base -- bash -c "echo file{1,2,3}.txt"
qxub --env base -- bash -c "find . -name \"*.[ch]\" -exec grep \"main\" {} +"
qxub --env base -- bash -c "for f in *.txt; do echo \$f; done"

# Bracket patterns:
qxub --env base -- bash -c "ls file[0-9].txt"
qxub --env base -- bash -c "echo file[!0-9].txt"
```

### 3. Variable Expansion Edge Cases

#### 3.1 Submission vs Execution Time Expansion
**Current Coverage:** ✅ Basic variables tested
**Gaps:** Systematic analysis of expansion timing

```bash
# Variables that should expand at SUBMISSION time:
qxub --env base -- echo "Submitted by: $USER"     # Current user
qxub --env base -- echo "Submit host: $HOSTNAME"  # Current host

# Variables that should expand at EXECUTION time:
qxub --env base -- echo "Job user: \$USER"        # Job execution user
qxub --env base -- echo "Compute node: \$HOSTNAME" # Compute node
qxub --env base -- echo "PBS job ID: \$PBS_JOBID" # Only available in job

# Mixed scenarios:
qxub --env base -- echo "Job \$PBS_JOBID submitted by $USER on \$HOSTNAME"
```

#### 3.2 Special Variable Cases
**Current Coverage:** ❌ Not tested
**Gaps:** Special variables, parameter expansion

```bash
# Test cases:
qxub --env base -- bash -c "echo \"Args: \$@\""
qxub --env base -- bash -c "echo \"Script name: \$0\""
qxub --env base -- bash -c "echo \"Exit code: \$?\""
qxub --env base -- bash -c "echo \"Process ID: \$\$\""

# Parameter expansion:
qxub --env base -- bash -c "VAR=test; echo \"\${VAR:-default}\""
qxub --env base -- bash -c "FILE=/path/to/file.txt; echo \"\${FILE##*/}\""
qxub --env base -- bash -c "echo \"\${PWD%/*}\""
```

#### 3.3 Array and Associative Array Variables
**Current Coverage:** ❌ Not tested
**Gaps:** Complex variable structures

```bash
# Test cases:
qxub --env base -- bash -c "ARRAY=(a b c); echo \"\${ARRAY[@]}\""
qxub --env base -- bash -c "declare -A MAP; MAP[key]=value; echo \"\${MAP[key]}\""
qxub --env base -- bash -c "IFS=':' read -ra PATHS <<< \"\$PATH\"; echo \"\${PATHS[0]}\""
```

### 4. Encoding and Character Set Edge Cases

#### 4.1 Unicode and Multi-byte Characters
**Current Coverage:** ❌ Not tested
**Gaps:** UTF-8 handling through base64 chain

```bash
# Test cases:
qxub --env base -- echo "Unicode: 你好世界"
qxub --env base -- echo "Emoji: 🚀🔬📊"
qxub --env base -- python -c "print('Greek: αβγδε')"
qxub --env base -- echo "Accents: café résumé naïve"

# Unicode in filenames:
qxub --env base -- bash -c "touch 'файл.txt' && ls *.txt"
```

#### 4.2 Control Characters and Special Sequences
**Current Coverage:** ❌ Not tested
**Gaps:** Non-printable characters

```bash
# Test cases:
qxub --env base -- printf "Tab:\\t Newline:\\n Null:\\0"
qxub --env base -- echo -e "\\x1b[31mRed text\\x1b[0m"  # ANSI colors
qxub --env base -- echo "Bell: \a"
qxub --env base -- python -c "print('Vertical tab: \\v')"
```

#### 4.3 Binary Data and Large Payloads
**Current Coverage:** ❌ Not tested
**Gaps:** Base64 size limits, binary handling

```bash
# Test large command payload:
LARGE_CMD="python -c \"import sys; data='x'*10000; print(len(data))\""
qxub --env base -- "$LARGE_CMD"

# Test command with embedded binary-like data:
qxub --env base -- python -c "print(bytes([0,1,2,255]).hex())"
```

### 5. Interactive and Signal Handling Edge Cases

#### 5.1 Interactive Command Simulation
**Current Coverage:** ❌ Not tested
**Gaps:** Commands expecting stdin

```bash
# Test cases:
qxub --env base -- bash -c "echo 'test input' | python -c \"import sys; print(sys.stdin.read())\""
qxub --env base -- bash -c "yes | head -5"
qxub --env base -- bash -c "echo -e 'line1\\nline2' | while read line; do echo \"Read: \$line\"; done"
```

#### 5.2 Long-Running and Background Commands
**Current Coverage:** ❌ Not tested
**Gaps:** Commands with different execution patterns

```bash
# Test cases:
qxub --env base -- bash -c "sleep 5 && echo 'done'"
qxub --env base -- bash -c "python -c \"import time; [print(i) or time.sleep(1) for i in range(3)]\""
qxub --env base -- bash -c "(sleep 2; echo 'background') &; echo 'foreground'; wait"
```

### 6. Environment and Context Edge Cases

#### 6.1 Environment Variable Inheritance
**Current Coverage:** ✅ Basic env vars tested
**Gaps:** Complex environment scenarios

```bash
# Test cases:
qxub --env base -- bash -c "env | grep PBS"
qxub --env base -- bash -c "echo \$CONDA_DEFAULT_ENV"
qxub --env base -- bash -c "echo \$PATH | tr ':' '\\n' | head -3"
qxub --env base -- python -c "import os; print(os.environ.get('PYTHONPATH', 'Not set'))"
```

#### 6.2 Working Directory and Path Edge Cases
**Current Coverage:** ✅ Basic execdir tested
**Gaps:** Complex path scenarios

```bash
# Test cases:
qxub --env base -- bash -c "pwd && ls -la"
qxub --env base -- bash -c "cd /tmp && pwd && cd - && pwd"
qxub --env base -- python -c "import os; print(f'CWD: {os.getcwd()}')"

# Relative path handling:
qxub --execdir ./subdir --env base -- bash -c "pwd && ls ../"
```

### 7. Cross-Execution Mode Edge Cases

#### 7.1 Mode-Specific Environment Conflicts
**Current Coverage:** ❌ Not systematically tested
**Gaps:** Environment interactions

```bash
# Conda vs system Python:
qxub --env myenv -- python -c "import sys; print(sys.executable)"
qxub --default -- python -c "import sys; print(sys.executable)"

# Module loading conflicts:
qxub --mod python/3.9 --mod python/3.11 -- python --version  # Should fail
qxub --mod python/3.9 -- python -c "import sys; print(sys.version)"

# Singularity path isolation:
qxub --sif container.sif -- bash -c "echo \$PATH"
qxub --default -- bash -c "echo \$PATH"
```

#### 7.2 Pre/Post Command Interactions
**Current Coverage:** ✅ Basic pre/post tested
**Gaps:** Complex interactions

```bash
# Test cases:
qxub --env base --pre "export MYVAR=pre_value" --post "echo Final: \$MYVAR" -- echo "Main: \$MYVAR"
qxub --env base --pre "mkdir -p /tmp/testdir" --post "rm -rf /tmp/testdir" -- bash -c "cd /tmp/testdir && pwd"

# Error handling in pre/post:
qxub --env base --pre "false" -- echo "Should this run?"
qxub --env base --post "exit 42" -- echo "Main command"
```

## Testing Methodology

### Test Categories

1. **Dry Run Tests** - Quick validation without job submission
2. **Integration Tests** - Full job submission and monitoring
3. **Comparison Tests** - Cross-execution mode consistency
4. **Failure Tests** - Expected failure scenarios

### Expected Behaviors

**Transparent Pass-Through:** Commands should execute identically to direct bash execution
**Consistent Encoding:** All execution modes should handle commands identically
**Variable Timing:** Clear distinction between submission-time vs execution-time expansion
**Error Propagation:** Failed commands should properly propagate exit codes

### Test Implementation Strategy

1. **Extend existing test files** with new edge case categories
2. **Create focused test scripts** for specific areas (unicode, variables, etc.)
3. **Add comparison framework** to verify consistency across execution modes
4. **Document expected behaviors** for each edge case category

## Priority Recommendations

**High Priority:**
- Variable expansion timing (submission vs execution)
- Unicode and encoding edge cases
- Complex quoting scenarios
- Large command payloads

**Medium Priority:**
- Advanced shell metacharacters
- Cross-mode consistency verification
- Interactive command simulation

**Low Priority:**
- Pathological edge cases
- Performance stress testing
- Binary data handling
