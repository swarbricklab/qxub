#!/bin/bash
#PBS -S /bin/bash
#PBS -j oe
#
# The job script for every qxub job. What it sets up comes from its
# settings: modules ($mods), a conda environment ($env) and/or a
# Singularity image ($sif, with $bind). The #PBS lines above are comments
# to any other scheduler.
#
# Each step (--pre, the command, --post) runs as it would in a plain shell:
# a failure part-way through doesn't stop the step, and its exit status is
# the last command's. --pre and --post run in this shell, so what --pre sets
# up (export, cd, source, module load) reaches the command. The command
# runs in a subshell (inside the container, for an image, with its bash or
# else its sh), so `exit` or `exec` in it ends only the command.

# qxub passes the job's settings in a file (QXUB_JOB_FILE) instead of as
# qsub -v variables; jobs from older qxub got them with -v and no file.
if [ -n "${QXUB_JOB_FILE:-}" ]; then
    if ! . "$QXUB_JOB_FILE"; then
        echo "❌ Cannot read the job settings file $QXUB_JOB_FILE" >&2
        exit 1
    fi
fi

# Older names for the image settings, still honoured if set.
sif_file="${sif_file:-${sif:-}}"

# The scheduler's ID for this job: qxub names the variable that holds it
# (QXUB_JOB_ID_VAR, e.g. PBS_JOBID); jobs from older qxub ran under PBS.
QXUB_JOB_ID="${QXUB_JOB_ID_VAR:+${!QXUB_JOB_ID_VAR}}"
QXUB_JOB_ID="${QXUB_JOB_ID:-${PBS_JOBID:-}}"

# Logging
echo "📁 Execution directory: $cwd"
[ -n "${mods:-}" ] && echo "📦 Modules to load: $mods"
[ -n "${env:-}" ] && echo "🐍 Conda environment: $env"
[ -n "$sif_file" ] && echo "📦 Singularity container: $sif_file"
[ -n "${bind:-}" ] && echo "🔗 Bind mounts: $bind"
[ -z "${mods:-}${env:-}$sif_file" ] && echo "🔧 Default execution (no environment activation)"
echo "🔖 qxub version: ${qxub_version:-unknown}"
echo "📦 qxub path: ${qxub_path:-unknown}"
echo "🖥️  Submitted from: ${qxub_submit_host:-unknown}"
echo "🏠 Hostname: $(hostname)"

# Decode base64 encoded commands
if [ -n "$cmd_b64" ]; then
    cmd=$(echo "$cmd_b64" | base64 -d)
fi
if [ -n "$pre_cmd_b64" ]; then
    pre_cmd=$(echo "$pre_cmd_b64" | base64 -d)
fi
if [ -n "$post_cmd_b64" ]; then
    post_cmd=$(echo "$post_cmd_b64" | base64 -d)
fi
# Decode user-defined environment variables
if [ -n "$user_vars_b64" ]; then
    user_vars=$(echo "$user_vars_b64" | base64 -d)
fi

echo "💻 Command: $cmd"
if [ -n "$pre_cmd" ]; then
    echo "⚡ Pre-command: $pre_cmd"
fi
if [ -n "$post_cmd" ]; then
    echo "🏁 Post-command: $post_cmd"
fi
echo "📄 STDOUT file: $out"
echo "📄 STDERR file: $err"
echo "🔇 Quiet mode: $quiet"
if [ "$quiet" = "true" ]; then
    echo "ℹ️  Note: In quiet mode, output will go to the job log instead of separate files"
fi
echo "👤 Current user: $(whoami)"
echo "📍 Current directory: $(pwd)"

# Strict mode for this script's own steps (not for the user's commands)
set -e

# Create status directory for file-based monitoring
# qxub passes the status directory; jobs from older qxub derive it from $out.
STATUS_DIR="${QXUB_STATUS_DIR:-${out%/*}/status}"
if ! mkdir -p "$STATUS_DIR" 2>/dev/null; then
    # e.g. another project's scratch, not mounted for this job: carry on
    # next to $out rather than failing the job.
    echo "⚠️  Cannot create the status directory $STATUS_DIR; using ${out%/*}/status instead" >&2
    STATUS_DIR="${out%/*}/status"
    mkdir -p "$STATUS_DIR"
fi
# Status files are named after the qxub run ID, not the scheduler's job ID.
# native_<job ID> lets readers that only know the job ID find the run ID.
STATUS_KEY="${QXUB_RUN_ID:-${QXUB_JOB_ID:-unknown}}"
if [ -n "${QXUB_RUN_ID:-}" ] && [ -n "$QXUB_JOB_ID" ]; then
    echo "$QXUB_RUN_ID" > "$STATUS_DIR/native_${QXUB_JOB_ID}" 2>/dev/null || true
fi

# Signal job started with timestamp and job ID
echo "$(date -Iseconds)" > "$STATUS_DIR/job_started_${STATUS_KEY}"
echo "${QXUB_JOB_ID:-unknown}" >> "$STATUS_DIR/job_started_${STATUS_KEY}"
echo "🚀 Job started at $(date -Iseconds)"

# Run record: one file per run ($QXUB_RUN_FILE), replaced by rename on each
# update. Jobs never open the job database: on Lustre with localflock,
# SQLite writes from many nodes corrupt it (docs/dev/sqlite-robustness.md).
qxub_run_update() {
    [ -n "${QXUB_RUN_FILE:-}" ] || return 0
    local tmp="${QXUB_RUN_FILE%/*}/.${QXUB_RUN_FILE##*/}.$$.tmp"
    mkdir -p "${QXUB_RUN_FILE%/*}" 2>/dev/null || true
    {
        if [ -f "$QXUB_RUN_FILE" ]; then cat "$QXUB_RUN_FILE"; fi
        printf '%s\n' "$@"
    } > "$tmp" 2>/dev/null && mv -f "$tmp" "$QXUB_RUN_FILE" 2>/dev/null || rm -f "$tmp" 2>/dev/null || true
    return 0
}

# Write a status file: a value and the time.
qxub_status() {
    { echo "$2"; date -Iseconds; } > "$STATUS_DIR/$1_${STATUS_KEY}" 2>/dev/null || true
}

# Record that the job is running
qxub_run_update state=running "started_at=$(date -Iseconds)" "host=$(hostname)" "native_job_id=$QXUB_JOB_ID"
# Every exit path (success, the ERR trap's exit, an explicit exit) records
# the final state. PBS stops a job (walltime, qdel, memory) with SIGTERM;
# without a TERM trap bash would run the EXIT trap with the last completed
# command's status, often 0, so a killed job looked successful. The signal
# traps exit with 128+signal instead. A job that gets SIGKILL records
# nothing and stays "running"; the joblog settles it. A --pre or --post
# that sets its own EXIT trap replaces this one, because it runs here.
#
# --pre and --post run in this shell, so an `exit` in one ends the job
# before any status file says so. Record the exit code then, so `qxub exec`
# doesn't wait for the joblog. Signal exits are left to the joblog, whose
# exit status says why PBS stopped the job.
qxub_on_exit() {
    local code=$?
    if [ -z "${QXUB_SIGNAL:-}" ] && [ -n "${STATUS_DIR:-}" ] && [ -n "${STATUS_KEY:-}" ] \
        && [ ! -f "$STATUS_DIR/final_exit_code_${STATUS_KEY}" ]; then
        qxub_status final_exit_code "$code"
        qxub_record_step "$code"
    fi
    qxub_run_update state=finished "exit_code=$code" "signal=${QXUB_SIGNAL:-}" "finished_at=$(date -Iseconds)"
}
trap qxub_on_exit EXIT
trap 'QXUB_SIGNAL=TERM; exit 143' TERM
trap 'QXUB_SIGNAL=INT; exit 130' INT
trap 'QXUB_SIGNAL=HUP; exit 129' HUP

# A setup step failed: say what failed and how to fix it, on stderr (the
# joblog) and in $err, which `qxub exec` and CI steps stream, then stop.
# Without this a bad execdir was an opaque "exit code 1" (#109).
qxub_setup_fail() {
    printf '%s\n' "$@" >&2
    if [ "$quiet" != "true" ] && [ -n "${err:-}" ]; then
        { mkdir -p "${err%/*}" && printf '%s\n' "$@" >> "$err"; } 2>/dev/null || true
    fi
    qxub_status final_exit_code 1
    exit 1
}

# Record the running step's exit code if the step didn't (it called exit),
# so the summary can say which step failed.
qxub_record_step() {
    [ -n "${QXUB_STEP:-}" ] && [ -n "${STATUS_DIR:-}" ] || return 0
    [ -f "$STATUS_DIR/${QXUB_STEP}_exit_code_${STATUS_KEY}" ] || qxub_status "${QXUB_STEP}_exit_code" "$1"
}

qxub_notify() {
    if [ "$qxub_notify" = "true" ]; then
        echo "📬 Queueing notification job..."
        qxub exec --internet --terse --default -- qxub notify send \
            --job-id "${QXUB_JOB_ID:-unknown}" \
            --exit-code "$1" \
            --job-name "${QXUB_JOB_NAME:-${PBS_JOBNAME:-unknown}}" \
            --output-dir "${qxub_notify_outdir:-$PWD}" \
            2>/dev/null || echo "⚠️  Failed to queue notification (non-fatal)"
    fi
}

# This script's own commands failing (not the user's): record and stop.
error_handler() {
    local exit_code=$?
    qxub_record_step "$exit_code"
    echo "❌ ERROR: Command failed with exit code $exit_code"
    echo "❌ Failed command: $BASH_COMMAND"
    qxub_status final_exit_code "$exit_code"
    qxub_notify "$exit_code"
    exit $exit_code
}
trap error_handler ERR

# A step ended: record its exit code, and stop the job if it failed.
qxub_step_done() {
    local step=$1 code=$2 label=$3
    qxub_status "${step}_exit_code" "$code"
    if [ "$code" -ne 0 ]; then
        echo "❌ $label failed with exit code $code"
        qxub_status final_exit_code "$code"
        qxub_notify "$code"
        exit "$code"
    fi
    echo "✅ $label completed"
}

# Check required variables
if [ -z "$cmd" ]; then
    qxub_setup_fail "❌ ERROR: no command was given"
fi

if [ -z "$cwd" ]; then
    qxub_setup_fail "❌ ERROR: no execution directory was given (cwd is empty)"
fi

# ==================== MODULES (loaded first) ====================
if [ -n "${mods:-}" ]; then
    echo "📦 Loading modules: $mods"
    for module in $mods; do
        echo "🔄 Loading module: $module"
        module load "$module" || qxub_setup_fail \
            "❌ ERROR: could not load module: $module" \
            "   (check the name with: module avail $module)"
    done
    echo "✅ All modules loaded"
fi

# ==================== CONDA (activated after modules) ====================
if [ -n "${env:-}" ]; then
    echo "🔄 Activating conda environment: $env"

    # Check if conda is available and properly initialized
    if ! command -v conda > /dev/null 2>&1; then
        qxub_setup_fail \
            "❌ ERROR: conda command not found" \
            "💡 Possible solutions:" \
            "   1. Run 'conda init bash' to initialize conda in your shell" \
            "   2. Add conda to your PATH in ~/.bashrc" \
            "   3. Ensure conda storage volumes are mounted (e.g., -l storage=gdata/conda)" \
            "   4. Check that your conda installation path is accessible from compute nodes"
    fi

    # Initialize conda shell integration
    echo "🔧 Initializing conda shell integration..."
    if ! eval "$(conda shell.bash hook)" 2>&1; then
        qxub_setup_fail \
            "❌ ERROR: Failed to initialize conda shell integration" \
            "💡 Possible solutions:" \
            "   1. Run 'conda init bash' in your login shell to set up conda properly" \
            "   2. Check that conda is properly installed and accessible" \
            "   3. Verify that your ~/.bashrc contains conda initialization code" \
            "   4. Ensure conda storage volumes are mounted in your PBS job"
    fi

    if ! conda activate "$env" 2>&1; then
        qxub_setup_fail \
            "❌ ERROR: Failed to activate conda environment: $env" \
            "💡 Possible solutions:" \
            "   1. Check that the environment exists: conda env list" \
            "   2. Create the environment: conda create -n $env" \
            "   3. Verify the environment name is spelled correctly"
    fi

    # CONDA_DEFAULT_ENV may be just the name or a full path
    ACTIVATED_ENV_NAME=$(basename "$CONDA_DEFAULT_ENV")
    if [ "$ACTIVATED_ENV_NAME" != "$(basename "$env")" ]; then
        qxub_setup_fail \
            "❌ ERROR: Conda environment activation failed - expected '$env', got '$CONDA_DEFAULT_ENV'"
    fi
    echo "✅ Successfully activated conda environment: $env"
fi

# ==================== SINGULARITY ====================
if [ -n "$sif_file" ]; then
    echo "📦 Loading singularity module"
    module load singularity || qxub_setup_fail \
        "❌ ERROR: could not load the singularity module" \
        "   (check it with: module avail singularity)"

    # Bind mounts are one argument, however many paths (/a,/b) they name.
    sing_args=()
    if [ -n "${bind:-}" ]; then
        sing_args+=(--bind "$bind")
    fi
    if [ -n "${sing_options:-}" ]; then
        # Older free-form options, split as before
        read -ra extra_sing_options <<< "$sing_options"
        sing_args+=("${extra_sing_options[@]}")
    fi
fi

# Create execution directory if requested
if [ "$create_execdir" = "true" ]; then
    echo "📁 Creating execution directory: $cwd"
    mkdir -p "$cwd" || qxub_setup_fail "❌ ERROR: could not create the execution directory: $cwd"
fi

# Change to execution directory
cd "$cwd" || qxub_setup_fail \
    "❌ ERROR: execution directory does not exist or can't be entered: $cwd" \
    "   (pass --create-execdir, or make sure it is created upstream)"

# A local image must exist (relative paths are from the execution directory);
# URIs such as docker:// or library:// are left to singularity.
if [ -n "$sif_file" ]; then
    case "$sif_file" in
        *://*) ;;
        *) [ -e "$sif_file" ] || qxub_setup_fail \
               "❌ ERROR: Singularity image not found: $sif_file" \
               "   (check the path; relative paths are from $cwd)" ;;
    esac
fi

# Export user-defined environment variables
if [ -n "$user_vars" ]; then
    echo "🔧 Setting user-defined environment variables:"
    IFS=';' read -ra VARS <<< "$user_vars"
    for var_def in "${VARS[@]}"; do
        if [ -n "$var_def" ]; then
            echo "   export $var_def"
            export "$var_def"
        fi
    done
fi

# Outside quiet mode the main command writes $out and $err. Pre- and
# post-commands, and this script's own messages, go to the joblog (#150).
if [ "$quiet" != "true" ]; then
    mkdir -p "$(dirname "$out")" "$(dirname "$err")"
    touch "$out" "$err" || qxub_setup_fail \
        "❌ ERROR: cannot write the output files: $out and $err"
fi

# The user's steps run without this script's strict mode and error trap,
# as they would in a plain shell; qxub_step_done then checks their status.
# (The trap is changed here, not in a function: bash restores a function's
# ERR trap when it returns.)
qxub_user_step() {
    QXUB_STEP=$1
    qxub_status "$1_started" ""
}

# The command, in a subshell: inside the container if there is one, with
# the image's bash, or its sh if it has no bash (minimal images).
qxub_main() {
    if [ -n "$sif_file" ]; then
        singularity exec "${sing_args[@]}" "$sif_file" /bin/sh -c \
            'if command -v bash >/dev/null 2>&1; then exec bash -c "$1"; fi
             exec sh -c "$1"' qxub "$cmd"
    else
        eval "$cmd"
    fi
}

# Execute pre-command if provided
if [ -n "$pre_cmd" ]; then
    echo "⚡ Executing pre-command..."
    qxub_user_step pre
    set +e; trap - ERR
    eval "$pre_cmd"
    PRE_EXIT_CODE=$?
    set -e; trap error_handler ERR
    qxub_step_done pre "$PRE_EXIT_CODE" "Pre-command"
fi

# Execute main command
echo "💻 Executing main command..."
qxub_user_step main
set +e; trap - ERR
if [ "$quiet" = "true" ]; then
    echo "================== COMMAND OUTPUT START =================="
    ( qxub_main )
    MAIN_EXIT_CODE=$?
    echo "=================== COMMAND OUTPUT END ==================="
else
    # Redirect only this command: an `exec` redirection would send the rest
    # of the script to $out and truncate it again at the next step.
    ( qxub_main ) > "$out" 2> "$err"
    MAIN_EXIT_CODE=$?
fi
set -e; trap error_handler ERR
qxub_step_done main "$MAIN_EXIT_CODE" "Main command"

# Execute post-command if provided
if [ -n "$post_cmd" ]; then
    echo "🏁 Executing post-command..."
    qxub_user_step post
    set +e; trap - ERR
    eval "$post_cmd"
    POST_EXIT_CODE=$?
    set -e; trap error_handler ERR
    qxub_step_done post "$POST_EXIT_CODE" "Post-command"
fi

# Every step succeeded
qxub_status final_exit_code 0
qxub_notify 0

echo "🎉 Job completed successfully"
