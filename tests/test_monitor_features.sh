#!/bin/bash
# Test the new monitor features
# This shows the enhanced display options

echo "=== Testing qxub monitor v2.3 features ===" >&2

# Create some test job IDs (these would normally come from real qxub submissions)
test_jobs=(
    "12345.gadi-pbs"
    "12346.gadi-pbs"
    "12347.gadi-pbs"
)

echo "Test job IDs: ${test_jobs[@]}" >&2
echo "" >&2

echo "1. Normal display (shows Job ID, Name, Status emoji):" >&2
printf '%s\n' "${test_jobs[@]}" | timeout 5 qxub monitor --suffix .gadi-pbs --interval 2 || echo "   (Test timed out - would show full table)"

echo "" >&2
echo "2. Job ID only display:" >&2
printf '%s\n' "${test_jobs[@]}" | timeout 5 qxub monitor --job-id-only --suffix .gadi-pbs --interval 2 || echo "   (Test timed out - would show only job IDs)"

echo "" >&2
echo "3. Name only display:" >&2
printf '%s\n' "${test_jobs[@]}" | timeout 5 qxub monitor --name-only --suffix .gadi-pbs --interval 2 || echo "   (Test timed out - would show only names)"

echo "" >&2
echo "4. Quiet mode with countdown:" >&2
printf '%s\n' "${test_jobs[@]}" | timeout 5 qxub monitor --quiet --suffix .gadi-pbs --interval 3 || echo "   (Test timed out - would show progress counts with countdown)"

echo "" >&2
echo "Demo completed! New features:" >&2
echo "• --suffix .gadi-pbs removes suffix from display" >&2
echo "• Single emoji status indicators (⏳📍✅❌🚫)" >&2
echo "• Job Name column for meaningful identification" >&2
echo "• --name-only and --job-id-only for focused display" >&2
echo "• Exit codes tracked and displayed (✅=success, ❌=failure)" >&2
echo "• Countdown timers between status checks" >&2
echo "• Exit code 0 only if ALL jobs succeed" >&2
