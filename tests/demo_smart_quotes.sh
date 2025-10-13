#!/bin/bash

# Demonstration script showing smart quote processing benefits
echo "🎯 Smart Quote Processing Demonstration"
echo

# Activate qxub environment
source venv/bin/activate

echo "1. BEFORE (Traditional ugly nested quotes):"
echo 'qxub --cmd '"'"'find /path -exec sh -c '"'"'"'"'"'"'"'"'echo "Found: $1"'"'"'"'"'"'"'"'"' _ {} \;'"'"
echo

echo "2. AFTER (Smart quotes - clean and readable):"
echo 'qxub --cmd "find /path -exec sh -c \"echo \\\"Found: \$1\\\"\" _ {} \;"'
echo

echo "3. Live demonstration:"
echo "   Traditional syntax result:"
qxub --env base --dry --cmd 'echo "Hello ${USER}"' | grep "Command to execute:"

echo "   Smart quote syntax result:"
qxub --env base --dry --cmd '"echo \"Hello ${USER}\""' | grep "Command to execute:"

echo
echo "4. Complex example with mixed variables:"
qxub --env base --dry --cmd '"echo \"User: ${USER}, Job: ${{PBS_JOBID}}, Cost: \$100\""' | grep "Command to execute:"

echo
echo "5. AWK field references preserved:"
qxub --env base --dry --cmd '"awk \"{print \\\$1, \\\$2}\" file.txt"' | grep "Command to execute:"

echo
echo "✨ Smart quotes make complex commands readable while preserving all functionality!"
