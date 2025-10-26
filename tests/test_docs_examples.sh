#!/usr/bin/env bash

set -euo pipefail

# qxub docs runnable examples tester
#
# Conventions:
# - Fenced code blocks must start with: ```bash runnable or ```sh runnable
# - Optional first-line comment with tags: "# tags: hpc, slow" etc.
# - Dangerous commands are rejected by a denylist.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
cd "$ROOT_DIR"

MODE="list"   # list | run
INCLUDE_TAGS=""  # comma-separated list; empty means only untagged
FILES=( "README.md" )

# Add default doc search paths
while IFS= read -r -d '' f; do FILES+=("$f"); done < <(find docs -type f -name "*.md" -not -path "docs/dev/*" -print0)

usage() {
  cat <<EOF
Usage: tests/test_docs_examples.sh [--list|--run] [--include-tag TAG] [--file PATH]

Options:
  --list                 List runnable snippets (default)
  --run                  Execute safe, runnable snippets
  --include-tag TAG      Include blocks with tag (can repeat). Example: --include-tag hpc
  --file PATH            Limit to a single Markdown file (can repeat)

Exit codes:
  0 on success, non-zero on failures.
EOF
}

DENY_PATTERNS=(
  "rm -rf"
  "sudo "
  "mkfs"
  "dd if="
  ":>"
  "> /dev/"
  "^:\\s*$"  # colons as no-op are fine; this is harmless but ignored
)

has_deny() {
  local content="$1"
  local pat
  for pat in "${DENY_PATTERNS[@]}"; do
    if grep -E -q "$pat" <<<"$content"; then
      return 0
    fi
  done
  return 1
}

BLOCK_ID=0
IN_BLOCK=0
BLOCK_LANG=""
BLOCK_CONTENT=""
BLOCK_FILE=""
BLOCK_START_LINE=0
BLOCK_TAGS=()

print_block_header() {
  echo "---"
  echo "file: $BLOCK_FILE"
  echo "line: $BLOCK_START_LINE"
  echo "id: $BLOCK_ID"
  if ((${#BLOCK_TAGS[@]})); then
    echo "tags: ${BLOCK_TAGS[*]}"
  fi
}

run_block() {
  local content="$1"
  if has_deny "$content"; then
    echo "[SKIP][denied] Potentially destructive commands detected" >&2
    return 0
  fi
  # Respect tags
  if ((${#BLOCK_TAGS[@]})); then
    local wanted=0
    IFS=',' read -r -a include <<<"$INCLUDE_TAGS"
    for t in "${BLOCK_TAGS[@]}"; do
      for it in "${include[@]}"; do
        if [[ "$t" == "$it" ]]; then wanted=1; fi
      done
    done
    if [[ "$INCLUDE_TAGS" == "" && $wanted -eq 0 ]]; then
      echo "[SKIP][tagged] Not included by default (tags: ${BLOCK_TAGS[*]})" >&2
      return 0
    fi
    if [[ "$INCLUDE_TAGS" != "" && $wanted -eq 0 ]]; then
      echo "[SKIP][tagged] Does not match include tags (tags: ${BLOCK_TAGS[*]})" >&2
      return 0
    fi
  fi
  echo "+ Running block $BLOCK_ID from $BLOCK_FILE:$BLOCK_START_LINE"
  bash -euo pipefail -c "$content"
}

ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --run) MODE="run"; shift ;;
    --list) MODE="list"; shift ;;
    --include-tag) INCLUDE_TAGS="${INCLUDE_TAGS:+$INCLUDE_TAGS,}$2"; shift 2 ;;
    --file) FILES=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 2 ;;
  esac
done

for FILE in "${FILES[@]}"; do
  [[ -f "$FILE" ]] || continue
  LINE=0
  while IFS= read -r raw || [[ -n "$raw" ]]; do
    LINE=$((LINE+1))
    if [[ $IN_BLOCK -eq 0 ]]; then
      if [[ "$raw" =~ ^\`\`\`(bash|sh)[[:space:]]+runnable[[:space:]]*$ ]]; then
        IN_BLOCK=1
        BLOCK_LANG="${BASH_REMATCH[1]}"
        BLOCK_CONTENT=""
        BLOCK_FILE="$FILE"
        BLOCK_START_LINE=$LINE
        BLOCK_TAGS=()
        BLOCK_ID=$((BLOCK_ID+1))
        FIRST_CONTENT_LINE=1
        continue
      fi
    else
      if [[ "$raw" == '```' ]]; then
        # close block
        print_block_header
        if [[ "$MODE" == "run" ]]; then
          run_block "$BLOCK_CONTENT"
        else
          echo "$BLOCK_CONTENT"
        fi
        IN_BLOCK=0
        BLOCK_LANG=""
        BLOCK_CONTENT=""
        continue
      fi
      if [[ $FIRST_CONTENT_LINE -eq 1 ]]; then
        # Optional tags line
        if [[ "$raw" =~ ^\#\ tags:\ (.*)$ ]]; then
          IFS=',' read -r -a BLOCK_TAGS <<<"${BASH_REMATCH[1]// /}"
          FIRST_CONTENT_LINE=0
          continue
        fi
        FIRST_CONTENT_LINE=0
      fi
      BLOCK_CONTENT+="$raw"$'\n'
    fi
  done < "$FILE"
done

echo "Done. Mode=$MODE"
