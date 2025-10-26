#!/usr/bin/env python3
"""
Import testing utility for qxub package migration.

This script tests whether Python packages/modules can be imported successfully.
Useful for validating package migrations and identifying import issues.

Usage:
    python tests/check_imports.py qxub.execution
    python tests/check_imports.py qxub.config qxub.history qxub.resources
    python tests/check_imports.py --all  # Test all qxub packages
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List


def test_import(package: str) -> bool:
    """
    Test if a package can be imported successfully.

    Args:
        package: Package name to test (e.g., 'qxub.execution')

    Returns:
        True if import successful, False otherwise
    """
    try:
        # Use subprocess to isolate import test
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                f"import {package}; print('✅ {package} imports successfully')",
            ],
            capture_output=True,
            text=True,
            timeout=10,
        )

        if result.returncode == 0:
            print(result.stdout.strip())
            return True
        else:
            print(f"❌ {package} has import errors:")
            if result.stderr:
                # Print first few lines of error for brevity
                error_lines = result.stderr.strip().split("\n")
                for line in error_lines[:5]:  # Show first 5 lines
                    print(f"   {line}")
                if len(error_lines) > 5:
                    print(f"   ... ({len(error_lines) - 5} more lines)")
            return False

    except subprocess.TimeoutExpired:
        print(f"❌ {package} import timed out")
        return False
    except Exception as e:
        print(f"❌ {package} failed with exception: {e}")
        return False


def get_qxub_packages() -> List[str]:
    """Get list of all qxub sub-packages."""
    qxub_dir = Path(__file__).parent.parent / "qxub"
    packages = ["qxub"]

    # Find all sub-packages (directories with __init__.py)
    for item in qxub_dir.iterdir():
        if item.is_dir() and (item / "__init__.py").exists():
            packages.append(f"qxub.{item.name}")

    return sorted(packages)


def main():
    parser = argparse.ArgumentParser(
        description="Test Python package imports",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python tests/check_imports.py qxub.execution
  python tests/check_imports.py qxub.config qxub.history
  python tests/check_imports.py --all
  python tests/check_imports.py --qxub  # Test all qxub packages
        """,
    )

    parser.add_argument(
        "packages", nargs="*", help="Package names to test (e.g., qxub.execution)"
    )

    parser.add_argument(
        "--all", action="store_true", help="Test all discoverable qxub packages"
    )

    parser.add_argument(
        "--qxub", action="store_true", help="Test all qxub sub-packages"
    )

    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Show verbose output"
    )

    args = parser.parse_args()

    # Determine packages to test
    if args.all or args.qxub:
        packages = get_qxub_packages()
    elif args.packages:
        packages = args.packages
    else:
        parser.print_help()
        sys.exit(1)

    print(f"🧪 Testing {len(packages)} package(s)...")
    print()

    # Test each package
    success_count = 0
    for package in packages:
        if test_import(package):
            success_count += 1

    print()
    print(f"📊 Results: {success_count}/{len(packages)} packages imported successfully")

    if success_count == len(packages):
        print("🎉 All imports successful!")
        sys.exit(0)
    else:
        print(f"⚠️  {len(packages) - success_count} package(s) had import issues")
        sys.exit(1)


if __name__ == "__main__":
    main()
